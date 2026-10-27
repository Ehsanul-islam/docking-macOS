import Foundation

public enum HelperClientError: Error, LocalizedError, Sendable {
    case notRunning
    case failedToStart(String)
    case helperExited(status: Int32, stderrTail: String)
    case malformedResponse(String)

    public var errorDescription: String? {
        switch self {
        case .notRunning: return "The chemistry helper is not running."
        case .failedToStart(let why): return "The chemistry helper could not be started: \(why)"
        case .helperExited(let status, let tail):
            return "The chemistry helper stopped unexpectedly (exit status \(status))." + (tail.isEmpty ? "" : "\n\(tail)")
        case .malformedResponse(let why): return "The chemistry helper sent an unreadable reply: \(why)"
        }
    }
}

public struct HelperProgress: Sendable, Equatable {
    public let fraction: Double
    public let message: String
}

/// JSON-RPC 2.0 client for the bundled Python chemistry helper (plan §4.4): one child process, requests are
/// answered in order by the helper, a `cancel` notification takes effect at once, and every failure is typed.
///
/// Cancelling the calling Task sends `cancel` for that job; the call then ends with an `RPCError` whose kind is
/// `.cancelled`. If the helper process dies, every pending call fails with `HelperClientError.helperExited`, and the
/// next call starts a fresh helper.
public final class HelperClient: @unchecked Sendable {
    public struct Launch: Sendable {
        public var executable: URL
        public var arguments: [String]
        public var environment: [String: String]
        public init(executable: URL, arguments: [String] = ["-m", "dockchem.server"], environment: [String: String] = [:]) {
            self.executable = executable; self.arguments = arguments; self.environment = environment
        }
        /// Development default: `python3 -m dockchem.server` found through /usr/bin/env.
        public static func development(python: String = "python3", environment: [String: String] = [:]) -> Launch {
            Launch(executable: URL(fileURLWithPath: "/usr/bin/env"), arguments: [python, "-m", "dockchem.server"],
                   environment: environment)
        }
    }

    private struct Pending {
        var cont: CheckedContinuation<HelperResult, Error>
        var onProgress: (@Sendable (HelperProgress) -> Void)?
    }

    private let launch: Launch
    private let lock = NSLock()
    private var process: Process?
    private var stdin: FileHandle?
    private var nextID = 1
    private var pending: [Int: Pending] = [:]
    private var stderrTail = Data()
    private static let stderrKeep = 4000

    public init(launch: Launch) { self.launch = launch }

    deinit { stop() }

    // MARK: lifecycle

    public func start() throws {
        lock.lock(); defer { lock.unlock() }
        try startLocked()
    }

    private func startLocked() throws {
        if let p = process, p.isRunning { return }
        let p = Process()
        p.executableURL = launch.executable
        p.arguments = launch.arguments
        var env = ProcessInfo.processInfo.environment
        launch.environment.forEach { env[$0.key] = $0.value }
        p.environment = env
        let inPipe = Pipe(), outPipe = Pipe(), errPipe = Pipe()
        p.standardInput = inPipe
        p.standardOutput = outPipe
        p.standardError = errPipe
        let splitter = LineSplitter { [weak self] line in self?.handleLine(line) }
        outPipe.fileHandleForReading.readabilityHandler = { h in
            let d = h.availableData
            if !d.isEmpty { splitter.feed(d) }
        }
        errPipe.fileHandleForReading.readabilityHandler = { [weak self] h in
            let d = h.availableData
            if !d.isEmpty { self?.appendStderr(d) }
        }
        p.terminationHandler = { [weak self] proc in
            outPipe.fileHandleForReading.readabilityHandler = nil
            errPipe.fileHandleForReading.readabilityHandler = nil
            splitter.feed(outPipe.fileHandleForReading.availableData)
            splitter.flush()
            self?.processEnded(proc)
        }
        do { try p.run() } catch { throw HelperClientError.failedToStart(error.localizedDescription) }
        process = p
        stdin = inPipe.fileHandleForWriting
    }

    /// Ends the helper. Pending calls fail with `helperExited`.
    public func stop() {
        lock.lock()
        let p = process
        let h = stdin
        lock.unlock()
        try? h?.close()
        if let p, p.isRunning { p.terminate() }
    }

    public var isRunning: Bool { lock.lock(); defer { lock.unlock() }; return process?.isRunning ?? false }

    // MARK: calls

    /// Calls one helper method. `onProgress` receives the helper's progress notifications for this call.
    public func call(_ method: String, _ params: [String: JSONValue] = [:],
                     onProgress: (@Sendable (HelperProgress) -> Void)? = nil) async throws -> HelperResult {
        let id = try reserveID()
        return try await withTaskCancellationHandler {
            try await withCheckedThrowingContinuation { (cont: CheckedContinuation<HelperResult, Error>) in
                var p = params
                p["job_id"] = .int(id)
                let req = RPCRequest(id: id, method: method, params: p)
                lock.lock()
                if Task.isCancelled {
                    lock.unlock()
                    cont.resume(throwing: RPCError.cancelledBeforeStart)
                    return
                }
                pending[id] = Pending(cont: cont, onProgress: onProgress)
                lock.unlock()
                do { try send(req) } catch {
                    if let pend = takePending(id) { pend.cont.resume(throwing: error) }
                }
            }
        } onCancel: {
            self.sendCancel(jobID: id)
        }
    }

    private func reserveID() throws -> Int {
        lock.lock(); defer { lock.unlock() }
        try startLocked()
        defer { nextID += 1 }
        return nextID
    }

    private func send<T: Encodable>(_ message: T) throws {
        let data = try JSONEncoder().encode(message) + Data([0x0A])
        lock.lock(); let h = stdin; lock.unlock()
        guard let h else { throw HelperClientError.notRunning }
        do { try h.write(contentsOf: data) } catch { throw HelperClientError.notRunning }
    }

    private func sendCancel(jobID: Int) {
        struct Cancel: Encodable { let jsonrpc = "2.0"; let method = "cancel"; let params: [String: Int] }
        try? send(Cancel(params: ["job_id": jobID]))
    }

    private func takePending(_ id: Int) -> Pending? {
        lock.lock(); defer { lock.unlock() }
        return pending.removeValue(forKey: id)
    }

    // MARK: incoming

    private func handleLine(_ line: String) {
        guard let data = line.data(using: .utf8), let obj = try? JSONDecoder().decode([String: JSONValue].self, from: data) else {
            return          // not protocol traffic (the helper keeps stdout for JSON only; ignore strays)
        }
        if obj["method"]?.string == "progress", let params = obj["params"] {
            guard let job = params["job_id"]?.int else { return }
            lock.lock(); let handler = pending[job]?.onProgress; lock.unlock()
            handler?(HelperProgress(fraction: params["fraction"]?.double ?? 0, message: params["message"]?.string ?? ""))
            return
        }
        guard let id = obj["id"]?.int, let pend = takePending(id) else { return }
        do {
            let resp = try JSONDecoder().decode(RPCResponse.self, from: data)
            if let err = resp.error { pend.cont.resume(throwing: err); return }
            guard let env = resp.result else {
                pend.cont.resume(throwing: HelperClientError.malformedResponse("reply has neither result nor error")); return
            }
            pend.cont.resume(returning: HelperResult(
                result: env["result"] ?? .null,
                warnings: env["warnings"]?.array?.compactMap(\.string) ?? [],
                provenance: env["provenance"] ?? .null,
                artifacts: env["artifacts"]?.array?.compactMap(\.string) ?? []))
        } catch {
            pend.cont.resume(throwing: HelperClientError.malformedResponse("\(error)"))
        }
    }

    private func appendStderr(_ d: Data) {
        lock.lock(); defer { lock.unlock() }
        stderrTail.append(d)
        if stderrTail.count > Self.stderrKeep { stderrTail = stderrTail.suffix(Self.stderrKeep) }
    }

    private func processEnded(_ proc: Process) {
        lock.lock()
        if process === proc { process = nil; stdin = nil }
        let tail = String(decoding: stderrTail, as: UTF8.self)
        let orphans = pending
        pending.removeAll()
        lock.unlock()
        for (_, p) in orphans {
            p.cont.resume(throwing: HelperClientError.helperExited(status: proc.terminationStatus, stderrTail: tail))
        }
    }
}

extension RPCError {
    /// A call whose Task was already cancelled never reaches the helper.
    static var cancelledBeforeStart: RPCError {
        let json = #"{"code":-32004,"message":"cancelled","data":{"type":"Cancelled","log":""}}"#
        return try! JSONDecoder().decode(RPCError.self, from: Data(json.utf8))
    }
}
