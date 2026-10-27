import Foundation

public enum EngineRunnerError: Error, Equatable {
    case cancelled
    case alreadyRunning
}

/// Runs the Vina CLI (or any helper executable) as a child process so cancel is a clean terminate (plan §4.2).
///
/// One runner runs one process at a time. `cancel()` is safe from any thread, may be called before the process has
/// started (the run then never starts), sends SIGTERM and follows with SIGKILL if the process is still alive after
/// `killGrace` seconds.
public final class EngineRunner: @unchecked Sendable {
    private let lock = NSLock()
    private var process: Process?
    private var cancelled = false
    private var running = false
    public var killGrace: TimeInterval

    public init(killGrace: TimeInterval = 1.0) { self.killGrace = killGrace }

    /// Returns the exit status. Throws `EngineRunnerError.cancelled` if the run was cancelled (explicitly or by
    /// cancelling the surrounding Task) and `EngineRunnerError.alreadyRunning` on overlap.
    public func run(executable: URL, arguments: [String],
                    onLine: @escaping @Sendable (String) -> Void) async throws -> Int32 {
        try begin()
        defer { finish() }
        return try await withTaskCancellationHandler {
            try await launch(executable: executable, arguments: arguments, onLine: onLine)
        } onCancel: {
            self.cancel()
        }
    }

    /// Target: cancel within 2 s (N-05).
    public func cancel() {
        lock.lock()
        cancelled = true
        let p = process
        lock.unlock()
        guard let p, p.isRunning else { return }
        p.terminate()
        let pid = p.processIdentifier
        DispatchQueue.global().asyncAfter(deadline: .now() + killGrace) { [weak p] in
            if let p, p.isRunning { kill(pid, SIGKILL) }
        }
    }

    private func begin() throws {
        lock.lock(); defer { lock.unlock() }
        if running { throw EngineRunnerError.alreadyRunning }
        running = true
        cancelled = Task.isCancelled
    }

    private func wasCancelled() -> Bool {
        lock.lock(); defer { lock.unlock() }
        return cancelled
    }

    private func finish() {
        lock.lock(); defer { lock.unlock() }
        running = false
        process = nil
    }

    private func launch(executable: URL, arguments: [String],
                        onLine: @escaping @Sendable (String) -> Void) async throws -> Int32 {
        let p = Process()
        p.executableURL = executable
        p.arguments = arguments
        let pipe = Pipe()
        p.standardOutput = pipe
        p.standardError = pipe
        let splitter = LineSplitter(onLine)
        pipe.fileHandleForReading.readabilityHandler = { h in
            let data = h.availableData
            if !data.isEmpty { splitter.feed(data) }
        }
        let status: Int32 = try await withCheckedThrowingContinuation { cont in
            p.terminationHandler = { proc in
                pipe.fileHandleForReading.readabilityHandler = nil
                splitter.feed(pipe.fileHandleForReading.availableData)
                splitter.flush()
                cont.resume(returning: proc.terminationStatus)
            }
            lock.lock()
            if cancelled {
                lock.unlock()
                cont.resume(throwing: EngineRunnerError.cancelled)
                return
            }
            do {
                try p.run()
                process = p
                lock.unlock()
            } catch {
                lock.unlock()
                pipe.fileHandleForReading.readabilityHandler = nil
                cont.resume(throwing: error)
            }
        }
        if wasCancelled() { throw EngineRunnerError.cancelled }
        return status
    }
}

/// Turns arbitrary chunks of output into whole lines (a chunk may end mid-line, and Vina's progress bar uses `\r`).
final class LineSplitter: @unchecked Sendable {
    private var buffer = Data()
    private let lock = NSLock()
    private let emit: @Sendable (String) -> Void

    init(_ emit: @escaping @Sendable (String) -> Void) { self.emit = emit }

    func feed(_ data: Data) {
        guard !data.isEmpty else { return }
        lock.lock()
        buffer.append(data)
        var lines: [String] = []
        while let i = buffer.firstIndex(where: { $0 == 0x0A || $0 == 0x0D }) {
            let line = buffer[buffer.startIndex..<i]
            buffer.removeSubrange(buffer.startIndex...i)
            if !line.isEmpty { lines.append(String(decoding: line, as: UTF8.self)) }
        }
        lock.unlock()
        lines.forEach(emit)
    }

    func flush() {
        lock.lock()
        let rest = buffer
        buffer.removeAll()
        lock.unlock()
        if !rest.isEmpty { emit(String(decoding: rest, as: UTF8.self)) }
    }
}
