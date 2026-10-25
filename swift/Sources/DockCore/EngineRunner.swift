import Foundation

/// Runs the Vina CLI as a child process so cancel is a clean terminate (plan §4.2).
public final class EngineRunner {
    private var process: Process?

    public init() {}

    public func run(executable: URL, arguments: [String],
                    onLine: @escaping (String) -> Void) async throws -> Int32 {
        let p = Process()
        p.executableURL = executable
        p.arguments = arguments
        let pipe = Pipe()
        p.standardOutput = pipe
        p.standardError = pipe
        pipe.fileHandleForReading.readabilityHandler = { h in
            if let s = String(data: h.availableData, encoding: .utf8), !s.isEmpty {
                s.split(whereSeparator: \.isNewline).forEach { onLine(String($0)) }
            }
        }
        process = p
        return try await withCheckedThrowingContinuation { cont in
            p.terminationHandler = { proc in
                pipe.fileHandleForReading.readabilityHandler = nil
                cont.resume(returning: proc.terminationStatus)
            }
            do { try p.run() } catch { cont.resume(throwing: error) }
        }
    }

    /// Target: cancel within 2 s (N-05).
    public func cancel() { process?.terminate() }
}
