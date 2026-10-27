import XCTest
@testable import DockCore

final class DockCoreTests: XCTestCase {
    func testManifestDecodes() throws {
        let json = #"{"schema_version":1,"id":"x","name":"n","inputs":[{"path":"inputs/a","sha256":"00"}],"parameters":{},"seeds":{}}"#
        let m = try JSONDecoder().decode(ProjectManifest.self, from: Data(json.utf8))
        XCTAssertEqual(m.schemaVersion, 1)
        XCTAssertEqual(m.inputs.first?.path, "inputs/a")
    }

    func testRPCErrorDecodes() throws {
        let json = #"{"id":1,"error":{"code":-32001,"message":"bad","data":{"type":"InvalidInput","log":""}}}"#
        let r = try JSONDecoder().decode(RPCResponse.self, from: Data(json.utf8))
        XCTAssertEqual(r.error?.data?.type, .invalidInput)
    }

    func testCancelTerminatesProcessQuickly() async throws {
        let runner = EngineRunner()
        let start = Date()
        let task = Task { try await runner.run(executable: URL(fileURLWithPath: "/bin/sleep"), arguments: ["30"]) { _ in } }
        try await Task.sleep(nanoseconds: 300_000_000)
        runner.cancel()
        do { _ = try await task.value; XCTFail("expected cancellation") } catch { XCTAssertEqual(error as? EngineRunnerError, .cancelled) }
        XCTAssertLessThan(Date().timeIntervalSince(start), 2.0)
    }

    func testCancellingTheTaskStopsTheProcess() async throws {
        let runner = EngineRunner()
        let start = Date()
        let task = Task { try await runner.run(executable: URL(fileURLWithPath: "/bin/sleep"), arguments: ["30"]) { _ in } }
        try await Task.sleep(nanoseconds: 300_000_000)
        task.cancel()
        do { _ = try await task.value; XCTFail("expected cancellation") } catch { XCTAssertEqual(error as? EngineRunnerError, .cancelled) }
        XCTAssertLessThan(Date().timeIntervalSince(start), 2.0)
    }

    func testCancelBeforeStartNeverLaunches() async throws {
        let runner = EngineRunner()
        runner.cancel()          // nothing running yet: harmless
        let status = try await runner.run(executable: URL(fileURLWithPath: "/bin/echo"), arguments: ["hi"]) { _ in }
        XCTAssertEqual(status, 0)   // a stale cancel from before the run must not poison the next run
    }

    func testOutputArrivesAsWholeLinesAndExitStatusIsReturned() async throws {
        final class Box: @unchecked Sendable { var lines: [String] = []; let l = NSLock()
            func add(_ s: String) { l.lock(); lines.append(s); l.unlock() } }
        let box = Box()
        let status = try await EngineRunner().run(executable: URL(fileURLWithPath: "/bin/sh"),
            arguments: ["-c", "printf 'one\ntwo\rthree'; exit 3"]) { box.add($0) }
        XCTAssertEqual(status, 3)
        XCTAssertEqual(box.lines, ["one", "two", "three"])
    }

    func testMissingExecutableThrows() async {
        do { _ = try await EngineRunner().run(executable: URL(fileURLWithPath: "/nonexistent/vina"), arguments: []) { _ in }
             XCTFail("expected an error") } catch { }
    }
}
