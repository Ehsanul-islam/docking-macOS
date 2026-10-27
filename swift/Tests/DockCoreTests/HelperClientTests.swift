import XCTest
@testable import DockCore

/// Talks to the real Python helper. Skipped (with the reason) when `python3 -m dockchem.server` is not available,
/// e.g. on a machine where the chemistry environment is not installed.
final class HelperClientTests: XCTestCase {
    func makeClient() throws -> HelperClient {
        let probe = Process()
        probe.executableURL = URL(fileURLWithPath: "/usr/bin/env")
        probe.arguments = ["python3", "-c", "import dockchem.server"]
        probe.standardError = Pipe(); probe.standardOutput = Pipe()
        do { try probe.run() } catch { throw XCTSkip("python3 not found") }
        probe.waitUntilExit()
        if probe.terminationStatus != 0 { throw XCTSkip("the dockchem package is not installed for python3") }
        return HelperClient(launch: .development())
    }

    func testPingReturnsAnEnvelope() async throws {
        let c = try makeClient(); defer { c.stop() }
        let r = try await c.call("system.ping")
        XCTAssertFalse(r.result.isNull)
        XCTAssertNotNil(r.provenance["cores"]?.int)
    }

    func testUnknownMethodIsAnRPCError() async throws {
        let c = try makeClient(); defer { c.stop() }
        do { _ = try await c.call("no.such.method"); XCTFail("expected an error") }
        catch let e as RPCError { XCTAssertEqual(e.code, -32601) }
    }

    func testTypedErrorCarriesKindAndMessage() async throws {
        let c = try makeClient(); defer { c.stop() }
        do { _ = try await c.call("structure.assess", ["path": "/definitely/not/here.pdb"]); XCTFail("expected an error") }
        catch let e as RPCError {
            XCTAssertEqual(e.kind, .invalidInput)
            XCTAssertFalse(e.message.isEmpty)
        }
    }

    func testProgressArrivesBeforeTheResult() async throws {
        let c = try makeClient(); defer { c.stop() }
        final class Seen: @unchecked Sendable { var f: [Double] = []; let l = NSLock()
            func add(_ x: Double) { l.lock(); f.append(x); l.unlock() } }
        let seen = Seen()
        let r = try await c.call("test.sleep", ["steps": 4, "interval": 0.02]) { seen.add($0.fraction) }
        XCTAssertEqual(r.result.string, "done")
        XCTAssertEqual(seen.f.last, 1.0)
        XCTAssertEqual(seen.f.count, 4)
    }

    func testCancellingTheTaskCancelsTheJob() async throws {
        let c = try makeClient(); defer { c.stop() }
        let start = Date()
        let t = Task { try await c.call("test.sleep", ["steps": 1000, "interval": 0.05]) }
        try await Task.sleep(nanoseconds: 500_000_000)
        t.cancel()
        do { _ = try await t.value; XCTFail("expected cancellation") }
        catch let e as RPCError { XCTAssertTrue(e.isCancelled) }
        XCTAssertLessThan(Date().timeIntervalSince(start), 3.0)
        let again = try await c.call("system.ping")          // the helper is still usable afterwards
        XCTAssertFalse(again.result.isNull)
    }

    func testSequentialCallsAndAHelperThatDiesIsReplaced() async throws {
        let c = try makeClient(); defer { c.stop() }
        _ = try await c.call("system.ping")
        c.stop()
        try await Task.sleep(nanoseconds: 300_000_000)
        XCTAssertFalse(c.isRunning)
        let r = try await c.call("system.ping")              // next call starts a fresh helper
        XCTAssertFalse(r.result.isNull)
    }

    func testCrashingHelperFailsPendingCallsWithTheTail() async throws {
        let c = HelperClient(launch: HelperClient.Launch(executable: URL(fileURLWithPath: "/bin/sh"),
                             arguments: ["-c", "echo boom >&2; exit 7"]))
        defer { c.stop() }
        do { _ = try await c.call("system.ping"); XCTFail("expected an error") }
        catch HelperClientError.helperExited(let status, let tail) {
            XCTAssertEqual(status, 7)
            XCTAssertTrue(tail.contains("boom"))
        } catch HelperClientError.notRunning { /* the write raced the exit; also a failure the app can show */ }
    }

    func testRPCErrorKeepsUnknownKindsAndExtraFields() throws {
        let json = #"{"code":-32007,"message":"x","data":{"type":"UnusableStructure","log":"","assessment":{"verdict":"unusable"}}}"#
        let e = try JSONDecoder().decode(RPCError.self, from: Data(json.utf8))
        XCTAssertEqual(e.kind, .unusableStructure)
        XCTAssertEqual(e.data?.extra["assessment"]?["verdict"]?.string, "unusable")
        let future = #"{"code":-1,"message":"x","data":{"type":"SomethingNew"}}"#
        XCTAssertEqual(try JSONDecoder().decode(RPCError.self, from: Data(future.utf8)).kind?.rawValue, "SomethingNew")
    }
}
