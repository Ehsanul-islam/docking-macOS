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
        async let status = runner.run(executable: URL(fileURLWithPath: "/bin/sleep"), arguments: ["30"]) { _ in }
        try await Task.sleep(nanoseconds: 200_000_000)
        runner.cancel()
        _ = try await status
        XCTAssertLessThan(Date().timeIntervalSince(start), 2.0)
    }
}
