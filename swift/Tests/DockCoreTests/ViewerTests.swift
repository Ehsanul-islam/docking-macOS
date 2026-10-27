import XCTest
@testable import DockCore

@MainActor
final class MockTransport: ViewerTransport {
    var sent: [ViewerCommand] = []
    func send(_ command: ViewerCommand) { sent.append(command) }
    func decoded(_ i: Int) throws -> JSONValue { try JSONDecoder().decode(JSONValue.self, from: Data(sent[i].json.utf8)) }
}

@MainActor
final class ViewerTests: XCTestCase {
    var dataDir: URL {
        URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("viewer/test/data")
    }

    func testCommandsAreValidJSONWithTheNamesThePageExpects() async throws {
        let c = ViewerCommand.loadPose(sdf: Data("abc".utf8), index: 3, score: -7.5, poses: [(1, -8.0), (2, nil)])
        let j = try JSONDecoder().decode(JSONValue.self, from: Data(c.json.utf8))
        XCTAssertEqual(j["cmd"]?.string, "loadPose")
        XCTAssertEqual(j["index"]?.int, 3)
        XCTAssertEqual(j["format"]?.string, "sdf")
        XCTAssertEqual(Data(base64Encoded: j["data"]?.string ?? ""), Data("abc".utf8))
        XCTAssertEqual(j["poses"]?[1]?["score"]?.isNull, true)
        let box = try JSONDecoder().decode(JSONValue.self, from: Data(ViewerCommand.setBox(SearchBox(center: [1, 2, 3], size: [20, 20, 20])).json.utf8))
        XCTAssertEqual(box["center"]?[2]?.double, 3)
        XCTAssertEqual(try JSONDecoder().decode(JSONValue.self, from: Data(ViewerCommand.setBox(nil).json.utf8))["center"]?.isNull, true)
        XCTAssertEqual(try JSONDecoder().decode(JSONValue.self, from: Data(ViewerCommand.setOptions(["surface": true]).json.utf8))["surface"]?.bool, true)
        XCTAssertEqual(ViewerCommand.clear.name, "clear")
    }

    func testEventsFromThePage() async {
        XCTAssertEqual(ViewerEvent(["type": "ready"]), .ready)
        XCTAssertEqual(ViewerEvent(["type": "poseSelected", "index": 4]), .poseSelected(4))
        XCTAssertNil(ViewerEvent(["type": "poseSelected"]))
        XCTAssertNil(ViewerEvent(["type": "nonsense"]))
        XCTAssertEqual(ViewerEvent(["type": "error", "message": "boom"]), .error("boom"))
        let png = Data([0x89, 0x50, 0x4E, 0x47])
        XCTAssertEqual(ViewerEvent(["type": "screenshot", "dataUri": .string("data:image/png;base64," + png.base64EncodedString())]), .screenshot(png))
        XCTAssertNil(ViewerEvent(["type": "screenshot", "dataUri": "garbage"]))
    }

    func testCommandsWaitForReadyAndKeepTheirOrder() async throws {
        let t = MockTransport()
        let v = ViewerController(helper: nil)
        v.attach(t)
        v.send(.setBox(nil)); v.send(.focusLigand())
        XCTAssertTrue(t.sent.isEmpty)
        v.handle(.ready)
        XCTAssertEqual(t.sent.map(\.name), ["setTheme", "setBox", "focus"])
        v.send(.clear)
        XCTAssertEqual(t.sent.last?.name, "clear")
        v.theme = .light
        XCTAssertEqual(try t.decoded(t.sent.count - 1)["theme"]?.string, "light")
    }

    func testAReloadedPageGetsTheContentAgain() async throws {
        let t = MockTransport()
        let v = ViewerController(helper: nil)
        v.attach(t); v.handle(.ready)
        v.showStructure(at: dataDir.appendingPathComponent("receptor.pdb"))
        XCTAssertEqual(t.sent.last?.name, "loadReceptor")
        v.pageReloaded()
        XCTAssertFalse(v.isReady)
        v.showStructure(at: dataDir.appendingPathComponent("receptor.pdb"))     // queued until the page is ready again
        let before = t.sent.count
        v.handle(.ready)
        XCTAssertGreaterThan(t.sent.count, before)
        v.showStructure(at: dataDir.appendingPathComponent("missing.pdb"))
        XCTAssertNotNil(v.lastError)
    }

    func testShowingAResultSendsReceptorBoxPoseAndRealInteractions() async throws {
        let probe = Process()
        probe.executableURL = URL(fileURLWithPath: "/usr/bin/env"); probe.arguments = ["python3", "-c", "import dockchem.server"]
        probe.standardError = Pipe(); probe.standardOutput = Pipe()
        do { try probe.run() } catch { throw XCTSkip("python3 not found") }
        probe.waitUntilExit()
        if probe.terminationStatus != 0 { throw XCTSkip("the dockchem package is not installed for python3") }

        let helper = HelperClient(launch: .development()); defer { helper.stop() }
        let prepared = try XCTUnwrap(QuickPrepResult(json: [
            "site": ["how": "x", "center": [9.7, 2.7, 0.4], "size": [40.0, 30.0, 45.0]],
            "receptor": ["status": "ok", "pdbqt": "", "pdb": .string(dataDir.appendingPathComponent("receptor.pdb").path)],
            "engine": ["scoring": "vina"], "ready_to_dock": true]))
        let ligand = DockedLigand(json: ["name": "biotin", "status": "docked", "affinity": -4.575,
                                         "scores": [-4.575, -4.564, -4.495], "sdf": .string(dataDir.appendingPathComponent("poses.sdf").path)])
        let t = MockTransport()
        let v = ViewerController(helper: helper)
        v.attach(t); v.handle(.ready)
        await v.showResult(prepared: prepared, ligand: ligand, pose: 2, structureName: "1STP")
        XCTAssertEqual(t.sent.map(\.name), ["setTheme", "loadReceptor", "setBox", "loadPose", "setInteractions", "setInteractions"])
        XCTAssertEqual(try t.decoded(3)["index"]?.int, 2)
        XCTAssertEqual(try t.decoded(3)["score"]?.double, -4.564)
        let ints = try t.decoded(5)["interactions"]?.array ?? []
        XCTAssertFalse(ints.isEmpty)
        XCTAssertNotNil(ints[0]["lig"]?.array?.count == 3 ? true : nil)
        // switching pose keeps the receptor: only pose + interactions are sent
        let n = t.sent.count
        await v.showPose(1, of: ligand, prepared: prepared)
        XCTAssertEqual(t.sent.dropFirst(n).map(\.name), ["loadPose", "setInteractions", "setInteractions"])
        // a pose selected on the page
        let before = t.sent.count
        v.handle(.poseSelected(3))
        for _ in 0..<200 where t.sent.count < before + 3 { try await Task.sleep(nanoseconds: 100_000_000) }   // up to 20 s
        XCTAssertEqual(t.sent.count, before + 3)
        XCTAssertEqual(v.poseIndex, 3)
        XCTAssertNil(v.lastError)
    }
}
