import XCTest
@testable import DockCore

final class ModelTests: XCTestCase {
    func fixture(_ name: String) throws -> JSONValue {
        let url = try XCTUnwrap(Bundle.module.url(forResource: name, withExtension: "json", subdirectory: "Fixtures"))
        let env = try JSONDecoder().decode(JSONValue.self, from: Data(contentsOf: url))
        return try XCTUnwrap(env["result"]?["result"])
    }

    func testAssessmentFromARealHelperReply() throws {
        let a = try XCTUnwrap(Assessment(json: fixture("assess_ready")))
        XCTAssertEqual(a.verdict, .ready)
        XCTAssertEqual(a.proteinResidues, 121)
        XCTAssertEqual(a.resolution, 2.6)
        XCTAssertEqual(a.findings.first?.code, "ligands_present")
        XCTAssertTrue(a.attention.isEmpty)
    }

    func testQuickPrepFromARealHelperReply() throws {
        let q = try XCTUnwrap(QuickPrepResult(json: fixture("quickprep_ok")))
        XCTAssertTrue(q.readyToDock)
        XCTAssertTrue(q.isBlind)
        XCTAssertEqual(q.box.center.count, 3)
        XCTAssertEqual(q.engine.exhaustiveness, 8)
        XCTAssertEqual(q.ligandsPrepared, 1)
        XCTAssertEqual(q.ligandFiles.count, 1)
        XCTAssertEqual(q.warnings.count, 5)
        XCTAssertEqual(q.assessment?.verdict, .ready)
        XCTAssertGreaterThan(q.box.volume, Guardrails.largeBoxVolume)
    }

    func testDockResultRanksBestFirstAndPutsFailuresLast() throws {
        let d = try XCTUnwrap(DockRunResult(json: fixture("dock_ok")))
        XCTAssertEqual(d.docked, 1)
        XCTAssertEqual(d.results.first?.affinity, -4.575)
        XCTAssertEqual(d.results.first?.scores.count, 9)
        let mixed: JSONValue = ["results": [
            ["name": "bad", "status": "failed", "error": "x"],
            ["name": "weak", "status": "docked", "affinity": -3.0],
            ["name": "strong", "status": "docked", "affinity": -9.5]]]
        let r = try XCTUnwrap(DockRunResult(json: mixed))
        XCTAssertEqual(r.ranked.map(\.name), ["strong", "weak", "bad"])
        XCTAssertEqual(r.ranked.last?.status, .failed)
    }

    func testMissingFieldsGiveNilNotACrash() {
        XCTAssertNil(QuickPrepResult(json: [:]))
        XCTAssertNil(DockRunResult(json: ["x": 1]))
        XCTAssertNil(Assessment(json: [:]))
        XCTAssertEqual(Assessment(json: ["verdict": "from_the_future"])?.verdict, .unknown)
    }

    func testGuardrailsFollowThePlan() {
        let big = SearchBox(center: [0, 0, 0], size: [40, 40, 40])
        let small = SearchBox(center: [0, 0, 0], size: [20, 20, 20])
        XCTAssertEqual(Guardrails.check(box: big, exhaustiveness: 32, blind: false).map(\.id), ["large-box"])
        XCTAssertEqual(Guardrails.check(box: small, exhaustiveness: 8, blind: true).map(\.id), ["blind-undersampled"])
        XCTAssertTrue(Guardrails.check(box: small, exhaustiveness: 32, blind: true).isEmpty)
        let many = Guardrails.check(box: small, exhaustiveness: 32, blind: false, maxTorsions: 25, ligandsFailed: 2, receptorLeftOut: 1)
        XCTAssertEqual(Set(many.map(\.id)), ["many-torsions", "ligands-failed", "receptor-left-out"])
        XCTAssertEqual(many.last?.id, "receptor-left-out")        // information sorts after cautions
    }

    func testSettingsClampCores() {
        XCTAssertEqual(AppSettings().cores, 2)
        XCTAssertEqual(AppSettings(cores: 0).effectiveCores(available: 8), 1)
        XCTAssertEqual(AppSettings(cores: 64).effectiveCores(available: 8), 8)
        XCTAssertEqual(AppSettings(cores: 4).effectiveCores(available: 8), 4)
        let data = try! JSONEncoder().encode(AppSettings(cores: 3, enginePath: "/x/vina"))
        XCTAssertEqual(try! JSONDecoder().decode(AppSettings.self, from: data), AppSettings(cores: 3, enginePath: "/x/vina"))
    }

    func rpc(_ type: String, message: String = "m", extra: String = "") throws -> RPCError {
        let json = #"{"code":-32001,"message":"\#(message)","data":{"type":"\#(type)","log":"L"\#(extra)}}"#
        return try JSONDecoder().decode(RPCError.self, from: Data(json.utf8))
    }

    func testUserFacingErrorsTellTheUserWhatToDo() throws {
        let unusable = UserFacingError.from(try rpc("UnusableStructure", message: "only CA atoms"))
        XCTAssertEqual(unusable.action, .prepareAnyway)
        XCTAssertTrue(unusable.detail.contains("only CA atoms"))
        let tpl = UserFacingError.from(try rpc("NeedsTemplate", extra: #","resname":"XYZ""#))
        XCTAssertEqual(tpl.action, .provideBondOrders("XYZ"))
        XCTAssertTrue(UserFacingError.from(try rpc("Cancelled")).isCancellation)
        XCTAssertEqual(UserFacingError.from(try rpc("UnsupportedFeature", message: "the docking engine (AutoDock Vina) was not found")).action, .chooseEngine)
        XCTAssertEqual(UserFacingError.from(try rpc("InternalError")).action, .copyDetails)
        XCTAssertEqual(UserFacingError.from(try rpc("SomethingNew")).action, .copyDetails)
        XCTAssertEqual(UserFacingError.from(try rpc("InvalidInput", message: "bad file")).detail, "bad file")
        XCTAssertTrue(UserFacingError.from(CancellationError()).isCancellation)
        let died = UserFacingError.from(HelperClientError.helperExited(status: 9, stderrTail: "oom"))
        XCTAssertEqual(died.log, "oom")
    }
}

final class HelperLocatorTests: XCTestCase {
    func testBundledPythonAndEngineWinOverEverythingElse() {
        let res = URL(fileURLWithPath: "/App/Contents/Resources")
        let l = HelperLocator(resourceDirectory: res, environment: ["DOCKCHEM_PYTHON": "/dev/py", "DOCKCHEM_VINA": "/dev/vina"],
                              fileExists: { _ in true })
        let launch = l.launch()
        XCTAssertEqual(launch.executable.path, "/App/Contents/Resources/python/bin/python3")
        XCTAssertEqual(launch.environment["DOCKCHEM_VINA"], "/App/Contents/Resources/engine/vina")
        XCTAssertEqual(launch.environment["PYTHONHOME"], "/App/Contents/Resources/python")
        XCTAssertEqual(launch.arguments, ["-m", "dockchem.server"])
    }

    func testFallsBackToEnvironmentThenToSystemPython() {
        let env = HelperLocator(resourceDirectory: nil, environment: ["DOCKCHEM_PYTHON": "/opt/py", "DOCKCHEM_VINA": "/opt/vina"],
                                fileExists: { $0.hasPrefix("/opt") })
        XCTAssertEqual(env.launch().executable.path, "/opt/py")
        XCTAssertEqual(env.enginePath(), "/opt/vina")
        let none = HelperLocator(resourceDirectory: nil, environment: [:], fileExists: { _ in false })
        XCTAssertEqual(none.launch().executable.path, "/usr/bin/env")
        XCTAssertNil(none.enginePath())
    }
}
