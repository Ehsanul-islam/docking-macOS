import XCTest
@testable import DockCore

/// End to end: the real Python helper and the real Vina binary. Skipped, with the reason, when either is missing
/// (set DOCKCHEM_VINA to a Vina binary to enable the docking tests).
@MainActor
final class QuickDockSessionTests: XCTestCase {
    var work: URL!
    var client: HelperClient!

    override func setUp() async throws {
        try await super.setUp()
        let probe = Process()
        probe.executableURL = URL(fileURLWithPath: "/usr/bin/env")
        probe.arguments = ["python3", "-c", "import dockchem.server"]
        probe.standardError = Pipe(); probe.standardOutput = Pipe()
        do { try probe.run() } catch { throw XCTSkip("python3 not found") }
        probe.waitUntilExit()
        if probe.terminationStatus != 0 { throw XCTSkip("the dockchem package is not installed for python3") }
        work = FileManager.default.temporaryDirectory.appendingPathComponent("dock-\(UUID().uuidString)")
        try FileManager.default.createDirectory(at: work, withIntermediateDirectories: true)
        client = HelperClient(launch: .development())
    }

    override func tearDown() async throws {
        client?.stop()
        if let work { try? FileManager.default.removeItem(at: work) }
        try await super.tearDown()
    }

    var receptor: URL {
        URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("pdb_structures_50/1STP.pdb")
    }

    func requireEngine() throws {
        let env = ProcessInfo.processInfo.environment
        guard let p = env["DOCKCHEM_VINA"], FileManager.default.isExecutableFile(atPath: p) else {
            throw XCTSkip("set DOCKCHEM_VINA to a Vina binary to run the docking tests")
        }
    }

    func session() -> QuickDockSession {
        let s = QuickDockSession(helper: client, workDirectory: work)
        s.receptor = receptor
        s.ligand = "CC1C2C(CCCCC(=O)O)SC1C(=O)N2"
        s.preset = .fast
        return s
    }

    func testAssessShowsAReadyStructure() async throws {
        let s = session()
        s.assess()
        await s.waitUntilIdle()
        guard case .assessed(let a) = s.phase else { return XCTFail("phase \(s.phase)") }
        XCTAssertEqual(a.verdict, .ready)
    }

    func testEmptyLigandAndMissingFileEndInReadableFailures() async throws {
        let s = session()
        s.ligand = "  "
        s.prepare()
        await s.waitUntilIdle()
        guard case .failed(let e1) = s.phase else { return XCTFail("phase \(s.phase)") }
        XCTAssertEqual(e1.title, "No ligand")

        let s2 = session()
        s2.receptor = work.appendingPathComponent("missing.pdb")
        s2.assess()
        await s2.waitUntilIdle()
        guard case .failed(let e2) = s2.phase else { return XCTFail("phase \(s2.phase)") }
        XCTAssertEqual(e2.title, "Check the input")
    }

    func testFullRunProducesRankedResultsAndSessionStaysUsable() async throws {
        try requireEngine()
        let s = session()
        s.runAll()
        await s.waitUntilIdle()
        guard case .finished(let d) = s.phase else { return XCTFail("phase \(s.phase)") }
        XCTAssertEqual(d.docked, 1)
        let best = try XCTUnwrap(d.ranked.first?.affinity)
        XCTAssertLessThan(best, 0)
        XCTAssertTrue(FileManager.default.fileExists(atPath: d.csvPath))
        XCTAssertFalse(s.guardrails.isEmpty)                         // blind docking at exhaustiveness 8
        XCTAssertTrue(s.guardrails.contains { $0.id == "blind-undersampled" })
        s.reset()
        XCTAssertEqual(s.phase, .idle)
    }

    func testCancelWhileDockingEndsCancelledAndHelperStillWorks() async throws {
        try requireEngine()
        let s = session()
        s.preset = .thorough
        s.runAll()
        for _ in 0..<600 {                                            // wait until docking has started
            if case .docking = s.phase { break }
            try await Task.sleep(nanoseconds: 100_000_000)
        }
        guard case .docking = s.phase else { return XCTFail("never reached docking: \(s.phase)") }
        try await Task.sleep(nanoseconds: 1_000_000_000)
        let t0 = Date()
        s.cancel()
        await s.waitUntilIdle()
        XCTAssertEqual(s.phase, .cancelled)
        XCTAssertLessThan(Date().timeIntervalSince(t0), 15)
        let again = try await client.call("system.ping")
        XCTAssertFalse(again.result.isNull)
    }
}
