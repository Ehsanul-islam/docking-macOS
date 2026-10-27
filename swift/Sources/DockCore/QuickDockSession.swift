import Foundation
import Observation

/// Where to dock (quick prep, plan §9.2 step 3). `auto` uses the largest drug-like ligand in the structure, or
/// blind docking when there is none; the result says which.
public enum SiteChoice: Sendable, Equatable {
    case auto
    case ligand(String)                          // "A:OHT:600" or "OHT"
    case box(center: [Double], size: [Double])
    case residues([String])                      // "A:394"

    var json: JSONValue? {
        switch self {
        case .auto: return nil
        case .ligand(let l): return ["ligand": .string(l)]
        case .box(let c, let s): return ["center": .array(c.map(JSONValue.double)), "size": .array(s.map(JSONValue.double))]
        case .residues(let r): return ["residues": .array(r.map(JSONValue.string))]
        }
    }
}

public enum DockPreset: String, Sendable, CaseIterable { case fast, standard, thorough }

/// The one-panel flow: structure in, ligand(s) in, ranked results out. UI-free and `@Observable`, so a SwiftUI view
/// can bind to it and a test can drive it. Every step is cancellable and every failure ends in `.failed` with a
/// message the user can act on.
@MainActor @Observable
public final class QuickDockSession {
    public enum Phase: Equatable {
        case idle
        case assessing
        case assessed(Assessment)
        case preparing
        case prepared(QuickPrepResult)
        case docking(fraction: Double, message: String)
        case finished(DockRunResult)
        case failed(UserFacingError)
        case cancelled
    }

    public private(set) var phase: Phase = .idle
    public private(set) var assessment: Assessment?
    public private(set) var prepared: QuickPrepResult?
    public private(set) var result: DockRunResult?
    public private(set) var warnings: [String] = []
    public private(set) var guardrails: [Guardrail] = []

    public var receptor: URL?
    /// A SMILES string, or the path of an SDF / MOL / SMI / PDB file (several molecules = a library).
    public var ligand: String = ""
    public var site: SiteChoice = .auto
    public var preset: DockPreset = .standard
    public var allowUnusable = false
    public var settings: AppSettings

    public let workDirectory: URL
    private let helper: HelperClient
    private var task: Task<Void, Never>?

    public init(helper: HelperClient, workDirectory: URL, settings: AppSettings = AppSettings()) {
        self.helper = helper; self.workDirectory = workDirectory; self.settings = settings
    }

    public var isBusy: Bool {
        switch phase { case .assessing, .preparing, .docking: return true; default: return false }
    }

    private var prepDir: URL { workDirectory.appendingPathComponent("prep") }
    private var runDir: URL { workDirectory.appendingPathComponent("run") }

    private func common() -> [String: JSONValue] {
        var p: [String: JSONValue] = ["cores": .int(settings.effectiveCores())]
        if let e = settings.enginePath { p["engine"] = .string(e) }
        return p
    }

    // MARK: steps

    /// Step 1: look at the structure. Cheap (about a second) and safe to call as soon as a file is chosen.
    public func assess() {
        guard receptor != nil, !isBusy else { return }
        start(.assessing) { [self] in
            let a = try await assessStructure()
            phase = .assessed(a)
        }
    }

    /// Step 2: prepare receptor and ligand(s) and choose the box.
    public func prepare() {
        guard receptor != nil, !isBusy else { return }
        start(.preparing) { [self] in
            phase = .prepared(try await prepareInputs())
        }
    }

    /// Step 3: dock everything that was prepared.
    public func dock() {
        guard prepared != nil, !isBusy else { return }
        start(.docking(fraction: 0, message: "Starting")) { [self] in
            phase = .finished(try await runDocking())
        }
    }

    /// All three steps; stops at the first problem. A structure the assessment calls unusable is refused by the
    /// preparation step unless `allowUnusable` was set by the user.
    public func runAll() {
        guard receptor != nil, !isBusy else { return }
        start(.assessing) { [self] in
            _ = try await assessStructure()
            phase = .preparing
            _ = try await prepareInputs()
            phase = .docking(fraction: 0, message: "Starting")
            phase = .finished(try await runDocking())
        }
    }

    private func assessStructure() async throws -> Assessment {
        let r = try await helper.call("structure.assess", ["path": .string(receptor!.path)])
        let a = try r.decode(Assessment.init(json:), as: "assessment")
        assessment = a
        warnings = r.warnings
        return a
    }

    private func prepareInputs() async throws -> QuickPrepResult {
        let smiles = ligand.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !smiles.isEmpty else {
            throw UserFacingError(title: "No ligand", detail: "Enter a SMILES string or choose a ligand file.")
        }
        var p = common()
        p["receptor"] = .string(receptor!.path)
        p["ligand"] = .string(smiles)
        p["outdir"] = .string(prepDir.path)
        p["preset"] = .string(preset.rawValue)
        p["allow_unusable"] = .bool(allowUnusable)
        if let s = site.json { p["site"] = s }
        let q = try await helper.call("quickprep.run", p).decode(QuickPrepResult.init(json:), as: "preparation")
        prepared = q
        warnings = q.warnings
        guardrails = Guardrails.check(box: q.box, exhaustiveness: q.engine.exhaustiveness, blind: q.isBlind,
                                      ligandsFailed: q.ligandsFailed.count, receptorLeftOut: q.leftOutResidues.count)
        return q
    }

    private func runDocking() async throws -> DockRunResult {
        var p = common()
        p["quickprep_dir"] = .string(prepDir.path)
        p["outdir"] = .string(runDir.path)
        let r = try await helper.call("dock.run", p) { [weak self] prog in
            Task { @MainActor in
                guard let self, case .docking = self.phase else { return }
                self.phase = .docking(fraction: prog.fraction, message: prog.message)
            }
        }
        let d = try r.decode(DockRunResult.init(json:), as: "docking result")
        result = d
        return d
    }

    /// 2D interaction diagram for one pose of one docked ligand (poses are numbered from 1). Returns the PNG path;
    /// the same call is cached on disk, so asking again for the same pose and style is instant.
    public func interactionImage(for ligand: DockedLigand, pose: Int = 1, style: String = "reference") async throws -> URL {
        guard let prepared else { throw UserFacingError(title: "Nothing to draw", detail: "Prepare and dock first.") }
        let safe = ligand.name.map { $0.isLetter || $0.isNumber || $0 == "-" ? $0 : "_" }
        let out = workDirectory.appendingPathComponent("images/\(String(safe))_pose\(pose)_\(style).png")
        if FileManager.default.fileExists(atPath: out.path) { return out }
        var p: [String: JSONValue] = ["receptor": .string(prepared.receptorPDB), "pose": .string(ligand.posesPDBQT),
                                      "pose_index": .int(pose), "out": .string(out.path), "style": .string(style),
                                      "compound": .string(ligand.name), "dpi": .int(200)]
        if pose >= 1, pose <= ligand.scores.count { p["affinity"] = .double(ligand.scores[pose - 1]) }
        _ = try await helper.call("complex.image2d", p)
        return out
    }

    /// Stops whatever is running. The helper stops its tools within a few seconds; partial results are kept
    /// on disk and a later run resumes from them.
    public func cancel() { task?.cancel() }

    public func reset() {
        task?.cancel()
        phase = .idle; assessment = nil; prepared = nil; result = nil; warnings = []; guardrails = []
    }

    // MARK: plumbing

    private func start(_ first: Phase, _ body: @escaping @MainActor () async throws -> Void) {
        phase = first
        task = Task { @MainActor [weak self] in
            do { try await body() }
            catch {
                let e = UserFacingError.from(error)
                self?.phase = e.isCancellation ? .cancelled : .failed(e)
            }
        }
    }

    /// Resolves when the current step has ended (used by tests and by callers that chain steps).
    public func waitUntilIdle() async { await task?.value }
}
