import Foundation

// Typed views of the helper's results. They read fields tolerantly: a field a newer or older helper does not send
// becomes a default, never a crash, and unknown fields are ignored.

public enum Verdict: String, Sendable, Equatable {
    case ready, fixable, unusable, unknown
}

public struct Finding: Sendable, Equatable {
    public enum Severity: String, Sendable { case note, warning, problem, unknown }
    public let severity: Severity
    public let code: String
    public let message: String
}

/// `structure.assess`: is this structure fit to dock into, and what is wrong with it.
public struct Assessment: Sendable, Equatable {
    public let verdict: Verdict
    public let summary: String
    public let findings: [Finding]
    public let proteinResidues: Int
    public let atoms: Int
    public let resolution: Double?

    public init?(json: JSONValue) {
        guard let v = json["verdict"]?.string else { return nil }
        verdict = Verdict(rawValue: v) ?? .unknown
        summary = json["summary"]?.string ?? ""
        findings = (json["findings"]?.array ?? []).map {
            Finding(severity: Finding.Severity(rawValue: $0["severity"]?.string ?? "") ?? .unknown,
                    code: $0["code"]?.string ?? "", message: $0["message"]?.string ?? "")
        }
        proteinResidues = json["stats"]?["protein_residues"]?.int ?? 0
        atoms = json["stats"]?["atoms"]?.int ?? 0
        resolution = json["stats"]?["resolution"]?.double
    }

    /// Findings the user should read first.
    public var attention: [Finding] { findings.filter { $0.severity == .problem || $0.severity == .warning } }
}

public struct SearchBox: Sendable, Equatable {
    public var center: [Double]
    public var size: [Double]
    public init(center: [Double], size: [Double]) { self.center = center; self.size = size }
    public init?(json: JSONValue) {
        guard let c = json["center"]?.array?.compactMap(\.double), c.count == 3,
              let s = json["size"]?.array?.compactMap(\.double), s.count == 3 else { return nil }
        self.init(center: c, size: s)
    }
    public var volume: Double { size.reduce(1, *) }
    public var longestSide: Double { size.max() ?? 0 }
}

public struct EngineSettings: Sendable, Equatable {
    public let scoring: String
    public let cores: Int
    public let exhaustiveness: Int
    public let numModes: Int
    public let seed: Int
    public let presetLabel: String
}

/// `quickprep.run`: everything the quick-prep panel decided, so the UI can show every decision.
public struct QuickPrepResult: Sendable, Equatable {
    public let preset: String
    public let warnings: [String]
    public let assessment: Assessment?
    public let receptorStatus: String
    public let receptorPDBQT: String
    public let receptorPDB: String
    public let leftOutResidues: [String]
    public let siteDescription: String
    public let box: SearchBox
    public let ligandsPrepared: Int
    public let ligandsFailed: [String]
    public let ligandFiles: [String]
    public let engine: EngineSettings
    public let readyToDock: Bool
    public let seconds: Double

    /// True when no binding site was given or found.
    public var isBlind: Bool { siteDescription.lowercased().hasPrefix("no binding site") }

    public init?(json: JSONValue) {
        guard let box = json["site"].flatMap(SearchBox.init(json:)), let rec = json["receptor"],
              let eng = json["engine"] else { return nil }
        self.box = box
        preset = json["preset"]?.string ?? "standard"
        warnings = (json["warnings"]?.array ?? []).compactMap(\.string)
        assessment = json["assessment"].flatMap(Assessment.init(json:))
        receptorStatus = rec["status"]?.string ?? "unknown"
        receptorPDBQT = rec["pdbqt"]?.string ?? ""
        receptorPDB = rec["pdb"]?.string ?? ""
        leftOutResidues = (rec["left_out"]?.array ?? []).map { $0["residue"]?.string ?? $0.string ?? "" }.filter { !$0.isEmpty }
        siteDescription = json["site"]?["how"]?.string ?? ""
        ligandsPrepared = json["ligands"]?["prepared"]?.int ?? 0
        ligandsFailed = (json["ligands"]?["failed"]?.array ?? []).map { $0["name"]?.string ?? $0.string ?? "?" }
        ligandFiles = (json["ligands"]?["files"]?.array ?? []).compactMap(\.string)
        engine = EngineSettings(scoring: eng["scoring"]?.string ?? "vina", cores: eng["cores"]?.int ?? 1,
                                exhaustiveness: eng["exhaustiveness"]?.int ?? 32, numModes: eng["num_modes"]?.int ?? 9,
                                seed: eng["seed"]?.int ?? 42, presetLabel: eng["preset_label"]?.string ?? "")
        readyToDock = json["ready_to_dock"]?.bool ?? false
        seconds = json["seconds"]?.double ?? 0
    }
}

/// One ligand's docking outcome (a row of results.json).
public struct DockedLigand: Sendable, Equatable, Identifiable {
    public enum Status: String, Sendable { case docked, failed, notRun = "not_run", unknown }
    public var id: String { name }
    public let name: String
    public let status: Status
    public let error: String
    public let seconds: Double
    public let affinity: Double?
    public let scores: [Double]
    public let posesPDBQT: String
    public let posesSDF: String?
    public let rank: Int?

    init(json: JSONValue) {
        name = json["name"]?.string ?? "?"
        status = Status(rawValue: json["status"]?.string ?? "") ?? .unknown
        error = json["error"]?.string ?? ""
        seconds = json["seconds"]?.double ?? 0
        affinity = json["affinity"]?.double
        scores = (json["scores"]?.array ?? []).compactMap(\.double)
        posesPDBQT = json["poses"]?.string ?? ""
        posesSDF = json["sdf"]?.string
        rank = json["rank"]?.int
    }
}

/// `dock.run`: the ranked outcome of a batch.
public struct DockRunResult: Sendable, Equatable {
    public let ligands: Int
    public let docked: Int
    public let failed: Int
    public let notRun: Int
    public let cores: Int
    public let results: [DockedLigand]
    public let csvPath: String
    public let outdir: String
    public let scoresAre: String

    public init?(json: JSONValue) {
        guard let rows = json["results"]?.array else { return nil }
        results = rows.map(DockedLigand.init(json:))
        ligands = json["ligands"]?.int ?? rows.count
        docked = json["docked"]?.int ?? results.filter { $0.status == .docked }.count
        failed = json["failed"]?.int ?? 0
        notRun = json["not_run"]?.int ?? 0
        cores = json["cores"]?.int ?? 1
        csvPath = json["csv"]?.string ?? ""
        outdir = json["outdir"]?.string ?? ""
        scoresAre = json["scores_are"]?.string ?? ""
    }

    /// Ranked, best (most negative) first; ligands that did not dock come last.
    public var ranked: [DockedLigand] {
        results.sorted {
            switch ($0.affinity, $1.affinity) {
            case let (a?, b?): return a < b
            case (_?, nil): return true
            case (nil, _?): return false
            default: return $0.name < $1.name
            }
        }
    }
}

extension HelperResult {
    public func decode<T>(_ make: (JSONValue) -> T?, as name: String) throws -> T {
        guard let v = make(result) else {
            throw HelperClientError.malformedResponse("the \(name) result is missing expected fields")
        }
        return v
    }
}
