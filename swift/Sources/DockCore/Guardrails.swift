import Foundation

/// Warnings shown before a run (plan §9.3). They never block: the user can always go ahead.
public struct Guardrail: Sendable, Equatable, Identifiable {
    public enum Severity: Int, Sendable, Comparable {
        case info, caution
        public static func < (a: Severity, b: Severity) -> Bool { a.rawValue < b.rawValue }
    }
    public let id: String
    public let severity: Severity
    public let message: String
}

public enum Guardrails {
    public static let largeBoxVolume = 27_000.0     // A^3, Vina prints the same warning
    public static let blindMinExhaustiveness = 32
    public static let manyTorsions = 20

    public static func check(box: SearchBox, exhaustiveness: Int, blind: Bool, maxTorsions: Int? = nil,
                             ligandsFailed: Int = 0, receptorLeftOut: Int = 0) -> [Guardrail] {
        var out: [Guardrail] = []
        if box.volume > largeBoxVolume {
            out.append(Guardrail(id: "large-box", severity: .caution, message:
                "Large search box (\(Int(box.volume.rounded())) Å³): increase exhaustiveness or narrow the site."))
        }
        if blind && exhaustiveness < blindMinExhaustiveness {
            out.append(Guardrail(id: "blind-undersampled", severity: .caution, message:
                "Blind docking is under-sampled at exhaustiveness \(exhaustiveness); 32 or more is advised."))
        }
        if let t = maxTorsions, t > manyTorsions {
            out.append(Guardrail(id: "many-torsions", severity: .caution, message:
                "A ligand has \(t) rotatable bonds: the search may be unreliable; consider more exhaustiveness."))
        }
        if ligandsFailed > 0 {
            out.append(Guardrail(id: "ligands-failed", severity: .caution, message:
                "\(ligandsFailed) ligand(s) could not be prepared and will not be docked."))
        }
        if receptorLeftOut > 0 {
            out.append(Guardrail(id: "receptor-left-out", severity: .info, message:
                "\(receptorLeftOut) receptor residue(s) were left out of the docking model (see the prep report)."))
        }
        return out.sorted { $0.severity > $1.severity }
    }
}
