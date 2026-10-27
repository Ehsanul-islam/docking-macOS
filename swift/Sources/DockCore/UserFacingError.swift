import Foundation

/// An error as the user should see it: what happened, why, and what they can do. Built from the helper's typed
/// errors so every failure path ends in a sentence a non-expert can act on.
public struct UserFacingError: Error, Sendable, Equatable, LocalizedError {
    public enum Action: Sendable, Equatable {
        case none
        case prepareAnyway               // UnusableStructure: re-run with allow_unusable after the user confirms
        case provideBondOrders(String)   // NeedsTemplate: ask for a SMILES for this residue
        case chooseEngine                // engine not found: open Settings
        case copyDetails                 // a bug: let the user copy the log
    }
    public let title: String
    public let detail: String
    public let suggestion: String
    public let action: Action
    public let log: String
    public var errorDescription: String? { "\(title). \(detail)" }

    public init(title: String, detail: String, suggestion: String = "", action: Action = .none, log: String = "") {
        self.title = title; self.detail = detail; self.suggestion = suggestion; self.action = action; self.log = log
    }

    public static func from(_ error: Error) -> UserFacingError {
        if let e = error as? UserFacingError { return e }
        if let e = error as? RPCError { return from(rpc: e) }
        if let e = error as? HelperClientError {
            switch e {
            case .helperExited(_, let tail):
                return UserFacingError(title: "The chemistry helper stopped", detail: "It ended unexpectedly while working.",
                                       suggestion: "Try again. If it keeps happening, copy the details and report it.",
                                       action: .copyDetails, log: tail)
            case .failedToStart(let why):
                return UserFacingError(title: "The chemistry helper could not start", detail: why,
                                       suggestion: "Reinstall the app, or check Settings for the helper location.", action: .copyDetails)
            case .notRunning, .malformedResponse:
                return UserFacingError(title: "The chemistry helper did not answer properly", detail: e.localizedDescription,
                                       suggestion: "Try again.", action: .copyDetails)
            }
        }
        if let e = error as? EngineRunnerError, e == .cancelled { return cancelled }
        if error is CancellationError { return cancelled }
        return UserFacingError(title: "Something went wrong", detail: error.localizedDescription, action: .copyDetails)
    }

    public static let cancelled = UserFacingError(title: "Cancelled", detail: "The run was stopped.")

    public var isCancellation: Bool { self == Self.cancelled }

    private static func from(rpc e: RPCError) -> UserFacingError {
        let log = e.data?.log ?? ""
        switch e.kind {
        case HelperErrorKind.cancelled:
            return .cancelled
        case HelperErrorKind.unusableStructure:
            return UserFacingError(title: "This structure can't be docked into", detail: e.message,
                                   suggestion: "Choose another structure, or prepare it anyway if you understand the risk.",
                                   action: .prepareAnyway, log: log)
        case HelperErrorKind.needsTemplate:
            let res = e.data?.extra["resname"]?.string ?? "a molecule"
            return UserFacingError(title: "Bond orders are needed for \(res)", detail: e.message,
                                   suggestion: "Provide its SMILES so the chemistry can be set up correctly.",
                                   action: .provideBondOrders(res), log: log)
        case HelperErrorKind.formatLimit:
            return UserFacingError(title: "This structure is too large for the PDB format", detail: e.message,
                                   suggestion: "Dock a single chain or a region instead.", log: log)
        case HelperErrorKind.invalidInput, HelperErrorKind.selectionError:
            return UserFacingError(title: "Check the input", detail: e.message, log: log)
        case HelperErrorKind.unsupportedFeature:
            let engine = e.message.lowercased().contains("engine")
            return UserFacingError(title: engine ? "The docking engine was not found" : "Not supported yet", detail: e.message,
                                   suggestion: engine ? "Choose the engine in Settings." : "", action: engine ? .chooseEngine : .none, log: log)
        case HelperErrorKind.toolFailure:
            return UserFacingError(title: "A preparation tool failed", detail: e.message,
                                   suggestion: "The details can be copied for a bug report.", action: .copyDetails, log: log)
        default:
            return UserFacingError(title: "Internal error", detail: e.message,
                                   suggestion: "This is a bug in the app. Copy the details and report it.",
                                   action: .copyDetails, log: log)
        }
    }
}
