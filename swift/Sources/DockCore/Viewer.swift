import Foundation
import Observation

/// A command for the 3D viewer page (`window.dock.command`, see viewer/README.md). Built here as JSON so the
/// transport (a WKWebView in the app, a recording mock in tests) only has to deliver a string.
public struct ViewerCommand: Sendable, Equatable {
    public let name: String
    public let json: String

    init(_ name: String, _ fields: [String: JSONValue] = [:]) {
        var f = fields
        f["cmd"] = .string(name)
        self.name = name
        let data = (try? JSONEncoder().encode(JSONValue.object(f))) ?? Data("{}".utf8)
        json = String(decoding: data, as: UTF8.self)
    }

    public enum Theme: String, Sendable { case dark, light }

    public static func setTheme(_ t: Theme) -> ViewerCommand { .init("setTheme", ["theme": .string(t.rawValue)]) }
    public static func loadReceptor(pdb: Data, name: String) -> ViewerCommand {
        .init("loadReceptor", ["data": .string(pdb.base64EncodedString()), "format": "pdb", "name": .string(name)])
    }
    public static func loadPose(sdf: Data, index: Int, score: Double?, poses: [(index: Int, score: Double?)]) -> ViewerCommand {
        .init("loadPose", ["data": .string(sdf.base64EncodedString()), "format": "sdf", "index": .int(index),
                           "score": score.map(JSONValue.double) ?? .null,
                           "poses": .array(poses.map { ["index": .int($0.index), "score": $0.score.map(JSONValue.double) ?? .null] })])
    }
    public static func setInteractions(_ records: [JSONValue]) -> ViewerCommand { .init("setInteractions", ["interactions": .array(records)]) }
    public static func setBox(_ box: SearchBox?) -> ViewerCommand {
        guard let box else { return .init("setBox", ["center": .null, "size": .null]) }
        return .init("setBox", ["center": .array(box.center.map(JSONValue.double)), "size": .array(box.size.map(JSONValue.double))])
    }
    public static func setOptions(_ options: [String: Bool]) -> ViewerCommand { .init("setOptions", options.mapValues(JSONValue.bool)) }
    public static func focusLigand() -> ViewerCommand { .init("focus", ["target": "ligand"]) }
    public static let screenshot = ViewerCommand("screenshot")
    public static let clear = ViewerCommand("clear")
}

/// What the viewer page reports back.
public enum ViewerEvent: Sendable, Equatable {
    case ready
    case rendered
    case poseSelected(Int)
    case screenshot(Data)
    case warning(String)
    case error(String)

    /// From the message object posted by the page (already converted to JSON by the web view).
    public init?(_ json: JSONValue) {
        switch json["type"]?.string {
        case "ready": self = .ready
        case "rendered": self = .rendered
        case "poseSelected": guard let i = json["index"]?.int else { return nil }; self = .poseSelected(i)
        case "screenshot":
            guard let uri = json["dataUri"]?.string, let comma = uri.firstIndex(of: ","),
                  let data = Data(base64Encoded: String(uri[uri.index(after: comma)...])) else { return nil }
            self = .screenshot(data)
        case "warning": self = .warning(json["message"]?.string ?? "")
        case "error": self = .error(json["message"]?.string ?? "The 3D view reported an error.")
        default: return nil
        }
    }
}

/// Delivers commands to the page. Implemented by the WKWebView wrapper in the app.
@MainActor public protocol ViewerTransport: AnyObject {
    func send(_ command: ViewerCommand)
}

/// Drives the viewer from docking results: loads the prepared receptor once, shows one pose at a time with its
/// interactions, and queues everything until the page has reported `ready`.
@MainActor @Observable
public final class ViewerController {
    public private(set) var isReady = false
    public private(set) var lastError: String?
    public private(set) var poseIndex: Int?
    public private(set) var isBusy = false
    public var theme: ViewerCommand.Theme = .dark { didSet { send(.setTheme(theme)) } }
    public var onScreenshot: (@MainActor (Data) -> Void)?

    private weak var transport: ViewerTransport?
    private var queue: [ViewerCommand] = []
    private var loadedReceptor: String?
    private var shown: (prepared: QuickPrepResult, ligand: DockedLigand)?
    private let helper: HelperClient?

    public init(helper: HelperClient?) { self.helper = helper }

    public func attach(_ transport: ViewerTransport) { self.transport = transport }

    // MARK: sending

    public func send(_ command: ViewerCommand) {
        guard isReady, let transport else { queue.append(command); return }
        transport.send(command)
    }

    public func handle(_ event: ViewerEvent) {
        switch event {
        case .ready:
            isReady = true
            if let t = transport {
                t.send(.setTheme(theme))
                queue.forEach { t.send($0) }
            }
            queue.removeAll()
        case .rendered: isBusy = false
        case .poseSelected(let i):
            if let shown { Task { await showPose(i, of: shown.ligand, prepared: shown.prepared) } }
        case .screenshot(let data): onScreenshot?(data)
        case .warning: break
        case .error(let m): lastError = m; isBusy = false
        }
    }

    /// The page crashed or was reloaded: everything must be sent again.
    public func pageReloaded() {
        isReady = false; loadedReceptor = nil
    }

    // MARK: content

    /// Just a structure (before docking).
    public func showStructure(at url: URL) {
        guard let data = try? Data(contentsOf: url) else { lastError = "The structure file could not be read."; return }
        loadedReceptor = nil
        send(.loadReceptor(pdb: data, name: url.deletingPathExtension().lastPathComponent))
    }

    /// The prepared receptor, the search box, and the chosen pose of a docked ligand with its interactions.
    public func showResult(prepared: QuickPrepResult, ligand: DockedLigand, pose: Int = 1, structureName: String = "") async {
        shown = (prepared, ligand)
        if loadedReceptor != prepared.receptorPDB {
            guard let data = try? Data(contentsOf: URL(fileURLWithPath: prepared.receptorPDB)) else {
                lastError = "The prepared receptor could not be read."; return
            }
            loadedReceptor = prepared.receptorPDB
            send(.loadReceptor(pdb: data, name: structureName))
            send(.setBox(prepared.box))
        }
        await showPose(pose, of: ligand, prepared: prepared)
    }

    public func showPose(_ index: Int, of ligand: DockedLigand, prepared: QuickPrepResult) async {
        guard let sdfPath = ligand.posesSDF, let sdf = try? Data(contentsOf: URL(fileURLWithPath: sdfPath)) else {
            lastError = "The docked poses of \(ligand.name) could not be read."; return
        }
        isBusy = true
        poseIndex = index
        let poses = ligand.scores.enumerated().map { (index: $0.offset + 1, score: Optional($0.element)) }
        let score = index >= 1 && index <= ligand.scores.count ? ligand.scores[index - 1] : nil
        send(.loadPose(sdf: sdf, index: index, score: score, poses: poses))
        send(.setInteractions([]))
        guard let helper else { return }
        do {
            let r = try await helper.call("complex.interactions", [
                "receptor": .string(prepared.receptorPDB), "pose": .string(sdfPath), "pose_index": .int(index),
                "affinity": score.map(JSONValue.double) ?? .null, "compound": .string(ligand.name)])
            // a newer selection may have arrived while the helper worked
            guard poseIndex == index, shown?.ligand.name == ligand.name else { return }
            send(.setInteractions(r.result["interactions"]?.array ?? []))
        } catch {
            lastError = UserFacingError.from(error).detail
        }
    }

    public func clear() {
        shown = nil; poseIndex = nil; loadedReceptor = nil
        send(.clear)
    }
}

extension ViewerController {
    /// After preparation, before docking: the prepared receptor (with its added hydrogens) and the search box.
    public func showPrepared(_ prepared: QuickPrepResult, structureName: String) {
        guard let data = try? Data(contentsOf: URL(fileURLWithPath: prepared.receptorPDB)) else { return }
        loadedReceptor = prepared.receptorPDB
        send(.loadReceptor(pdb: data, name: structureName))
        send(.setBox(prepared.box))
    }
}
