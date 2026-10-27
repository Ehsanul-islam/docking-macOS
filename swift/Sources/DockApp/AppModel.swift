#if canImport(SwiftUI)
import SwiftUI
import Observation
import DockCore

/// App-wide state: the chemistry helper, the current run and the user's settings.
@MainActor @Observable
final class AppModel {
    var settings: AppSettings {
        didSet { Self.save(settings) }
    }
    private(set) var session: QuickDockSession
    private let client: HelperClient
    private let locator: HelperLocator

    init() {
        let locator = HelperLocator(resourceDirectory: Bundle.main.resourceURL)
        let client = HelperClient(launch: locator.launch())
        let settings = Self.load()
        self.locator = locator
        self.client = client
        self.settings = settings
        self.session = QuickDockSession(helper: client, workDirectory: Self.newRunDirectory(), settings: settings)
    }

    /// Every run gets its own folder under Application Support, so nothing is overwritten and a run can be reopened.
    private static func newRunDirectory() -> URL {
        let base = FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
            .appendingPathComponent("Dock/Runs", isDirectory: true)
        let stamp = ISO8601DateFormatter().string(from: Date()).replacingOccurrences(of: ":", with: "-")
        return base.appendingPathComponent(stamp, isDirectory: true)
    }

    func newRun() {
        session.reset()
        session = QuickDockSession(helper: client, workDirectory: Self.newRunDirectory(), settings: settings)
    }

    func run() {
        session.settings = settings
        session.runAll()
    }

    var engineFound: Bool { settings.enginePath != nil || locator.enginePath() != nil }

    private static let key = "dock.settings.v1"
    private static func load() -> AppSettings {
        guard let data = UserDefaults.standard.data(forKey: key),
              let s = try? JSONDecoder().decode(AppSettings.self, from: data) else { return AppSettings() }
        return s
    }
    private static func save(_ s: AppSettings) {
        if let data = try? JSONEncoder().encode(s) { UserDefaults.standard.set(data, forKey: key) }
    }
}
#endif
