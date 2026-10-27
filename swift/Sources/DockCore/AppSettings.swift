import Foundation

/// User settings that affect a run. Stored by the app; validated here so a bad value can never reach the helper.
public struct AppSettings: Codable, Equatable, Sendable {
    public static let defaultCores = 2

    public var cores: Int
    public var enginePath: String?

    public init(cores: Int = AppSettings.defaultCores, enginePath: String? = nil) {
        self.cores = cores; self.enginePath = enginePath
    }

    /// Cores to ask for: at least 1 and at most what the machine has.
    public func effectiveCores(available: Int = ProcessInfo.processInfo.activeProcessorCount) -> Int {
        min(max(cores, 1), max(available, 1))
    }
}
