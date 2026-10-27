import Foundation

/// Finds the bundled chemistry environment and docking engine (plan §15). In a shipped app they live inside the
/// bundle; during development, environment variables or the system Python are used.
public struct HelperLocator: Sendable {
    public var resourceDirectory: URL?
    public var environment: [String: String]
    /// Writable folder for caches (matplotlib fonts etc.); the bundle itself must never be written to.
    public var cacheDirectory: URL?
    public var fileExists: @Sendable (String) -> Bool

    public init(resourceDirectory: URL? = nil, environment: [String: String] = ProcessInfo.processInfo.environment,
                cacheDirectory: URL? = nil,
                fileExists: @escaping @Sendable (String) -> Bool = { FileManager.default.isExecutableFile(atPath: $0) }) {
        self.resourceDirectory = resourceDirectory; self.environment = environment
        self.cacheDirectory = cacheDirectory; self.fileExists = fileExists
    }

    /// Bundled Python first, then $DOCKCHEM_PYTHON, then `python3` from the path (development).
    public func launch() -> HelperClient.Launch {
        var env: [String: String] = ["PYTHONDONTWRITEBYTECODE": "1"]     // never write .pyc into a signed bundle
        if let cache = cacheDirectory { env["MPLCONFIGDIR"] = cache.appendingPathComponent("matplotlib").path }
        if let engine = enginePath() { env["DOCKCHEM_VINA"] = engine }
        if let r = resourceDirectory {
            let py = r.appendingPathComponent("python/bin/python3").path
            if fileExists(py) {
                env["PYTHONHOME"] = r.appendingPathComponent("python").path
                env["PYTHONNOUSERSITE"] = "1"
                return HelperClient.Launch(executable: URL(fileURLWithPath: py), arguments: ["-m", "dockchem.server"], environment: env)
            }
        }
        if let py = environment["DOCKCHEM_PYTHON"], fileExists(py) {
            return HelperClient.Launch(executable: URL(fileURLWithPath: py), arguments: ["-m", "dockchem.server"], environment: env)
        }
        return .development(environment: env)
    }

    /// Bundled engine, then $DOCKCHEM_VINA; nil lets the helper search (its own order) and report if it finds none.
    public func enginePath() -> String? {
        if let r = resourceDirectory {
            let p = r.appendingPathComponent("engine/vina").path
            if fileExists(p) { return p }
        }
        if let p = environment["DOCKCHEM_VINA"], fileExists(p) { return p }
        return nil
    }
}

/// Finds the 3D viewer page: inside the app bundle, or (development) next to the sources.
public struct ViewerLocator: Sendable {
    public var resourceDirectory: URL?
    public var environment: [String: String]
    public var workingDirectory: URL
    public var fileExists: @Sendable (String) -> Bool

    public init(resourceDirectory: URL? = nil, environment: [String: String] = ProcessInfo.processInfo.environment,
                workingDirectory: URL = URL(fileURLWithPath: FileManager.default.currentDirectoryPath),
                fileExists: @escaping @Sendable (String) -> Bool = { FileManager.default.fileExists(atPath: $0) }) {
        self.resourceDirectory = resourceDirectory; self.environment = environment
        self.workingDirectory = workingDirectory; self.fileExists = fileExists
    }

    /// The directory holding index.html, or nil if the viewer is not installed.
    public func directory() -> URL? {
        var candidates: [URL] = []
        if let r = resourceDirectory { candidates.append(r.appendingPathComponent("viewer")) }
        if let e = environment["DOCK_VIEWER_DIR"] { candidates.append(URL(fileURLWithPath: e)) }
        candidates.append(workingDirectory.appendingPathComponent("viewer"))
        return candidates.first { fileExists($0.appendingPathComponent("index.html").path) }
    }
}
