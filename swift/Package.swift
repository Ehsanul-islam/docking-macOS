// swift-tools-version:5.9
import PackageDescription

let package = Package(
    name: "DockCore",
    platforms: [.macOS(.v14)],
    products: [.library(name: "DockCore", targets: ["DockCore"]), .executable(name: "DockApp", targets: ["DockApp"])],
    targets: [
        .target(name: "DockCore"),
        // The SwiftUI app. Its sources are guarded with canImport(SwiftUI), so the package still builds on Linux CI.
        .executableTarget(name: "DockApp", dependencies: ["DockCore"]),
        .testTarget(name: "DockCoreTests", dependencies: ["DockCore"], resources: [.copy("Fixtures")]),
    ]
)
