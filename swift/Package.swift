// swift-tools-version:5.9
import PackageDescription

let package = Package(
    name: "DockCore",
    platforms: [.macOS(.v14)],
    products: [.library(name: "DockCore", targets: ["DockCore"])],
    targets: [
        .target(name: "DockCore"),
        .testTarget(name: "DockCoreTests", dependencies: ["DockCore"]),
    ]
)
