#if canImport(SwiftUI)
import SwiftUI
import DockCore

@main @MainActor
struct DockApp: App {
    @State private var model = AppModel()

    var body: some Scene {
        WindowGroup("Dock") {
            ContentView()
                .environment(model)
                .frame(minWidth: 760, minHeight: 640)
        }
        Settings {
            SettingsView()
                .environment(model)
        }
    }
}
#else
@main
struct DockApp {
    static func main() { print("DockApp needs macOS (SwiftUI).") }
}
#endif
