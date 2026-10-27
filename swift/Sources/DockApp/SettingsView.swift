#if canImport(SwiftUI)
import SwiftUI
import DockCore

@MainActor struct SettingsView: View {
    @Environment(AppModel.self) private var model

    var body: some View {
        @Bindable var model = model
        Form {
            Section("Performance") {
                Stepper(value: $model.settings.cores, in: 1...ProcessInfo.processInfo.activeProcessorCount) {
                    Text("CPU cores to use: \(model.settings.effectiveCores())")
                }
                Text("This Mac has \(ProcessInfo.processInfo.activeProcessorCount) cores. Several ligands run side by side, one core each; a single ligand uses all chosen cores.")
                    .font(.callout).foregroundStyle(.secondary)
            }
            Section("Docking engine") {
                TextField("Path to the Vina binary (leave empty for the bundled engine)",
                          text: Binding(get: { model.settings.enginePath ?? "" },
                                        set: { model.settings.enginePath = $0.isEmpty ? nil : $0 }))
                Label(model.engineFound ? "Engine found" : "No engine found: docking is not possible until one is chosen",
                      systemImage: model.engineFound ? "checkmark.circle" : "exclamationmark.triangle")
                    .foregroundStyle(model.engineFound ? Color.green : Color.orange)
            }
        }
        .formStyle(.grouped)
        .frame(width: 520)
        .padding()
    }
}
#endif
