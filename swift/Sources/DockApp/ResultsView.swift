#if canImport(SwiftUI)
import SwiftUI
import AppKit
import DockCore

struct ResultsView: View {
    let session: QuickDockSession
    let result: DockRunResult
    @State private var selection: DockedLigand.ID?
    @State private var pose = 1

    private var rows: [DockedLigand] { result.ranked }
    private var selected: DockedLigand? { rows.first { $0.id == selection } ?? rows.first }

    var body: some View {
        GroupBox("Results") {
            VStack(alignment: .leading, spacing: 10) {
                HStack {
                    Text("\(result.docked) of \(result.ligands) docked" + (result.failed > 0 ? " · \(result.failed) failed" : ""))
                    Spacer()
                    Button("Show files") {
                        NSWorkspace.shared.activateFileViewerSelecting([URL(fileURLWithPath: result.csvPath)])
                    }
                }
                Table(rows, selection: $selection) {
                    TableColumn("Rank") { Text($0.rank.map(String.init) ?? "–") }.width(50)
                    TableColumn("Ligand", value: \.name)
                    TableColumn("Score (kcal/mol)") { Text($0.affinity.map { String(format: "%.2f", $0) } ?? "–") }.width(110)
                    TableColumn("Poses") { Text(String($0.scores.count)) }.width(50)
                    TableColumn("Status") { Text($0.status == .docked ? "docked" : ($0.error.isEmpty ? $0.status.rawValue : $0.error)) }
                }
                .frame(minHeight: min(CGFloat(rows.count) * 28 + 40, 260))

                if let sel = selected, sel.status == .docked {
                    PoseView(session: session, ligand: sel, pose: $pose)
                        .id(sel.id)
                }
                Text(result.scoresAre).font(.footnote).foregroundStyle(.secondary)
            }
        }
        .onChange(of: selection) { _, _ in pose = 1 }
    }
}

struct PoseView: View {
    let session: QuickDockSession
    let ligand: DockedLigand
    @Binding var pose: Int
    @State private var image: NSImage?
    @State private var failure: String?

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack {
                Stepper("Pose \(pose) of \(max(ligand.scores.count, 1))", value: $pose, in: 1...max(ligand.scores.count, 1))
                if pose <= ligand.scores.count {
                    Text(String(format: "%.2f kcal/mol", ligand.scores[pose - 1])).foregroundStyle(.secondary)
                }
            }
            ZStack {
                if let image {
                    Image(nsImage: image).resizable().scaledToFit()
                } else if let failure {
                    Text(failure).foregroundStyle(.secondary).textSelection(.enabled)
                } else {
                    ProgressView()
                }
            }
            .frame(maxWidth: .infinity, minHeight: 280)
        }
        .task(id: pose) {
            image = nil; failure = nil
            do {
                let url = try await session.interactionImage(for: ligand, pose: pose)
                image = NSImage(contentsOf: url)
                if image == nil { failure = "The interaction image could not be read." }
            } catch {
                failure = UserFacingError.from(error).detail
            }
        }
    }
}
#endif
