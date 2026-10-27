#if canImport(SwiftUI)
import SwiftUI
import AppKit
import UniformTypeIdentifiers
import DockCore

struct ContentView: View {
    @Environment(AppModel.self) private var model

    var body: some View {
        let session = model.session
        HSplitView {
            ScrollView {
                VStack(alignment: .leading, spacing: 16) {
                    ReceptorCard(session: session)
                    LigandCard(session: session)
                    OptionsCard(session: session)
                    RunBar(session: session)
                    MessagesView(session: session)
                    if let result = session.result {
                        ResultsView(session: session, result: result)
                    }
                }
                .padding(20)
            }
            .frame(minWidth: 380, idealWidth: 440, maxWidth: 560)

            ViewerPanel()
                .frame(minWidth: 380)
        }
        .toolbar {
            ToolbarItem {
                Button("New run") { model.newRun() }
                    .disabled(session.isBusy)
            }
        }
        .onChange(of: session.receptor) { _, url in
            if let url { model.viewer.showStructure(at: url) }
        }
        .onChange(of: session.prepared) { _, prepared in
            guard let prepared else { return }
            model.viewer.showPrepared(prepared, structureName: session.receptor?.deletingPathExtension().lastPathComponent ?? "")
        }
    }
}

// MARK: structure

struct ReceptorCard: View {
    @Bindable var session: QuickDockSession
    @State private var choosing = false
    @State private var targeted = false

    var body: some View {
        GroupBox("1. Structure") {
            VStack(alignment: .leading, spacing: 10) {
                HStack {
                    Image(systemName: "doc.badge.plus")
                    Text(session.receptor?.lastPathComponent ?? "Drop a PDB or mmCIF file here")
                        .foregroundStyle(session.receptor == nil ? .secondary : .primary)
                    Spacer()
                    Button("Choose…") { choosing = true }
                        .disabled(session.isBusy)
                }
                .padding(10)
                .frame(maxWidth: .infinity)
                .background(RoundedRectangle(cornerRadius: 8).strokeBorder(targeted ? Color.accentColor : Color.secondary.opacity(0.4),
                                                                             style: StrokeStyle(lineWidth: 1, dash: [5])))
                .dropDestination(for: URL.self) { urls, _ in
                    guard let url = urls.first, !session.isBusy else { return false }
                    set(url)
                    return true
                } isTargeted: { targeted = $0 }

                if let a = session.assessment { AssessmentView(assessment: a) }
            }
        }
        .fileImporter(isPresented: $choosing, allowedContentTypes: [.item]) { result in
            if case .success(let url) = result { set(url) }
        }
    }

    private func set(_ url: URL) {
        session.receptor = url
        session.assess()
    }
}

struct AssessmentView: View {
    let assessment: Assessment

    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack(spacing: 8) {
                Label(title, systemImage: symbol)
                    .foregroundStyle(color)
                    .font(.headline)
                Text("\(assessment.proteinResidues) residues" + (assessment.resolution.map { String(format: " · %.1f Å", $0) } ?? ""))
                    .foregroundStyle(.secondary)
            }
            if !assessment.summary.isEmpty { Text(assessment.summary).foregroundStyle(.secondary) }
            ForEach(Array(assessment.attention.enumerated()), id: \.offset) { _, f in
                Label(f.message, systemImage: f.severity == .problem ? "exclamationmark.octagon" : "exclamationmark.triangle")
                    .font(.callout)
            }
        }
    }

    private var title: String {
        switch assessment.verdict {
        case .ready: return "Ready to dock"
        case .fixable: return "Needs repairs (done automatically)"
        case .unusable: return "Can't be docked into"
        case .unknown: return "Checked"
        }
    }
    private var symbol: String {
        switch assessment.verdict {
        case .ready: return "checkmark.seal"
        case .fixable: return "wrench.and.screwdriver"
        case .unusable: return "xmark.octagon"
        case .unknown: return "questionmark.circle"
        }
    }
    private var color: Color {
        switch assessment.verdict {
        case .ready: return .green
        case .fixable: return .orange
        case .unusable: return .red
        case .unknown: return .secondary
        }
    }
}

// MARK: ligand

struct LigandCard: View {
    @Bindable var session: QuickDockSession
    @State private var choosing = false

    var body: some View {
        GroupBox("2. Ligand(s)") {
            HStack {
                TextField("SMILES, or choose a file (SDF, MOL, SMI, PDB; several molecules = a library)", text: $session.ligand)
                    .textFieldStyle(.roundedBorder)
                    .disabled(session.isBusy)
                Button("Choose…") { choosing = true }
                    .disabled(session.isBusy)
            }
            .padding(.vertical, 4)
            .dropDestination(for: URL.self) { urls, _ in
                guard let url = urls.first, !session.isBusy else { return false }
                session.ligand = url.path
                return true
            }
        }
        .fileImporter(isPresented: $choosing, allowedContentTypes: [.item]) { result in
            if case .success(let url) = result { session.ligand = url.path }
        }
    }
}

// MARK: options

struct OptionsCard: View {
    @Bindable var session: QuickDockSession
    @State private var siteKind = SiteKind.auto
    @State private var ligandName = ""
    @State private var residues = ""
    @State private var center = ["", "", ""]
    @State private var size = ["20", "20", "20"]

    enum SiteKind: String, CaseIterable, Identifiable {
        case auto = "Automatic", ligand = "Around a ligand", residues = "Around residues", box = "Box"
        var id: String { rawValue }
    }

    var body: some View {
        GroupBox("3. Settings") {
            VStack(alignment: .leading, spacing: 10) {
                Picker("Search effort", selection: $session.preset) {
                    Text("Fast").tag(DockPreset.fast)
                    Text("Standard (recommended)").tag(DockPreset.standard)
                    Text("Thorough").tag(DockPreset.thorough)
                }
                .pickerStyle(.segmented)
                .disabled(session.isBusy)

                DisclosureGroup("Binding site") {
                    VStack(alignment: .leading, spacing: 8) {
                        Picker("Where to dock", selection: $siteKind) {
                            ForEach(SiteKind.allCases) { Text($0.rawValue).tag($0) }
                        }
                        .onChange(of: siteKind) { _, _ in apply() }
                        switch siteKind {
                        case .auto:
                            Text("Uses the largest drug-like ligand in the structure; blind docking if there is none.")
                                .font(.callout).foregroundStyle(.secondary)
                        case .ligand:
                            TextField("Ligand, e.g. OHT or A:OHT:600", text: $ligandName).onSubmit(apply)
                        case .residues:
                            TextField("Residues, e.g. A:394 A:353", text: $residues).onSubmit(apply)
                        case .box:
                            HStack {
                                ForEach(0..<3, id: \.self) { i in TextField("c\(["x", "y", "z"][i])", text: $center[i]) }
                                ForEach(0..<3, id: \.self) { i in TextField("s\(["x", "y", "z"][i])", text: $size[i]) }
                            }
                            .onSubmit(apply)
                        }
                        Button("Use this site", action: apply)
                            .disabled(siteKind == .auto)
                    }
                    .padding(.top, 6)
                }

                Toggle("Prepare even if the structure is judged unusable", isOn: $session.allowUnusable)
                    .font(.callout)
            }
            .padding(.vertical, 4)
        }
    }

    private func apply() {
        switch siteKind {
        case .auto:
            session.site = .auto
        case .ligand:
            session.site = ligandName.isEmpty ? .auto : .ligand(ligandName.trimmingCharacters(in: .whitespaces))
        case .residues:
            let r = residues.split(whereSeparator: { $0 == " " || $0 == "," }).map(String.init)
            session.site = r.isEmpty ? .auto : .residues(r)
        case .box:
            let c = center.compactMap { Double($0) }, s = size.compactMap { Double($0) }
            session.site = (c.count == 3 && s.count == 3) ? .box(center: c, size: s) : .auto
        }
    }
}

// MARK: run

struct RunBar: View {
    @Environment(AppModel.self) private var model
    let session: QuickDockSession

    var body: some View {
        HStack(spacing: 12) {
            switch session.phase {
            case .assessing: ProgressView().controlSize(.small); Text("Checking the structure…")
            case .preparing: ProgressView().controlSize(.small); Text("Preparing receptor and ligand…")
            case .docking(let fraction, let message):
                ProgressView(value: fraction).frame(maxWidth: 240)
                Text(message.isEmpty ? "Docking…" : message).lineLimit(1)
            case .cancelled: Label("Cancelled", systemImage: "stop.circle").foregroundStyle(.secondary)
            case .finished: Label("Done", systemImage: "checkmark.circle").foregroundStyle(.green)
            default: EmptyView()
            }
            Spacer()
            if session.isBusy {
                Button("Cancel", role: .cancel) { session.cancel() }
            } else {
                Button("Dock") { model.run() }
                    .keyboardShortcut(.defaultAction)
                    .disabled(session.receptor == nil || session.ligand.trimmingCharacters(in: .whitespaces).isEmpty)
            }
        }
    }
}

// MARK: messages

struct MessagesView: View {
    @Environment(AppModel.self) private var model
    @Bindable var session: QuickDockSession
    @State private var showWarnings = false

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            if case .failed(let e) = session.phase { ErrorView(error: e, session: session) }
            ForEach(session.guardrails) { g in
                Label(g.message, systemImage: g.severity == .caution ? "exclamationmark.triangle" : "info.circle")
                    .foregroundStyle(g.severity == .caution ? Color.orange : Color.secondary)
            }
            if let q = session.prepared {
                Text(q.siteDescription).font(.callout).foregroundStyle(.secondary)
                Text(String(format: "Box %.0f × %.0f × %.0f Å · exhaustiveness %d · %d core(s)", q.box.size[0], q.box.size[1], q.box.size[2],
                            q.engine.exhaustiveness, q.engine.cores))
                    .font(.callout).foregroundStyle(.secondary)
            }
            if !session.warnings.isEmpty {
                DisclosureGroup("\(session.warnings.count) note(s) from preparation", isExpanded: $showWarnings) {
                    VStack(alignment: .leading, spacing: 4) {
                        ForEach(Array(session.warnings.enumerated()), id: \.offset) { _, w in
                            Text(w).font(.callout).textSelection(.enabled)
                        }
                    }
                }
            }
        }
    }
}

struct ErrorView: View {
    let error: UserFacingError
    @Bindable var session: QuickDockSession

    var body: some View {
        VStack(alignment: .leading, spacing: 6) {
            Label(error.title, systemImage: "exclamationmark.octagon.fill").foregroundStyle(.red).font(.headline)
            Text(error.detail).textSelection(.enabled)
            if !error.suggestion.isEmpty { Text(error.suggestion).foregroundStyle(.secondary) }
            switch error.action {
            case .prepareAnyway:
                Button("Prepare anyway") { session.allowUnusable = true; session.runAll() }
            case .copyDetails:
                Button("Copy details") {
                    NSPasteboard.general.clearContents()
                    NSPasteboard.general.setString("\(error.title)\n\(error.detail)\n\(error.log)", forType: .string)
                }
            default: EmptyView()
            }
        }
        .padding(12)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(RoundedRectangle(cornerRadius: 8).fill(Color.red.opacity(0.08)))
    }
}
#endif
