#if canImport(SwiftUI)
import SwiftUI
import AppKit
import WebKit
import DockCore

/// The 3D viewer: a WKWebView showing viewer/index.html (Mol*), driven by a `ViewerController`.
struct MoleculeView: NSViewRepresentable {
    let controller: ViewerController
    let directory: URL

    func makeCoordinator() -> Coordinator { Coordinator(controller: controller, directory: directory) }

    func makeNSView(context: Context) -> WKWebView {
        let config = WKWebViewConfiguration()
        config.userContentController.add(context.coordinator, name: "dock")
        let web = WKWebView(frame: .zero, configuration: config)
        web.setValue(false, forKey: "drawsBackground")            // the page paints its own background
        web.navigationDelegate = context.coordinator
        #if DEBUG
        web.isInspectable = true
        #endif
        context.coordinator.web = web
        controller.attach(context.coordinator)
        web.loadFileURL(directory.appendingPathComponent("index.html"), allowingReadAccessTo: directory)
        return web
    }

    func updateNSView(_ nsView: WKWebView, context: Context) {}

    static func dismantleNSView(_ nsView: WKWebView, coordinator: Coordinator) {
        nsView.configuration.userContentController.removeScriptMessageHandler(forName: "dock")
    }

    @MainActor
    final class Coordinator: NSObject, WKScriptMessageHandler, WKNavigationDelegate, ViewerTransport {
        let controller: ViewerController
        let directory: URL
        weak var web: WKWebView?

        init(controller: ViewerController, directory: URL) { self.controller = controller; self.directory = directory }

        // page -> app
        func userContentController(_ ucc: WKUserContentController, didReceive message: WKScriptMessage) {
            guard JSONSerialization.isValidJSONObject(message.body),
                  let data = try? JSONSerialization.data(withJSONObject: message.body),
                  let json = try? JSONDecoder().decode(JSONValue.self, from: data),
                  let event = ViewerEvent(json) else { return }
            controller.handle(event)
        }

        // app -> page
        func send(_ command: ViewerCommand) {
            web?.evaluateJavaScript("window.dock.command(\(command.json));") { _, error in
                if let error { NSLog("viewer command \(command.name) failed: \(error.localizedDescription)") }
            }
        }

        // the viewer page is local: nothing else may be loaded into it
        func webView(_ webView: WKWebView, decidePolicyFor action: WKNavigationAction,
                     decisionHandler: @escaping (WKNavigationActionPolicy) -> Void) {
            if let url = action.request.url, url.isFileURL, url.path.hasPrefix(directory.path) { decisionHandler(.allow) }
            else { decisionHandler(.cancel) }
        }

        // the web content process can be killed under memory pressure: reload and let the controller resend
        func webViewWebContentProcessDidTerminate(_ webView: WKWebView) {
            controller.pageReloaded()
            webView.reload()
        }
    }
}

/// The viewer with its fallbacks: not installed, or the engine reported an error.
struct ViewerPanel: View {
    @Environment(AppModel.self) private var model

    var body: some View {
        ZStack(alignment: .bottom) {
            if let dir = model.viewerDirectory {
                MoleculeView(controller: model.viewer, directory: dir)
            } else {
                ContentUnavailableView("3D view unavailable", systemImage: "cube.transparent",
                                       description: Text("The viewer files are missing from the app."))
            }
            if let error = model.viewer.lastError {
                Label(error, systemImage: "exclamationmark.triangle.fill")
                    .font(.callout)
                    .padding(10)
                    .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 10))
                    .padding(.bottom, 70)
            }
        }
        .onAppear {
            model.viewer.onScreenshot = { data in
                let panel = NSSavePanel()
                panel.allowedContentTypes = [.png]
                panel.nameFieldStringValue = "dock-view.png"
                if panel.runModal() == .OK, let url = panel.url { try? data.write(to: url) }
            }
        }
    }
}
#endif
