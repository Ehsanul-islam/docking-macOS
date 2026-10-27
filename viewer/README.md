# 3D viewer (Mol* shell)

`index.html` + `app.css` + `app.js` wrap **Mol* 5.13.1** (MIT, `vendor/`, licence in `vendor/MOLSTAR-LICENSE`) in a styled,
dark/light interface. The app hosts it in a `WKWebView` and drives it with small JSON commands; nothing is loaded from the
network.

The scene is described with MolViewSpec and rebuilt from state on every change (camera preserved): protein cartoon
coloured N to C, pocket residues and optional surface, the docked pose with bond orders, interactions as dashed lines
(same records as the 2D diagram, from `complex.interactions`), the search box, residue labels.

## Commands (`window.dock.command(json)`)

| `cmd` | fields |
|---|---|
| `setTheme` | `theme`: `dark` or `light` |
| `loadReceptor` | `data` (base64 PDB text), `format` (`pdb`), `name` |
| `loadPose` | `data` (base64 SDF; a multi-pose file is split), `index` (1-based), `score`, `poses` (`[{index, score}]`) |
| `setInteractions` | `interactions` (records from `complex.interactions`) |
| `setBox` | `center`, `size` (or null to hide) |
| `setOptions` | any of `cartoon surface pocket interactions labels box spin` (booleans) |
| `focus` | `target`: `ligand` or `all` |
| `screenshot` | posts a `screenshot` event with a PNG data URI |
| `clear` | empties the view |

## Events (to `webkit.messageHandlers.dock`, or a `dock-event` DOM event)

`ready`, `rendered`, `poseSelected {index}`, `screenshot {dataUri}`, `warning {message}`, `error {message}`.

## Test

```bash
NODE_PATH=$(npm root -g) node viewer/test/shot.js /tmp/shots dark     # also: light
```

Loads real docking data (`test/data`, from 1STP) in headless Chromium, drives every command, writes screenshots and
exits non-zero on any console error. Examples: `docs/examples/viewer_dark.png`, `viewer_light.png`.
