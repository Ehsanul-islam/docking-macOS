# docking-macOS

Native macOS molecular docking workstation built around unmodified AutoDock Vina 1.2.7.
See [docs/vina_macos_plan_v2.md](docs/vina_macos_plan_v2.md) for the full plan.

## Layout

| Path | Contents | Plan |
|---|---|---|
| `chemhelper/` | Python chemistry helper: JSON-RPC server, `.dockproj` store | §4.4, §4.6, §14 |
| `engine/` | Vina 1.2.7 build script (universal binary) | §15.1 |
| `tools/equivalence/` | Engine-equivalence harness | §12.1 |
| `swift/` | `DockCore` Swift package: RPC types, manifest, engine runner | §4 |
| `schemas/` | JSON schema for `project.json` | §4.6 |

## Status

Phase 0 (foundations) scaffold. Python parts are tested on Linux in CI. The Swift package,
`engine/build_vina.sh` and the equivalence run against real Vina binaries have **not yet been
run**: they need macOS. See docs/PHASE0_STATUS.md.

## Development

```
python -m venv .venv && . .venv/bin/activate
pip install pytest -e chemhelper
pytest chemhelper tools/equivalence
```
