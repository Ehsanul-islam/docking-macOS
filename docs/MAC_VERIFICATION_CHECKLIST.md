# Mac verification checklist

Everything below could **not** be checked in the Linux development container. Run it on a real Mac
(Apple silicon first, then Intel if you ship x86_64) and send back the artifacts listed in §9.
Stop at the first failing step and report it; later steps depend on earlier ones.

## 0. What is and is not verified already

| Verified on Linux (CI + local) | Not verified anywhere |
|---|---|
| Python helper (`chemhelper`), all RPCs, preparation, docking runner, 2D images | Swift package build and tests on macOS |
| Vina 1.2.x engine behaviour via a Linux binary | macOS Vina build (arm64 / x86_64 / universal), signing, notarization |
| Redocking accuracy (Astex 85, PoseBusters 215) | Cross-architecture numerics (Mac vs Linux) |
| Fuzzing, cancellation, limits | Python helper packaged inside the app bundle |
| | The SwiftUI app itself (not written yet beyond `DockCore`) |

## 1. Prerequisites

- [ ] macOS 14+ (the Swift package targets `.macOS(.v14)`), Xcode 15+ with command-line tools: `xcodebuild -version`
- [ ] Homebrew; `brew install cmake boost python@3.11` (Python 3.10–3.12)
- [ ] Static Boost built **per architecture** (see §3). Record `uname -m` and `sw_vers`.

## 2. Swift package

```bash
swift build --package-path swift
swift test  --package-path swift
```

- [ ] Both succeed with 0 failures. CI (`macos-14`) runs the same `swift test`; local result should match.
- [ ] Record the Swift version (`swift --version`).

## 2b. App UI smoke test (the SwiftUI code has never been compiled; expect fixes)

`swift build --package-path swift` must compile `Sources/DockApp` (it is excluded from the Linux build). Fix any
compiler errors first; send them back verbatim. Then, with `DOCKCHEM_VINA` set and the Python helper installed:

```bash
swift run --package-path swift DockApp
```

- [ ] Window opens; Settings (⌘,) shows cores (default 2) and the engine status.
- [ ] Drop `pdb_structures_50/1STP.pdb`: an assessment appears ("Ready to dock", 121 residues, 2.6 Å).
- [ ] Ligand `CC1C2C(CCCCC(=O)O)SC1C(=O)N2`, preset Fast, **Dock**: progress bar moves, then a result table with one row.
- [ ] Select the row: the 2D interaction image appears; the pose stepper changes it.
- [ ] **Cancel** during a Thorough run: it stops within a few seconds; **Dock** works again.
- [ ] A missing/garbage file shows a red message with a readable explanation, not a crash.
- [ ] Drag a ligand file onto the ligand field; "New run" clears the window.
- [ ] Dark mode and keyboard navigation (Tab, Return = Dock) behave.

## 3. Build the engine (Vina)

```bash
export BOOST_ROOT_PREFIX=/path/to/static/boost   # per arch; see the script header
engine/build_vina.sh arm64        # then x86_64, then universal if shipping both
```

- [ ] Produces `build/engine/vina` and `build-manifest.txt`.
- [ ] `file build/engine/vina` shows the expected architecture(s); `otool -L build/engine/vina` lists **only** system libraries (no Homebrew/Boost dylibs).
- [ ] `build/engine/vina --version` prints 1.2.x and the manifest records the exact source tag/commit.
- [ ] Resolve every `[VERIFY]` note in `engine/build_vina.sh` (Makefile path `build/mac/release`, per-arch static Boost). Write the outcome of each into the manifest or this file.
- [ ] Engine is unmodified upstream Vina (diff of the source tree against the tag is empty).

## 4. Engine equivalence

```bash
python3 tools/equivalence/run_equivalence.py --help
```

Use `tools/equivalence/cases.example.json` as the template (receptor, ligand, centre, size, seed, exhaustiveness).

- [ ] Same architecture (a reference Vina built the same way vs the candidate): `--mode exact`, outputs byte-identical.
- [ ] Cross architecture (Mac arm64 vs the Linux binary): `--mode tolerance`, top score within the tolerance for every case. Use ≥10 cases including a flexible ligand (>10 torsions) and a metal-free pocket.
- [ ] Rosetta/x86_64 vs native arm64 if both are shipped: same check.

## 5. Python helper on the Mac

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install pytest ./chemhelper      # same as CI; dependencies come from chemhelper/pyproject.toml
pytest chemhelper -q
```

- [ ] Every dependency installs as a **wheel** on arm64 (RDKit, gemmi, OpenMM, PDBFixer, PDB2PQR, PROPKA, Meeko, scipy, matplotlib, Pillow). Note any that compile from source.
- [ ] `pytest chemhelper` passes. Expected differences must be explained, not skipped.
- [ ] OpenMM uses the CPU platform and honours the cores option (`system.resources` RPC reports the same default of 2).

## 6. End-to-end with the real Mac engine

```bash
export DOCKCHEM_VINA=$PWD/build/engine/vina
python3 -m dockchem.server <<'EOF'   # one-shot smoke test; the app talks to it over stdio
{"jsonrpc":"2.0","id":1,"method":"dock.engine","params":{}}
EOF
```

- [ ] `dock.engine` reports the arm64 binary and its version.
- [ ] Run the PoseBusters/Astex sample on the Mac and compare with the Linux numbers in `docs/RELIABILITY.md`:
  `python3 tools/benchmark/redock.py --set astex --limit 20 --workers 2 --vina "$DOCKCHEM_VINA" --out /tmp/mac_redock --zip <dataset>`
- [ ] Success rate (RMSD ≤ 2 Å, top-1) within ±10 percentage points of the Linux run on the same subset; any structure whose preparation status differs is listed.
- [ ] `dock.run` with 1 ligand uses all configured cores; with many ligands, one core each (check Activity Monitor).
- [ ] Cancel mid-dock: the Vina process disappears within a few seconds and no orphan remains (`pgrep vina`).

## 7. Offline and sandbox behaviour

- [ ] Disconnect from the network (or set an unroutable proxy): preparation, docking and image rendering all still work. Nothing may need the internet at run time.
- [ ] If the app will be sandboxed: the helper and engine launch from inside the bundle, read/write only user-chosen folders and the project directory.

## 8. Packaging, signing, notarization (plan D-2)

Scripts and decisions: `docs/PACKAGING.md` (`packaging/build_app.sh`, `verify_bundle.sh`, `make_dmg.sh`). `verify_bundle.sh` covers the first four boxes below; its output is what to send back.

- [ ] The helper runs from the bundle with a **relocatable** Python (no dependence on `/usr/local`, Homebrew or the build venv): copy the `.app` to a clean user account / second Mac and run §6.
- [ ] Every Mach-O in the bundle (Python, `.so`, `.dylib`, `vina`) is signed with the hardened runtime: `codesign --verify --deep --strict --verbose=2 Dock.app`.
- [ ] `xcrun notarytool submit … --wait` returns *Accepted*; `spctl -a -vv Dock.app` says *accepted, source=Notarized Developer ID*; staple with `xcrun stapler staple`.
- [ ] First launch on a clean Mac shows no Gatekeeper warning.

## 9. What to send back

- `sw_vers`, `uname -m`, `swift --version`, `python3 --version`, `pip freeze`
- Output of §2 (`swift test`), §3 (`build-manifest.txt`, `otool -L`, `file`), §4 (equivalence report), §5 (`pytest` summary), §6 (redock `summary.json`/`results.csv`)
- Notarization log (`notarytool log`) if §8 fails
- For any failure: the exact command, full output, and the first step that failed

## 10. Pass criteria

All boxes in §2–§6 ticked, §7 offline works, §8 notarized and runs on a clean machine. Anything else is a defect to fix, not a note to file away. Items that cannot be satisfied (for example no Intel Mac available) must be written down as **not verified**, not omitted.
