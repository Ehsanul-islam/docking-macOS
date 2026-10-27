# Packaging the app (plan §15)

Scripts live in `packaging/`. They run on macOS only and **have not been run yet**: everything below marked
*checked* was verified from the Linux container, everything else is for the first Mac build
(see `MAC_VERIFICATION_CHECKLIST.md` §3 and §8).

## Decisions

| Plan item | Decision | Why |
|---|---|---|
| D-2 Python packaging | python-build-standalone CPython 3.12, one bundle **per architecture** | A universal Python needs every `.so` merged with `lipo`; two single-architecture apps are simpler and each is testable. The engine can still be universal. |
| §15.2 [VERIFY] wheels | *checked 2026-10-27*: macOS wheels exist for arm64 **and** x86_64 for matplotlib, pillow, numpy, scipy, rdkit, gemmi, meeko, openmm, pdb2pqr, propka, dimorphite_dl (cp312). `pdbfixer` has no wheel (pure Python): built once from its sdist. | `pip download --platform macosx_14_0_{arm64,x86_64}` |
| Unused dependency | `prody` removed from `chemhelper/pyproject.toml` | Nothing imports it, and it has no macOS cp312 wheel. |
| Entitlements | hardened runtime, **none** extra | All nested code is re-signed with one identity, so library validation passes. Add one only when a hardened run proves it needed. |
| Channel | notarized DMG (plan §15.3). Mac App Store deferred. | The sandbox and bundled interpreters need separate investigation. |

### Intel differs: RDKit

PyPI has no Intel-Mac wheel of the newest RDKit, so `requirements-x86_64.txt` pins RDKit 2025.9.3 while arm64 uses
2026.9.1. Run on both versions here: the ligand, interaction and format tests pass on 2025.9.3 except that the
hypervalent SF5 test ligand (CCD A1JWR) cannot be embedded sanely; the app **refuses** it with a clear message
instead of docking a distorted molecule, and the test now accepts either outcome. Docking scores on the two
architectures may differ slightly for other reasons too: treat Intel as a separate release with its own redock run
(`tools/benchmark/redock.py`, checklist §6).

## Build, verify, ship (on a Mac)

```bash
# 0. once: pick the python-build-standalone 3.12 "install_only" archive for each arch from its releases page,
#    put URL and SHA-256 in packaging/python-build-standalone.lock, and export them:
export PBS_URL=… PBS_SHA256=…

engine/build_vina.sh arm64                       # build/engine/vina-arm64  (see checklist §3)
SIGN_IDENTITY="Developer ID Application: Name (TEAMID)" packaging/build_app.sh arm64
packaging/verify_bundle.sh build/app/arm64/Dock.app
NOTARY_PROFILE=dock-notary SIGN_IDENTITY="Developer ID Application: …" packaging/make_dmg.sh build/app/arm64/Dock.app
```

What each script does:

- `build_app.sh <arch>` – release build of `DockApp`; `Info.plist`; unpacks Python into `Resources/python`;
  installs the pinned wheels (`--only-binary`, `--no-deps`, then `pip check`); deletes pip's console scripts (their
  shebangs hold build paths) and test/doc trees; byte-compiles (the app sets `PYTHONDONTWRITEBYTECODE`, because
  writing into a signed bundle breaks the seal); copies Vina and the 3D viewer (`viewer/`, tests and source maps excluded); writes an SBOM (`sbom-python.txt`); signs.
- `sign_bundle.sh` – signs every Mach-O inside out with the hardened runtime and secure timestamp, then the app;
  `SIGN_IDENTITY=-` signs ad hoc for local runs.
- `verify_bundle.sh` – deep/strict signature check; no link outside the system or bundle (`otool -L`); no build paths;
  helper answers from inside the bundle with an empty environment and an unroutable proxy; every chemistry package
  imports; Vina runs; the signature is still valid after the helper ran. Run it on a clean second Mac as well.
- `make_dmg.sh` – DMG with an Applications link; with `NOTARY_PROFILE` also `notarytool submit --wait`, `stapler`, `spctl`.

At run time the app looks for `Resources/python/bin/python3` and `Resources/engine/vina` (`HelperLocator`, unit-tested),
sets `PYTHONHOME`, `PYTHONNOUSERSITE`, `PYTHONDONTWRITEBYTECODE`, and puts matplotlib's cache under `~/Library/Caches/Dock`.
The helper was also run here with an empty environment and a read-only home directory: it and the 2D image renderer work.

## Not verified yet

Everything that touches macOS: the script run itself, signing, notarization, Gatekeeper on a clean Mac, bundle size
(the wheels alone download as 86 MB arm64 / 96 MB x86_64; installed size is unmeasured), start-up time of the
helper from the bundle, whether the hardened runtime needs any entitlement, the product name and bundle identifier
(plan D-1; `org.example.dock` is a placeholder), the Sparkle update channel, and the licence audit of what is bundled
(plan §16; the SBOM is the input for it).
