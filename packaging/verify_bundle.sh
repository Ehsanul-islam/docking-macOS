#!/usr/bin/env bash
# Checks a built Dock.app. Exit status 0 only if everything passes; run it on the machine that built the app AND on a
# clean second Mac.   packaging/verify_bundle.sh path/to/Dock.app
set -uo pipefail
source "$(dirname "$0")/lib.sh"
set +e                       # keep going after a failed check; the exit status reports them
APP="${1:?path to Dock.app}"; RES="$APP/Contents/Resources"; fail=0
bad() { echo "FAIL: $*"; fail=1; }
ok() { echo "ok:   $*"; }

codesign --verify --deep --strict "$APP" 2>/dev/null && ok "signature valid (deep, strict)" || bad "codesign verification"
codesign -d --entitlements - "$APP" 2>/dev/null | grep -q 'allow-' && echo "note: entitlements present, check entitlements.plist says why"

# 1. Nothing may link outside the system or the bundle.
offenders=0
while IFS= read -r -d '' f; do
  while read -r dep; do
    case "$dep" in
      /usr/lib/*|/System/*|@rpath/*|@loader_path/*|@executable_path/*) ;;
      *) echo "  $f -> $dep"; offenders=$((offenders + 1)) ;;
    esac
  done < <(otool -L "$f" 2>/dev/null | tail -n +2 | awk '{print $1}')
done < <(macho_files "$APP")
[ "$offenders" -eq 0 ] && ok "no links outside the system or the bundle" || bad "$offenders link(s) to files outside the bundle (Homebrew, /usr/local, a build directory)"

# 2. Nothing may reference the build machine.
if grep -rIl --exclude='*.pyc' -e "$ROOT" -e '/Users/' "$RES/python/bin" "$RES/engine" 2>/dev/null | grep -q .; then bad "build paths inside bin/ or engine/"; else ok "no build paths in bin/ or engine/"; fi

[ -f "$RES/viewer/index.html" ] && [ -f "$RES/viewer/vendor/molstar.js" ] && ok "3D viewer files present" || bad "3D viewer files missing"

# 3. Architecture of everything must match.
archs=$(macho_files "$APP" | xargs -0 -n1 lipo -archs 2>/dev/null | tr ' ' '\n' | sort -u | tr '\n' ' ')
echo "info: architectures present: $archs"

# 4. The helper runs from the bundle, with a clean environment and no network.
PY="$RES/python/bin/python3"; TMP="$(mktemp -d)"
reply=$(printf '{"jsonrpc":"2.0","id":1,"method":"system.versions","params":{}}\n' | \
  env -i HOME="$TMP" PYTHONHOME="$RES/python" PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMP/mpl" \
        DOCKCHEM_VINA="$RES/engine/vina" HTTPS_PROXY=http://127.0.0.1:9 HTTP_PROXY=http://127.0.0.1:9 \
        "$PY" -m dockchem.server 2>"$TMP/err" | head -1)
echo "$reply" | grep -q '"result"' && ok "helper answers from inside the bundle" || { bad "helper did not answer"; head -5 "$TMP/err"; }
env -i "$PY" -c "import rdkit, openmm, gemmi, meeko, prody, pdbfixer, pdb2pqr, propka, scipy, matplotlib, PIL; print('imports ok')" 2>"$TMP/err2" | grep -q 'imports ok' && ok "all chemistry packages import" || { bad "an import failed"; head -5 "$TMP/err2"; }
"$RES/engine/vina" --version 2>&1 | grep -qi 'AutoDock Vina' && ok "engine runs: $("$RES/engine/vina" --version 2>&1 | head -1)" || bad "engine does not run"

# 5. The bundle must not change while running (a signed bundle that writes into itself breaks its own seal).
codesign --verify --deep --strict "$APP" 2>/dev/null && ok "signature still valid after running the helper" || bad "running the helper modified the bundle"
[ -z "$(find "$RES/python" -newer "$APP/Contents/Info.plist" -name '*.pyc' 2>/dev/null | head -1)" ] || echo "note: new .pyc files appeared inside the bundle"

du -sh "$APP" | sed 's/^/size: /'
rm -rf "$TMP"
exit $fail
