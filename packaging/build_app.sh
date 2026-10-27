#!/usr/bin/env bash
# Build Dock.app for one architecture: SwiftUI executable + relocatable Python with the chemistry stack + Vina.
#
#   packaging/build_app.sh arm64          (or x86_64)
#
# Required environment
#   PBS_URL, PBS_SHA256   python-build-standalone CPython 3.12 "install_only" archive for that architecture and its
#                         SHA-256 (pick both from the python-build-standalone releases page and keep them in
#                         packaging/python-build-standalone.lock so builds are reproducible)
# Optional
#   VINA_BIN       engine binary (default build/engine/vina-<arch>, from engine/build_vina.sh)
#   SIGN_IDENTITY  "Developer ID Application: Name (TEAMID)", or "-" for ad-hoc (default: "-")
#   VERSION, BUILD bundle versions (default 0.1.0 / git commit count)
#   BUNDLE_ID      default org.example.dock  (set the real one once the product name is decided, plan D-1)
set -euo pipefail
source "$(dirname "$0")/lib.sh"
ARCH="${1:-}"; check_arch "$ARCH"
[ "$(uname -s)" = Darwin ] || die "run on macOS"
: "${PBS_URL:?see the header of this script}" "${PBS_SHA256:?see the header of this script}"
VERSION="${VERSION:-0.1.0}"; BUILD="${BUILD:-$(git -C "$ROOT" rev-list --count HEAD)}"
BUNDLE_ID="${BUNDLE_ID:-org.example.dock}"
VINA_BIN="${VINA_BIN:-$ROOT/build/engine/vina-$ARCH}"
[ -x "$VINA_BIN" ] || die "engine binary not found at $VINA_BIN (build it with engine/build_vina.sh $ARCH)"

OUT="$ROOT/build/app/$ARCH"; APP="$OUT/Dock.app"; RES="$APP/Contents/Resources"
rm -rf "$OUT"; mkdir -p "$APP/Contents/MacOS" "$RES/engine" "$OUT/dl"

echo "==> 1/6 SwiftUI executable ($ARCH)"
swift build -c release --arch "$ARCH" --product DockApp --package-path "$ROOT/swift"
BIN="$(swift build -c release --arch "$ARCH" --package-path "$ROOT/swift" --show-bin-path)/DockApp"
cp "$BIN" "$APP/Contents/MacOS/DockApp"

echo "==> 2/6 Info.plist"
sed -e "s|@BUNDLE_ID@|$BUNDLE_ID|" -e "s|@VERSION@|$VERSION|" -e "s|@BUILD@|$BUILD|" \
    "$ROOT/packaging/Info.plist.template" > "$APP/Contents/Info.plist"
plutil -lint "$APP/Contents/Info.plist"

echo "==> 3/6 Python (python-build-standalone)"
curl -fL --retry 3 -o "$OUT/dl/python.tar.gz" "$PBS_URL"
echo "$PBS_SHA256  $OUT/dl/python.tar.gz" | shasum -a 256 -c -
tar -xzf "$OUT/dl/python.tar.gz" -C "$RES"          # unpacks to Resources/python
PY="$RES/python/bin/python3"
as_arch "$ARCH" "$PY" -c 'import sys; assert sys.version_info[:2] == (3, 12), sys.version'

echo "==> 4/6 Chemistry packages (pinned, wheels only; pdbfixer is pure Python and built once)"
LOCK="$ROOT/packaging/requirements-$ARCH.txt"
grep -v '^pdbfixer' "$LOCK" | grep -v '^#' > "$OUT/dl/wheels.txt"
as_arch "$ARCH" "$PY" -m pip install --no-cache-dir --no-deps --only-binary=:all: -r "$OUT/dl/wheels.txt"
as_arch "$ARCH" "$PY" -m pip wheel --no-deps --no-cache-dir -w "$OUT/dl" "$(grep '^pdbfixer' "$LOCK" | sed 's/ .*//')"
as_arch "$ARCH" "$PY" -m pip install --no-cache-dir --no-deps "$OUT"/dl/pdbfixer-*.whl
as_arch "$ARCH" "$PY" -m pip install --no-cache-dir --no-deps "$ROOT/chemhelper"
as_arch "$ARCH" "$PY" -m pip check
as_arch "$ARCH" "$PY" -m pip freeze > "$OUT/sbom-python.txt"

echo "==> 5/6 Make the bundle sealable"
# pip's console scripts carry the build machine's absolute path in their shebang: drop everything but python itself.
find "$RES/python/bin" -type f ! -name 'python3*' ! -name 'python' -delete
rm -rf "$RES/python/lib/python3.12/test" "$RES/python/lib/python3.12/idlelib" "$RES/python/lib/python3.12/tkinter" \
       "$RES/python/lib/python3.12/turtledemo" "$RES/python/share" "$RES/python/include"
find "$RES/python" -type d \( -name tests -o -name test -o -name testing \) -path '*site-packages*' -prune -exec rm -rf {} + 2>/dev/null || true
# byte-compile now: Python must not write .pyc files into a signed bundle at run time (the app sets
# PYTHONDONTWRITEBYTECODE for the same reason)
as_arch "$ARCH" "$PY" -m compileall -q -j 0 "$RES/python/lib" || true
cp "$VINA_BIN" "$RES/engine/vina"; chmod +x "$RES/engine/vina"
lipo -archs "$RES/engine/vina" | grep -qw "$ARCH" || die "engine binary is not built for $ARCH"
cp "$OUT/sbom-python.txt" "$RES/sbom-python.txt"
du -sh "$APP" "$RES/python" "$RES/engine" | sed 's/^/size: /'

echo "==> 6/6 Sign"
SIGN_IDENTITY="${SIGN_IDENTITY:--}" "$ROOT/packaging/sign_bundle.sh" "$APP"
echo "built $APP"
