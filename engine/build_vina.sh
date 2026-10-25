#!/usr/bin/env bash
# Build AutoDock Vina 1.2.7 (unmodified upstream) for one architecture, or a universal binary.
# Usage: engine/build_vina.sh [arm64|x86_64|universal]   (run on macOS)
# [VERIFY] upstream CLI build is the Makefile in build/mac/release; plan §15.1 mentions CMake.
# [VERIFY] static Boost is built per-arch (Homebrew Boost is single-arch).
set -euo pipefail

VINA_TAG="v1.2.7"
ARCH="${1:-universal}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/build/engine"
SRC="$OUT/src"
mkdir -p "$OUT"

[ -d "$SRC" ] || git clone --depth 1 --branch "$VINA_TAG" https://github.com/ccsb-scripps/AutoDock-Vina "$SRC"

build_one() {
  local arch="$1" boost="${BOOST_ROOT_PREFIX:?set BOOST_ROOT_PREFIX to a static Boost prefix, one per arch}/$arch"
  local flags="-arch $arch -mmacosx-version-min=14.0 -O3"
  ( cd "$SRC/build/mac/release" && make clean >/dev/null 2>&1 || true
    make -j"$(sysctl -n hw.ncpu)" BASE="$boost" \
         C_OPTIONS="$flags" GPP="clang++ -arch $arch" >/dev/null )
  cp "$SRC/build/mac/release/vina" "$OUT/vina-$arch"
  echo "$arch: clang++ $(clang++ --version | head -1) | flags: $flags" >> "$OUT/build-manifest.txt"
}

: > "$OUT/build-manifest.txt"
echo "vina $VINA_TAG commit $(git -C "$SRC" rev-parse HEAD)" >> "$OUT/build-manifest.txt"
case "$ARCH" in
  arm64|x86_64) build_one "$ARCH"; cp "$OUT/vina-$ARCH" "$OUT/vina" ;;
  universal) build_one arm64; build_one x86_64
             lipo -create "$OUT/vina-arm64" "$OUT/vina-x86_64" -output "$OUT/vina" ;;
  *) echo "unknown arch $ARCH" >&2; exit 2 ;;
esac
shasum -a 256 "$OUT"/vina* >> "$OUT/build-manifest.txt"
echo "built $OUT/vina"
