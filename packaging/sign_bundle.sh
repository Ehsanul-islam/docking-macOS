#!/usr/bin/env bash
# Sign every Mach-O in the app, inside out, then the app itself, with the hardened runtime.
#   SIGN_IDENTITY="Developer ID Application: Name (TEAMID)" packaging/sign_bundle.sh path/to/Dock.app
# SIGN_IDENTITY="-" signs ad hoc (local testing only; Gatekeeper will not accept it elsewhere).
set -euo pipefail
source "$(dirname "$0")/lib.sh"
APP="${1:?path to Dock.app}"; ID="${SIGN_IDENTITY:?set SIGN_IDENTITY}"
ENT="$ROOT/packaging/entitlements.plist"
ts=(--timestamp); [ "$ID" = "-" ] && ts=()          # no secure timestamp for ad-hoc signatures

n=0
while IFS= read -r -d '' f; do
  codesign --force --options runtime "${ts[@]}" --sign "$ID" "$f"; n=$((n + 1))
done < <(macho_files "$APP/Contents/Resources")
echo "signed $n nested binaries"
codesign --force --options runtime "${ts[@]}" --entitlements "$ENT" --sign "$ID" "$APP/Contents/MacOS/DockApp"
codesign --force --options runtime "${ts[@]}" --entitlements "$ENT" --sign "$ID" "$APP"
codesign --verify --deep --strict --verbose=2 "$APP"
