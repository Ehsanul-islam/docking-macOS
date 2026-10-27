#!/usr/bin/env bash
# Disk image with the app and an Applications shortcut, then (optionally) notarize and staple it.
#   packaging/make_dmg.sh path/to/Dock.app            builds build/Dock-<version>.dmg
#   NOTARY_PROFILE=<keychain profile> packaging/make_dmg.sh path/to/Dock.app    also notarizes + staples
# Create the profile once with:  xcrun notarytool store-credentials <name> --apple-id … --team-id … --password <app-specific>
set -euo pipefail
source "$(dirname "$0")/lib.sh"
APP="${1:?path to Dock.app}"
VER="$(/usr/libexec/PlistBuddy -c 'Print :CFBundleShortVersionString' "$APP/Contents/Info.plist")"
ARCH="$(lipo -archs "$APP/Contents/MacOS/DockApp" | tr ' ' '-')"
DMG="$ROOT/build/Dock-$VER-$ARCH.dmg"; STAGE="$(mktemp -d)"
cp -R "$APP" "$STAGE/"; ln -s /Applications "$STAGE/Applications"
rm -f "$DMG"; hdiutil create -volname "Dock $VER" -srcfolder "$STAGE" -ov -format UDZO "$DMG"
rm -rf "$STAGE"
if [ -n "${SIGN_IDENTITY:-}" ] && [ "$SIGN_IDENTITY" != "-" ]; then codesign --force --timestamp --sign "$SIGN_IDENTITY" "$DMG"; fi
if [ -n "${NOTARY_PROFILE:-}" ]; then
  xcrun notarytool submit "$DMG" --keychain-profile "$NOTARY_PROFILE" --wait
  xcrun stapler staple "$DMG"
  xcrun stapler validate "$DMG"
  spctl -a -t open --context context:primary-signature -vv "$DMG"
fi
echo "$DMG"
