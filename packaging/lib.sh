# Shared helpers for the packaging scripts (sourced, not run).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() { echo "error: $*" >&2; exit 1; }

check_arch() {
  case "${1:-}" in arm64|x86_64) ;; *) die "usage: $0 <arm64|x86_64>" ;; esac
}

# Run a command as the target architecture (Rosetta for x86_64 on Apple silicon).
as_arch() {
  local arch="$1"; shift
  if [ "$(uname -m)" = "$arch" ]; then "$@"; else arch "-$arch" "$@"; fi
}

# Every Mach-O file under a directory (executables, .so, .dylib), NUL separated.
macho_files() {
  find "$1" -type f ! -name '*.py' ! -name '*.pyc' -print0 | while IFS= read -r -d '' f; do
    if file -b "$f" | grep -q 'Mach-O'; then printf '%s\0' "$f"; fi
  done
}
