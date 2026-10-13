#!/usr/bin/env bash
set -euo pipefail

binary=${1:?binary path required}

deps=$(otool -L "$binary")
printf '%s\n' "$deps"

nonportable_deps=$(
  printf '%s\n' "$deps" |
    awk 'NR > 1 { print $1 }' |
    grep -E -v '^(/usr/lib/|/System/Library/)' || true
)

if [[ -n "$nonportable_deps" ]]; then
  echo "Error: non-portable dynamic dependencies linked on macOS:" >&2
  printf '%s\n' "$nonportable_deps" >&2
  exit 1
fi
