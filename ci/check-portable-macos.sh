#!/usr/bin/env bash
set -euo pipefail

binary=${1:?binary path required}

deps=$(otool -L "$binary")
printf '%s\n' "$deps"
if printf '%s\n' "$deps" | grep -E -q "lib(hts|deflate|bz2|lzma|z).*dylib|/(opt/homebrew|usr/local)/"; then
  echo "Error: non-portable Homebrew or htslib dylib dependency linked on macOS" >&2
  exit 1
fi
