#!/usr/bin/env bash
set -euo pipefail

binary=${1:?binary path required}

dlls=$(objdump -p "$binary" | sed -n '/DLL Name/p')
printf '%s\n' "$dlls"
if printf '%s\n' "$dlls" | grep -Ei "hts|deflate|bz2|lzma|zlib|libgcc|libstdc\\+\\+|libwinpthread|libpcre|libtre|libregex"; then
  echo "Error: non-portable DLL dependency linked on Windows" >&2
  exit 1
fi
