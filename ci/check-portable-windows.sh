#!/usr/bin/env bash
set -euo pipefail

binary=${1:?binary path required}

objdump -p "$binary" | sed -n '/DLL Name/p' | tee dlls.txt
if grep -Ei "hts|deflate|bz2|lzma|zlib|libgcc|libstdc\\+\\+|libwinpthread|libpcre|libtre|libregex" dlls.txt; then
  echo "Error: non-portable DLL dependency linked on Windows" >&2
  exit 1
fi
