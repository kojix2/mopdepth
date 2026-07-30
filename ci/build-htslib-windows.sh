#!/usr/bin/env bash
set -euo pipefail

htslib_ref=${HTSLIB_REF:-develop}

git clone --depth 1 --branch "$htslib_ref" https://github.com/samtools/htslib.git htslib
(
  cd htslib
  git submodule update --init --recursive
  autoreconf -i
  ./configure --disable-libcurl --disable-shared --enable-static
  make -j2
  test -f libhts.a || { echo "libhts.a not found" >&2; ls -l; exit 1; }
  rm -f *.dll *.dll.a *.lib || true
  git rev-parse HEAD
)

for archive in \
  /mingw64/lib/libbz2.a \
  /mingw64/lib/liblzma.a \
  /mingw64/lib/libz.a \
  /mingw64/lib/libdeflate.a \
  /mingw64/lib/libregex.a \
  /mingw64/lib/libtre.a \
  /mingw64/lib/libintl.a; do
  test -f "$archive" || { echo "Error: $archive not found" >&2; exit 1; }
done

ci/write-htslib-pc.sh htslib
