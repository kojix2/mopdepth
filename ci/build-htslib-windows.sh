#!/usr/bin/env bash
set -euo pipefail

htslib_ref=${HTSLIB_REF:-develop}
mingw_prefix=${MINGW_PREFIX:?MINGW_PREFIX must be set by MSYS2}

git clone --depth 1 --branch "$htslib_ref" https://github.com/samtools/htslib.git htslib
(
  cd htslib
  git submodule update --init --recursive
  autoreconf -i
  ./configure --disable-libcurl --disable-shared --enable-static
  make -j2
  test -f libhts.a || { echo "libhts.a not found" >&2; ls -l; exit 1; }
  rm -f *.dll *.dll.a *.lib
  git rev-parse HEAD
)

for archive in \
  "$mingw_prefix/lib/libbz2.a" \
  "$mingw_prefix/lib/liblzma.a" \
  "$mingw_prefix/lib/libz.a" \
  "$mingw_prefix/lib/libdeflate.a" \
  "$mingw_prefix/lib/libregex.a" \
  "$mingw_prefix/lib/libtre.a" \
  "$mingw_prefix/lib/libintl.a"; do
  test -f "$archive" || { echo "Error: $archive not found" >&2; exit 1; }
done

ci/write-htslib-pc.sh htslib
