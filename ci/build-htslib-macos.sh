#!/usr/bin/env bash
set -euo pipefail

htslib_ref=${HTSLIB_REF:-develop}

git clone --depth 1 --branch "$htslib_ref" https://github.com/samtools/htslib.git htslib
(
  cd htslib
  git submodule update --init --recursive
  autoreconf -i
  ./configure --disable-libcurl --disable-shared --enable-static \
    --with-libdeflate \
    CPPFLAGS="-I$(brew --prefix libdeflate)/include" \
    LDFLAGS="-L$(brew --prefix libdeflate)/lib"
  make
  git rev-parse HEAD
)

for archive in \
  "$(brew --prefix libdeflate)/lib/libdeflate.a" \
  "$(brew --prefix bzip2)/lib/libbz2.a" \
  "$(brew --prefix xz)/lib/liblzma.a" \
  "$(brew --prefix zlib)/lib/libz.a"; do
  test -f "$archive" || { echo "Error: $archive not found" >&2; exit 1; }
done

ci/write-htslib-pc.sh htslib
