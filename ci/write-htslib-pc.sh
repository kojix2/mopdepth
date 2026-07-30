#!/usr/bin/env sh
set -eu

htslib_dir=${1:-htslib}
htslib_abs=$(cd "$htslib_dir" && pwd -P)
version=$(ci/htslib-version.sh "$htslib_dir")

cat > "$htslib_abs/htslib.pc" <<EOF
includedir=$htslib_abs
libdir=$htslib_abs

Name: htslib
Description: C library for high-throughput sequencing data formats
Version: $version
Cflags: -I\${includedir}
Libs: -L\${libdir} -lhts
EOF

PKG_CONFIG_PATH="$htslib_abs:${PKG_CONFIG_PATH:-}" pkg-config --modversion htslib
