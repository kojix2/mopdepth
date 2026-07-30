#!/usr/bin/env sh
set -eu

htslib_dir=${1:-htslib}
htslib_abs=$(cd "$htslib_dir" && pwd -P)
version=$(ci/htslib-version.sh "$htslib_dir")

rm -f "$htslib_abs/htslib-uninstalled.pc"

pc_abs=$htslib_abs
if command -v cygpath >/dev/null 2>&1; then
  pc_abs=$(cygpath -m "$htslib_abs")
fi

cat > "$htslib_abs/htslib.pc" <<EOF
includedir=$pc_abs
libdir=$pc_abs

Name: htslib
Description: C library for high-throughput sequencing data formats
Version: $version
Cflags: -I\${includedir}
Libs: -L\${libdir} -lhts
EOF

PKG_CONFIG_DISABLE_UNINSTALLED=1 PKG_CONFIG_PATH="$htslib_abs:${PKG_CONFIG_PATH:-}" pkg-config --modversion htslib
