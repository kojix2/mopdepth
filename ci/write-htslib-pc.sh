#!/usr/bin/env sh
set -eu

htslib_dir=${1:-htslib}
htslib_abs=$(cd "$htslib_dir" && pwd -P)
version=$(ci/htslib-version.sh "$htslib_dir")

rm -f "$htslib_abs/htslib-uninstalled.pc"

cat > "$htslib_abs/htslib.pc" <<EOF
includedir=$htslib_abs
libdir=$htslib_abs

Name: htslib
Description: C library for high-throughput sequencing data formats
Version: $version
Cflags: -I\${includedir}
Libs: -L\${libdir} -lhts
EOF

PKG_CONFIG_DISABLE_UNINSTALLED=1 PKG_CONFIG_PATH="$htslib_abs:${PKG_CONFIG_PATH:-}" pkg-config --modversion htslib
pc_path=$(PKG_CONFIG_DISABLE_UNINSTALLED=1 PKG_CONFIG_PATH="$htslib_abs:${PKG_CONFIG_PATH:-}" pkg-config --path htslib)
test "$pc_path" = "$htslib_abs/htslib.pc" || {
  echo "Error: pkg-config resolved $pc_path, expected $htslib_abs/htslib.pc" >&2
  exit 1
}
