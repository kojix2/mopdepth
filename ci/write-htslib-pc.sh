#!/usr/bin/env sh
set -eu

htslib_dir=${1:-htslib}
version=$(ci/htslib-version.sh "$htslib_dir")
prefix=$(pwd)

sed \
  -e "s#@-includedir@#$prefix/$htslib_dir#g" \
  -e "s#@-libdir@#$prefix/$htslib_dir#g" \
  -e "s#@-PACKAGE_VERSION@#$version#g" \
  -e "s#@.PACKAGE_VERSION@#$version#g" \
  "$htslib_dir/htslib.pc.tmp" > "$htslib_dir/htslib.pc"

PKG_CONFIG_PATH="$prefix/$htslib_dir:${PKG_CONFIG_PATH:-}" pkg-config --modversion htslib
