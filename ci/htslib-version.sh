#!/usr/bin/env sh
set -eu

version=$(sed -n 's/^VERSION=//p' "${1:-htslib}/version.sh")
case "$version" in
  *.*.*) ;;
  *.*) version="$version.0" ;;
  *) version="$version.0.0" ;;
esac

printf '%s\n' "$version"
