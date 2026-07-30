#!/usr/bin/env sh
set -eu

htslib_dir=${1:-htslib}

if git -C "$htslib_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  version=$(git -C "$htslib_dir" show HEAD:version.sh | sed -n 's/^VERSION=//p')
else
  version=$(sed -n 's/^VERSION=//p' "$htslib_dir/version.sh")
fi
case "$version" in
  *.*.*) ;;
  *.*) version="$version.0" ;;
  *) version="$version.0.0" ;;
esac

printf '%s\n' "$version"
