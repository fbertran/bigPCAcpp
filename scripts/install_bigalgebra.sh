#!/usr/bin/env bash
# Install the bigalgebra package and its CRAN dependencies from source tarballs.
# Requires that R is already installed (see scripts/install_r.sh).

set -euo pipefail

if ! command -v R >/dev/null 2>&1; then
  echo "R is not installed. Run scripts/install_r.sh first." >&2
  exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required to download CRAN packages." >&2
  exit 1
fi

CRAN_BASE="https://cran.r-project.org/src/contrib"

declare -A VERSIONS=(
  [BH]="1.87.0-1"
  [Rcpp]="1.1.0"
  [uuid]="1.2-1"
  [bigmemory.sri]="0.1.8"
  [bigmemory]="4.6.4"
)

order=(BH Rcpp uuid bigmemory.sri bigmemory bigalgebra)

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

for pkg in "${order[@]}"; do
  version="${VERSIONS[$pkg]}"
  tarball="${pkg}_${version}.tar.gz"
  url="${CRAN_BASE}/${tarball}"
  echo "Downloading ${pkg} ${version}..."
  curl -fsSL "$url" -o "$workdir/$tarball"
  echo "Installing ${pkg} ${version}..."
  R CMD INSTALL "$workdir/$tarball"
  echo
done

echo "All packages installed successfully."
