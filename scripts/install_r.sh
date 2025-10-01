#!/usr/bin/env bash
# This script installs a minimal R environment inside the container.
# It updates the package index, installs the base R runtime and development tools,
# and cleans up cached package data to keep the image lean.

set -euo pipefail

# Ensure non-interactive apt operations.
export DEBIAN_FRONTEND=noninteractive

# Update package lists and install R base packages.
apt-get update
apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    r-recommended

# Remove apt cache to reduce layer size.
apt-get clean
rm -rf /var/lib/apt/lists/*

echo "R installation complete. You can now run 'R' or 'Rscript'."
