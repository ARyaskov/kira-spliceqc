#!/bin/bash
# Bioconda build script: a plain cargo install with the lockfile. The
# repository's .cargo/config.toml targets the build host's CPU, which would
# produce non-portable binaries; override it with the conda baseline.
set -euxo pipefail
export RUSTFLAGS="${RUSTFLAGS:-} -C target-cpu=x86-64-v2"
case "$(uname -m)" in
  aarch64|arm64) export RUSTFLAGS="${RUSTFLAGS/-C target-cpu=x86-64-v2/}" ;;
esac
cargo install --locked --no-track --root "$PREFIX" --path .
mkdir -p "$PREFIX/share/kira-spliceqc"
cp -r resources "$PREFIX/share/kira-spliceqc/"
