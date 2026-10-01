# Installation

`kira-spliceqc` is a single static binary; HDF5 (for `.h5ad` input) is linked
in, so nothing else is needed at runtime.

## Binaries

Every tagged release ships archives for Linux (x86_64, aarch64), macOS (arm64,
x86_64) and Windows (x86_64) on the
[releases page](https://github.com/ARyaskov/kira-spliceqc/releases), each with
a SHA-256 sum, the documentation and the default geneset catalog.

```bash
tar -xzf kira-spliceqc-<version>-x86_64-unknown-linux-gnu.tar.gz
sudo install kira-spliceqc-<version>-x86_64-unknown-linux-gnu/kira-spliceqc /usr/local/bin/
```

## bioconda

```bash
conda install -c conda-forge -c bioconda kira-spliceqc
```

(The recipe lives in `packaging/bioconda/` and is submitted with the first
tagged release; until it is merged, use a binary or `cargo install`.)

## Container

```bash
docker build -f packaging/docker/Dockerfile -t kira-spliceqc .
docker run --rm -v "$PWD/data:/data" kira-spliceqc run --input /data/pbmc3k --out /data/out
```

BioContainers images (`biocontainers/kira-spliceqc`) are built automatically
from the bioconda package.

## From source

Rust 1.95 or newer and CMake (HDF5 is built from source once):

```bash
cargo install --locked kira-spliceqc
```

The repository's `.cargo/config.toml` builds for the host CPU; set
`RUSTFLAGS="-C target-cpu=x86-64-v2"` for a portable binary.

## Python

```bash
pip install kira-spliceqc        # wrapper only; needs the binary on PATH or KIRA_SPLICEQC_BIN
```

See `python/README.md`: `kira_spliceqc.run(...)` returns an `AnnData` with the
per-cell columns in `obs`.
