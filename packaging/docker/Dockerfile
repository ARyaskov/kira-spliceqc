# kira-spliceqc container (BioContainers-style: one tool, one binary).
#
#   docker build -f packaging/docker/Dockerfile -t kira-spliceqc .
#   docker run --rm -v $PWD/data:/data kira-spliceqc run --input /data/pbmc3k --out /data/out
#
# HDF5 is built from source and linked statically (hdf5-sys `static`), so
# the runtime image needs nothing beyond glibc and CA certificates.

FROM rust:1.95-bookworm AS builder
RUN apt-get update && apt-get install -y --no-install-recommends cmake && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY resources ./resources
# The repository's .cargo/config.toml targets the build host's CPU; use a
# portable baseline so the image runs on any x86-64-v2 / aarch64 machine.
ENV RUSTFLAGS="-C target-cpu=x86-64-v2"
RUN cargo build --release --locked && strip target/release/kira-spliceqc

FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="kira-spliceqc" \
      org.opencontainers.image.description="Deterministic, explainable splicing QC for single-cell RNA-seq" \
      org.opencontainers.image.source="https://github.com/ARyaskov/kira-spliceqc" \
      org.opencontainers.image.licenses="MIT"
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates procps && rm -rf /var/lib/apt/lists/*
COPY --from=builder /src/target/release/kira-spliceqc /usr/local/bin/kira-spliceqc
COPY --from=builder /src/resources /opt/kira-spliceqc/resources
WORKDIR /data
ENTRYPOINT ["kira-spliceqc"]
CMD ["--help"]
