# NGS-DOSE: counting engine (Rust) + modelling layer (Python) + the GRCh38 resource bundle.
FROM rust:1-bookworm AS build
# clang/libclang: rust-htslib generates its htslib bindings with bindgen at build time
RUN apt-get update && apt-get install -y --no-install-recommends \
      clang libclang-dev pkg-config cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked

FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates libcurl4 libbz2-1.0 liblzma5 zlib1g && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/ngs-dose /usr/local/bin/ngs-dose
WORKDIR /opt/ngs-dose
COPY pyproject.toml README.md ./
COPY ngsdose ./ngsdose
COPY resources ./resources
RUN pip install --no-cache-dir .
ENV NGSDOSE_RESOURCES=/opt/ngs-dose/resources/GRCh38 \
    NGSDOSE_SATELLITES=/opt/ngs-dose/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz
CMD ["ngs-dose", "--help"]
