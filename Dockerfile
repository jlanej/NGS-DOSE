# NGS-DOSE: counting engine (Rust) + modelling layer (Python) + the GRCh38 resource bundle.
# TODO: pin both base images by digest (FROM rust:1-bookworm@sha256:... and
# python:3.12-slim-bookworm@sha256:...), taking the digests from `docker buildx imagetools inspect <tag>`
# on a machine that can check them, and bump them deliberately. Until then a rebuild of the same
# commit can get a newer rustc or Python; /opt/ngs-dose/pip-freeze.txt records the Python packages
# each image has.
FROM rust:1-bookworm AS build
# clang/libclang: rust-htslib generates its htslib bindings with bindgen at build time
RUN apt-get update && apt-get install -y --no-install-recommends \
      clang libclang-dev pkg-config cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
# the commit this image is built from goes into every counts file it writes (engine_build)
ARG NGSDOSE_BUILD=dev
RUN NGSDOSE_BUILD="$NGSDOSE_BUILD" cargo build --release --locked

# Everything a cluster run needs besides Apptainer, SLURM and the cohort's own scripts (which live
# with the cohort, e.g. NGS-DOSE-1000G): the engine, the ngsdose package, the GRCh38 bundle and
# aria2c for staging CRAMs.
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates libcurl4 libbz2-1.0 liblzma5 zlib1g aria2 curl && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/ngs-dose /usr/local/bin/ngs-dose
WORKDIR /opt/ngs-dose
COPY pyproject.toml README.md ./
COPY ngsdose ./ngsdose
COPY resources ./resources
RUN pip install --no-cache-dir . && pip freeze > /opt/ngs-dose/pip-freeze.txt
COPY docs ./docs
ENV NGSDOSE_RESOURCES=/opt/ngs-dose/resources/GRCh38 \
    NGSDOSE_SATELLITES=/opt/ngs-dose/resources/experimental/satellites.CHM13v2.k31.panel.tsv.gz \
    NGSDOSE_TELOMERE=/opt/ngs-dose/resources/experimental/telomere.k31.panel.tsv.gz
CMD ["ngs-dose", "--help"]
