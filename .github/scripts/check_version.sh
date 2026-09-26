#!/usr/bin/env bash
# One version everywhere: Cargo.toml (the engine; every counts file records it as engine_version) and
# ngsdose/__init__.py (the Python package; pyproject.toml reads its version from there). Prints it;
# errors go to stderr as workflow annotations.
#
#   .github/scripts/check_version.sh [TAG]          TAG vX.Y.Z: the sources must say X.Y.Z
#   .github/scripts/check_version.sh [TAG] IMAGE    and so must the engine and the package inside IMAGE
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
cargo_v=$(grep -m1 '^version' Cargo.toml | cut -d'"' -f2 || true)
py_v=$(grep -m1 '^__version__' ngsdose/__init__.py | cut -d'"' -f2 || true)
[ -n "$cargo_v" ] && [ "$cargo_v" = "$py_v" ] || { echo "::error::Cargo.toml says $cargo_v, ngsdose/__init__.py says $py_v" >&2; exit 1; }
if [ -n "${1:-}" ] && [ "$1" != "v$cargo_v" ]; then
  echo "::error::tag $1, but the sources say $cargo_v: bump Cargo.toml (and Cargo.lock) and ngsdose/__init__.py first" >&2; exit 1
fi
if [ -n "${2:-}" ]; then
  engine=$(docker run --rm "$2" ngs-dose --version)
  package=$(docker run --rm "$2" ngsdose --version)
  [ "$engine" = "ngs-dose $cargo_v" ] && [ "$package" = "$cargo_v" ] \
    || { echo "::error::image $2 reports '$engine' and ngsdose '$package', the sources $cargo_v" >&2; exit 1; }
fi
echo "$cargo_v"
