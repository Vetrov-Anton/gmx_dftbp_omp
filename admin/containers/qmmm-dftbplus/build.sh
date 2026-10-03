#!/bin/bash
# Build the GROMACS + DFTB+ + PLUMED image on Intel MKL for an Intel or an AMD CPU.
#   ./build.sh amd|intel [image.sif] [more apptainer build options, e.g. --build-arg GMX_SIMD=AVX_512]
# The GROMACS sources are those of the current commit (git archive HEAD); uncommitted changes
# are not included. Needs apptainer >= 1.2 (--build-arg) and --fakeroot.
set -e
vendor=${1:?usage: $0 amd|intel [image.sif] [apptainer build options]}
case $vendor in amd|intel) ;; *) echo "CPU vendor must be amd or intel" >&2; exit 1 ;; esac
out=${2:-gmx-dftbplus-mkl-$vendor.sif}; shift; [ $# -gt 0 ] && shift
here=$(cd "$(dirname "$0")" && pwd)
top=$(git -C "$here" rev-parse --show-toplevel)
[ "${out:0:1}" = / ] || out=$PWD/$out
ctx=$(mktemp -d "${TMPDIR:-/tmp}/gmx-dftbplus-build.XXXXXX")
trap 'rm -rf "$ctx"' EXIT
git -C "$top" archive --format=tar.gz --prefix=gromacs/ -o "$ctx/gromacs-src.tar.gz" HEAD
cp "$here/gmx-dftbplus-mkl.def" "$ctx/"
echo "Building $out for CPU_VENDOR=$vendor from $(git -C "$top" rev-parse --short HEAD)"
cd "$ctx"
apptainer build --fakeroot --build-arg CPU_VENDOR=$vendor --build-arg SRC_TARBALL=gromacs-src.tar.gz "$@" \
    "$out" gmx-dftbplus-mkl.def
