#!/bin/bash
# Build the GROMACS + DFTB+ + PLUMED image on Intel MKL for an Intel or an AMD CPU from this checkout.
# (Without a checkout: apptainer build --fakeroot gmx-dftbplus-mkl.sif gmx-dftbplus.def clones the sources.)
#   ./build.sh [auto|amd|intel] [image.sif] [more apptainer build options, e.g. --build-arg GMX_SIMD=AVX_512]
#   ./build.sh auto gmx-dftbplus-openblas.sif --build-arg BLAS=openblas      # OpenBLAS instead of MKL
# auto (default): the CPU vendor is detected at run time, the image works on Intel and AMD.
# The GROMACS sources are those of the current commit (git archive HEAD); uncommitted changes
# are not included. Needs apptainer >= 1.2 (--build-arg) and --fakeroot.
set -e
vendor=${1:-auto}
case $vendor in auto|amd|intel) ;; *) echo "usage: $0 [auto|amd|intel] [image.sif] [apptainer build options]" >&2; exit 1 ;; esac
out=${2:-gmx-dftbplus-mkl.sif}; [ $# -gt 0 ] && shift; [ $# -gt 0 ] && shift
here=$(cd "$(dirname "$0")" && pwd)
top=$(git -C "$here" rev-parse --show-toplevel)
[ "${out:0:1}" = / ] || out=$PWD/$out
ctx=$(mktemp -d "${TMPDIR:-/tmp}/gmx-dftbplus-build.XXXXXX")
trap 'rm -rf "$ctx"' EXIT
git -C "$top" archive --format=tar.gz --prefix=gromacs/ -o "$ctx/gromacs-src.tar.gz" HEAD
# The definition file clones the sources from GitHub; with the tarball in the image it uses that.
cp "$here/gmx-dftbplus.def" "$ctx/"
printf '\n%%files\n    gromacs-src.tar.gz /opt/src/gromacs-src.tar.gz\n' >> "$ctx/gmx-dftbplus.def"
echo "Building $out for CPU_VENDOR=$vendor from $(git -C "$top" rev-parse --short HEAD)"
cd "$ctx"
apptainer build --fakeroot --build-arg CPU_VENDOR=$vendor "$@" \
    "$out" gmx-dftbplus.def
