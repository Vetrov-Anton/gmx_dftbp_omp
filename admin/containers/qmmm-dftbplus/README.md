# Container: GROMACS + DFTB+ + PLUMED on Intel MKL, for Intel or AMD CPUs

`gmx-dftbplus-mkl.def` builds, from `ubuntu:24.04`:

- GROMACS from this tree (the commit `build.sh` is run from), double precision, OpenMPI,
  OpenMP, with the PLUMED patch for GROMACS 2026 and MKL as its BLAS/LAPACK;
- DFTB+ 25.1 (C API, OpenMP, shared) on MKL (`mkl_gnu_thread`: the OpenMP pool of mdrun);
- PLUMED v2.10 with libtorch (runtime-loaded through `PLUMED_KERNEL`).

```bash
./build.sh amd   gmx-dftbplus-mkl-amd.sif
./build.sh intel gmx-dftbplus-mkl-intel.sif
./build.sh amd   gmx-zen4.sif --build-arg GMX_SIMD=AVX_512 --build-arg NJOBS=16
```

| build argument | values | default | effect |
|---|---|---|---|
| `CPU_VENDOR` | `amd`, `intel` | `amd` (set by `build.sh`) | `amd`: `libmklfix.so` (`mkl_serv_intel_cpu_true()` = 1) becomes a `DT_NEEDED` of `gmx` and `dftb+`, so MKL takes its fast code paths on AMD CPUs. `intel`: MKL is used as it is. |
| `GMX_SIMD` | `auto`, any GROMACS `GMX_SIMD` | `auto` = `AVX2_256` (amd), `AVX_512` (intel) | SIMD of GROMACS. Zen 4/5: `AVX_512`; Intel without AVX-512: `AVX2_256`. |
| `NJOBS` | integer | 8 | parallel build jobs |
| `DFTBPLUS_VERSION`, `PLUMED_COMMIT`, `PLUMED_PATCH`, `LIBTORCH_URL` | | 25.1, 53a1773, gromacs-2026.0, libtorch 2.2.0 cpu | versions |

`build.sh` passes the sources as `git archive HEAD`; commit before building.

## Why MKL needs care

- **One LAPACK per process.** GROMACS is built with MKL as its external BLAS/LAPACK. With
  OpenBLAS in libgromacs, the LAPACK calls of DFTB+ are bound to OpenBLAS (libgromacs comes
  first in the symbol lookup order); with the internal LAPACK of GROMACS, they are bound to
  that incomplete copy and DFTB+ crashes in `dlaswp_`. The recipe checks both.
- **AMD.** MKL chooses its kernels by a vendor check; `MKL_DEBUG_CPU_TYPE` has been ignored
  since MKL 2020.1. On a Threadripper 3990X, 600 QM atoms, 8 threads: 3474 ms/step without
  the override, 2495 with it, energies identical. mdrun prints at startup what MKL does:
  `MKL reports Intel CPU code paths (mkl_serv_intel_cpu_true() = 1)` or `generic (non-Intel)`.
  The override uses an internal MKL function, checked with MKL 2020.4 (Ubuntu 24.04); check
  that line again with another MKL.

## Running

```bash
apptainer exec -B /data --env LD_PRELOAD= gmx-dftbplus-mkl-amd.sif \
    gmx mdrun -deffnm md -ntomp 8 -pin on -plumed plumed.dat
```

Use `Solver = DivideAndConquer {}` in the `Hamiltonian = DFTB` block of `dftb_in.hsd`: with
MKL it diagonalises 1.5× faster than the default on 8 threads. See `README.QMMM-DFTB.md` at
the top of the tree for the QM/MM options, threads (`GMX_QMMM_DFTB_NTHREADS`) and timings
(`GMX_QMMM_TIMING`).
