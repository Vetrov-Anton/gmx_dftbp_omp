# Container: GROMACS + DFTB+ + PLUMED on Intel MKL (Intel or AMD CPUs) or OpenBLAS

`gmx-dftbplus.def` builds, from `ubuntu:24.04`:

- GROMACS from this repository, https://github.com/Vetrov-Anton/gmx_dftbp_omp (branch `main`,
  or the commit `build.sh` is run from), double precision, OpenMPI,
  OpenMP, with the PLUMED patch for GROMACS 2026 and MKL or OpenBLAS as its BLAS/LAPACK;
- DFTB+ 25.1 (C API, OpenMP, shared) on the same library: MKL (`mkl_gnu_thread`) or OpenBLAS
  built with OpenMP, so that its threads are those of the OpenMP pool of mdrun;
- PLUMED v2.10 with libtorch (runtime-loaded through `PLUMED_KERNEL`).

The definition file alone is enough; it clones the GROMACS sources (`GMX_REPO`, `GMX_REF`):

```bash
wget https://raw.githubusercontent.com/Vetrov-Anton/gmx_dftbp_omp/main/admin/containers/qmmm-dftbplus/gmx-dftbplus.def
apptainer build --fakeroot gmx-dftbplus-mkl.sif gmx-dftbplus.def       # MKL; works on Intel and AMD
apptainer build --fakeroot --build-arg BLAS=openblas gmx-dftbplus-openblas.sif gmx-dftbplus.def
apptainer build --fakeroot --build-arg GMX_SIMD=AVX_512 --build-arg NJOBS=16 gmx-zen4.sif gmx-dftbplus.def
```

To build from a local checkout (its current commit, e.g. with your own changes), use `build.sh`:

```bash
git clone https://github.com/Vetrov-Anton/gmx_dftbp_omp.git      # branch main
cd gmx_dftbp_omp/admin/containers/qmmm-dftbplus
./build.sh                                      # -> gmx-dftbplus-mkl.sif
./build.sh auto gmx-zen4.sif --build-arg GMX_SIMD=AVX_512 --build-arg NJOBS=16
./build.sh auto gmx-dftbplus-openblas.sif --build-arg BLAS=openblas
```

| build argument | values | default | effect |
|---|---|---|---|
| `BLAS` | `mkl`, `openblas` | `mkl` | BLAS/LAPACK of DFTB+ and GROMACS: Intel MKL (the fastest) or OpenBLAS with OpenMP threads |
| `CPU_VENDOR` | `auto`, `amd`, `intel` | `auto` | how MKL treats the CPU, see below (`BLAS=mkl` only) |
| `GMX_SIMD` | `auto`, any GROMACS `GMX_SIMD` | `auto` | SIMD of GROMACS; `auto` = detected on the build host. Set it when the image runs on another CPU generation (e.g. `AVX_512` for Zen 4/5 and Xeon, `AVX2_256` for Zen 1–3) |
| `NJOBS` | integer | 8 | parallel build jobs |
| `GMX_REPO`, `GMX_REF` | git URL, branch or tag | this repository, `main` | GROMACS sources to clone (not used by `build.sh`) |
| `DFTBPLUS_VERSION`, `PLUMED_COMMIT`, `PLUMED_PATCH`, `LIBTORCH_URL` | | 25.1, 53a1773, gromacs-2026.0, libtorch 2.2.0 cpu | versions |

`/opt/build_config.txt` in the image records the vendor mode, the SIMD used and the CPU of
the build host.

### CPU vendor (`CPU_VENDOR`)

MKL takes its fast code paths only on CPUs that its internal `mkl_serv_intel_cpu_true()`
identifies as Intel. `libmklfix.so` (`mklfix.c`) defines that function and is a `DT_NEEDED` of
`gmx` and `dftb+`, ahead of MKL, so MKL calls it instead (no `LD_PRELOAD` needed):

- `auto` (default): the vendor is read with CPUID when the program starts. On AMD (and Hygon)
  CPUs the function returns 1; on Intel and others it returns the answer of MKL's own function
  (`dlsym(RTLD_NEXT)`), so MKL behaves exactly as without the library. One image for both.
- `amd`: always 1. `intel`: the library is not built or linked.
- At run time, `MKL_VENDOR_OVERRIDE=off` switches the override off (MKL's own check),
  `MKL_VENDOR_OVERRIDE=on` forces it.
- The QM/MM report files (`qmmm_topology_report.txt` of grompp and
  `qmmm_redistribution_report.txt` of mdrun) are off by default in the image: its environment
  sets `GMX_QMMM_REPORTS=off`. `apptainer exec --env GMX_QMMM_REPORTS=on`,
  `APPTAINERENV_GMX_QMMM_REPORTS=on` on the host, or `env GMX_QMMM_REPORTS=on gmx ...` inside the
  container switches them on; a plain host variable of that name is not passed into the container.

`build.sh` puts the sources into the image as `git archive HEAD`; commit before building.

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
apptainer exec -B /data --env LD_PRELOAD= gmx-dftbplus-mkl.sif \
    gmx mdrun -deffnm md -ntomp 8 -pin on -plumed plumed.dat
```

Use `Solver = DivideAndConquer {}` in the `Hamiltonian = DFTB` block of `dftb_in.hsd`: with
MKL it diagonalises 1.5× faster than the default on 8 threads. See `README.md` at
the top of the tree for the QM/MM options, threads (`GMX_QMMM_DFTB_NTHREADS`) and timings
(`GMX_QMMM_TIMING`).
