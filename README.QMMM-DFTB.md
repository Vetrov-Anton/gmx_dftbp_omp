# QM/MM with DFTB+ in this GROMACS tree

GROMACS 2026 with the DFTB+ QM/MM interface of Kubař *et al.* The QM region is
described by DFTB+ called as a library through its C API; the electrostatic
QM/MM coupling (cut-off variants or PME) is computed by GROMACS.

## Building

DFTB+ has to be installed with its C API and as a shared library:

```bash
cmake -B _build -DWITH_API=ON -DBUILD_SHARED_LIBS=ON -DENABLE_DYNAMIC_LOADING=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/dftbplus
cmake --build _build -j && cmake --install _build
```

GROMACS is then configured with

```bash
cmake .. -DGMX_QMMM_PROGRAM=dftbplus -DGMX_QMMM_DFTBPLUS_LIB=/opt/dftbplus ...
```

`GMX_QMMM_DFTBPLUS_LIB` is the install prefix of DFTB+ (its `lib/` subdirectory
is accepted as well). DFTB+ is found through its exported CMake package
(`find_package(DftbPlus CONFIG)`), and linked to `libgromacs`.

### DFTB+ versions

The C API of DFTB+ is not the same in every version, so CMake reads the installed
`dftbplus.h` and adapts the interface to it (`cmake/gmxManageDftbPlus.cmake`):

| DFTB+ | QM atoms and species | `GMX_DFTB_ATOMIC_SHIFTS` |
|---|---|---|
| releases 21.x – 25.x | from the `Geometry` block of `dftb_in.hsd` | not available (ignored with a note) |
| forks whose `dftbp_process_input()` takes an atom list | from the GROMACS topology; `dftb_in.hsd` has no `Geometry` block | available if the fork provides `dftbp_get_atomic_shifts()` |

CMake prints which variant was found, e.g.

```
-- Found DFTB+ (C API 0.4.0) for QM/MM in /opt/dftbplus: QM atoms and species from the
   Geometry block of dftb_in.hsd, atomic shifts not available
```

With a DFTB+ release the `Geometry` block of `dftb_in.hsd` must list the QM atoms
**in the order of the QM group** of the run. mdrun stops if the number of atoms
differs. A `dftb_in.hsd` written for DFTB+ 21 needs one change for newer versions:
`Analysis { CalculateForces = Yes }` becomes `Analysis { PrintForces = Yes }`.
