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

## Topology at the QM/MM boundary (grompp)

The QM atoms are found from the `QMMM-grps` of the `.mdp`; their charges are set to
zero and the force-field terms that the QM calculation describes are removed.
Link atoms are QM virtual sites (`[ virtual_sites2 ]`) constructed from a QM atom
(QM1) and an MM atom (MM1).

### `GMX_QMMM_BONDED_SCHEME`

| value | bonded terms removed |
|---|---|
| `classic` (default) | every term in which all but one atom are QM (a QM–QM–MM angle, a QM–QM–QM–MM dihedral), because the QM calculation with the link atom describes it |
| `amber` | only terms whose atoms are all QM; every term with an MM atom is kept at the force-field level |

Bonds between two QM atoms are converted to connections, so the connectivity (and the
exclusions generated from it) stays complete. The scheme is stored in the `tpr`;
mdrun only reminds you of that if it sees the variable.

### `GMX_QMMM_LJ_SCHEME`

| value | Lennard-Jones between QM and MM atoms |
|---|---|
| `forcefield` (default) | by the exclusion rules of the force field (`nrexcl`, `[ pairs ]`); only the LJ within the QM region is excluded |
| `exclude` | in addition, the LJ and LJ-14 of every QM atom with the MM atoms bonded to the QM region are excluded (`classic` only; this was the behaviour of the code before the variable existed) |

No LJ exclusions are generated for the link atoms; give them zero LJ parameters.

### Restraints

Position, flat-bottomed position, distance, orientation, angle and dihedral restraints
and restraint potentials are never removed, also not on QM atoms.

### Output

grompp prints the schemes in use and a summary of the removed terms:

```
QM/MM: force-field terms removed with the 'classic' scheme, by interaction type:
  interaction                all-QM     QM--MM    MM-only
  Bond                            4          0          0
  Angle                           4          3          0
  Proper Dih.                     2          1          0
  LJ-14                           2          0          0
  total                          12          4          0
```

and writes every removed term, kept restraint, converted bond, removed pair and
generated exclusion, atom by atom, to `qmmm_topology_report.txt` (name set with
`GMX_QMMM_TOPOLOGY_REPORT`; `GMX_QMMM_REPORTS=off` switches the report files off).

## QM–MM electrostatics at the boundary (mdrun)

### Atoms at the boundary

| name | meaning |
|---|---|
| QM1 | QM atom covalently bonded to an MM atom |
| MM1 | that MM atom; a link atom is constructed from QM1 and MM1 |
| MM2 | MM atoms bonded to MM1 (not QM) |
| LA | link atom: a virtual site of the QM group constructed from one QM and one MM atom |

### `GMX_QMMM_POT_SCHEME`

With a boundary charge scheme, the charge of every MM1 atom is removed from the QM–MM
electrostatics of every QM atom and replaced as below; `q0 = q(MM1)/n`, with `n` the number
of MM2 atoms of that MM1.

| value | QM–MM electrostatics |
|---|---|
| `none` (default) | every MM charge in full |
| `RC` | MM1 removed; `q0` at the midpoint of every MM1–MM2 bond (Lin & Truhlar 2005) |
| `RCD` | MM1 removed; `2·q0` at every midpoint, and `q(MM2) − q0` on every MM2 (Lin & Truhlar 2005) |
| `CS` | MM1 removed; `q(MM2) + q0` on every MM2, and a pair `+q0/0.12` / `−q0/0.12` on the MM1→MM2 line at 0.94 and 1.06 of the bond length (charge shift, Sherwood *et al.* 2003) |
| `AMBER` | MM1 removed; the MM1 charges of a molecule spread evenly over all other MM atoms of the same molecule |

**The same charges are used for the external potential passed to DFTB+ and for the QM/MM
gradient**, so the forces are the gradient of the energy that DFTB+ returns. The
redistribution exists only in the QM–MM electrostatics: the topology and every MM–MM
interaction keep the charges of the force field.

The fictitious charges have no coordinates of their own: a charge at
`x = (1 − f)·x(MM1) + f·x(MM2)` is a two-atom virtual site, and the force on it is passed to
MM1 and MM2 with the weights `1 − f` and `f`. With PME the fictitious charges act on the QM
atoms with the full `1/r` (they are not on the grid), and the charge of MM1 is removed in the
central cell — its reciprocal-space term is subtracted as a pair term `erf(βr)/r` — while its
periodic images stay. With the cut-off variants the fictitious charges use the same kernel as
the MM atoms. `AMBER` shares are added to the charges of the receiving atoms, on the
short-range list and on the PME grid.

mdrun stops with an error if an MM1 atom has no MM2 atom, if an MM2 atom is itself an MM1
atom or is bonded to a QM atom, if one MM1 atom belongs to two link atoms, or if the molecule
of an MM1 atom has no other MM atom (`AMBER`). It prints the scheme in use,

```
QM/MM electrostatics with the boundary charge scheme CS: 1 MM1 charges removed, 9 fictitious point charges added.
  The external potential passed to DFTB+ and the QM/MM forces are built from the same charges;
  the forces on the point charges are passed to their MM1 and MM2 atoms as for two-atom virtual sites.
  The MM--MM interactions keep the charges of the topology.
```

and writes the link atoms, the removed MM1 charges, the `AMBER` shares and the fictitious
charges to `qmmm_exclusion_report.txt` (`GMX_QMMM_EXCLUSION_REPORT`, `GMX_QMMM_REPORTS=off`).

### `GMX_QMMM_VARIANT`

| value | QM–MM electrostatics |
|---|---|
| `0` | none (vacuum QM) — the default when unset |
| `1` | PME (requires a periodic system) |
| `2` | switched cut-off |
| `3` | reaction field |
| `4` | shifted cut-off |

All boundary charge schemes work with every variant.

### Periodic images of the QM charges (PME)

With `GMX_QMMM_VARIANT=1`, DFTB+ receives the potential of the periodic images of the QM
charges through a callback in every SCC iteration and counts `Σ_A q_A V_img(A)` in its energy
in full, while the Ewald energy of a charge distribution with its own images is
`½ Σ_A q_A V_img(A)` — and the image forces are those of the halved term. mdrun subtracts the
other half from the energy of DFTB+.

### Virial

The QM/MM forces are collected in a buffer of their own, so their virial is supplied by the
QM/MM code: every force is paired with the position it was computed from (the periodic image
nearest to the first QM atom), which is exact for everything computed in real space, the
forces of the fictitious boundary charges included; the reciprocal-space part of PME is
replaced by the exact virial of the grid energy (QM and MM charges together, minus the MM
charges alone). Before, the QM/MM forces were missing from the virial altogether. The extra
work is one more PME call on the steps where the virial is needed. Rectangular boxes only.

### Verification

Central finite differences of the total potential energy against the forces of the same
run (solvated cysteine with one link atom, `SCCTolerance = 1e-10`): with every scheme and
every variant, |F + dE/dx| ≤ 0.025 kJ mol⁻¹ nm⁻¹ for MM1, MM2, QM1 and their neighbours.
The virial, as `dU/dε = 2 tr(Ξ)` under isotropic scaling of all coordinates and of the box
(flexible water): agreement to 0.1 of 7.8·10⁴ kJ/mol with PME and with reaction field;
the previous code missed 270–280 kJ/mol there, about 75 bar on that box.

## DFTB output files (mdrun)

Each of the following is set to a **positive integer stride in steps**; the file is opened
in append mode in the run directory on step 0 and written every *N* steps. Unset means no
file; zero or anything else that is not a positive number is reported and ignored. The
variables that describe the MM environment are ignored when `GMX_QMMM_VARIANT=0`.

| variable | file | content |
|---|---|---|
| `GMX_DFTB_CHARGES=N` | `qm_dftb_charges.xvg` | Mulliken charges of the QM atoms |
| `GMX_DFTB_ATOMIC_SHIFTS=N` | `qm_dftb_atomic_shifts.xvg` | atomic shifts of the QM atoms (only with a DFTB+ that provides `dftbp_get_atomic_shifts()`) |
| `GMX_DFTB_ESP=N` | `qm_dftb_esp.xvg` | potential on each QM atom, MM and QM images summed (V) |
| `GMX_DFTB_ESP_SPLIT=N` | `qm_dftb_esp_split.xvg` | the same potential with the two contributions kept apart |
| `GMX_DFTB_QM_COORD=N` | `qm_dftb_qm.qxyz` | QM coordinates (Å) and charges |
| `GMX_DFTB_MM_COORD=N` | `qm_dftb_mm.qxyz` | charges and coordinates (Å) of the MM atoms **on the short-range list**, as they enter the QM–MM electrostatics, followed by the fictitious charges of the boundary scheme |
| `GMX_DFTB_MM_COORD_FULL=N` | `qm_dftb_mm_full.qxyz` | charges and coordinates of **all** MM atoms (the charges of the PME grid) |
| `GMX_DFTB_QMMM_GRAD=N` | `qm_dftb_grad.xvg` | gradients on the QM atoms and on the short-range MM atoms |
| `GMX_DFTB_QMMM_GRAD_FULL=N` | `qm_dftb_grad_full.xvg` | reciprocal-space gradients on **all** MM atoms (PME only) |
| `GMX_DFTB_ENERGY_CORR=N` | `qm_dftb_energy_corr.xvg` | energy of DFTB+, its correction for the periodic QM images, and the corrected QM energy (kJ/mol) |

`GMX_DFTB_MM_COORD_FULL` and `GMX_DFTB_QMMM_GRAD_FULL` on a solvated system write the whole
box every *N* steps — pick a large stride.

`qm_dftb_esp_split.xvg` has one row per written step: the step, then three numbers per QM
atom, in volts — the potential of the MM atoms (with the boundary scheme), that of the
periodic images of the QM charges (zero unless PME is used), and their sum, which is what
`qm_dftb_esp.xvg` contains. With the cut-off variants the first number can be recomputed
from `qm_dftb_qm.qxyz` and `qm_dftb_mm.qxyz` alone.

`qm_dftb_grad.xvg` has, for every written step, one row per QM atom with the running number
of the QM atom and three vectors — the gradient returned by DFTB+, the electrostatic
gradient due to the environment (with PME including the periodic QM images), and their
sum — and one row per MM atom of the short-range list with its running number, its global
atom number (1-based) and its electrostatic gradient, the forces of the fictitious boundary
charges included. Units are hartree/bohr; these are gradients, the force is their negative
(multiply by 4.96147·10⁴ for kJ mol⁻¹ nm⁻¹). An MM atom of the short-range list with PME
appears in both gradient files, and its gradient is the sum of the two entries.

## Parallel runs: one MPI rank, OpenMP threads

The QM/MM interface runs on **one MPI rank** (no domain decomposition), which is also what
PLUMED needs. All of the parallelism is OpenMP within that rank, in one thread pool:

| part | threads |
|---|---|
| MM forces (nonbonded, bonded, update) | `-ntomp`, as in plain GROMACS |
| PME of the MM system, and the QM/MM PME calls (potential of the MM charges on the QM atoms, potential of the periodic QM images in every SCC iteration, reciprocal-space gradients) | PME threads (`-ntomp`, or `-ntomp_pme`) |
| QM–MM real-space potential and gradients, boundary charge schemes, Ewald exclusion terms, search of the MM atoms near the QM zone, forces and virial | `-ntomp` |
| DFTB+ (SCC, forces) and its BLAS/LAPACK | `-ntomp`, or `GMX_QMMM_DFTB_NTHREADS` |

Run, for example:

```bash
gmx mdrun -deffnm md -ntomp 8 -pin on -pinoffset 100 -pinstride 2 [-plumed plumed.dat]
```

- **Reproducibility.** The loops over the MM atoms are split into contiguous chunks; the
  QM-atom sums of the threads are added in the order of the threads. With one thread the terms
  are added in the order of the serial code. Different thread counts differ only by rounding
  (≈10⁻⁹ kJ/mol on the potential energy of the test system).
- **DFTB+ threads.** DFTB+ runs in the process of mdrun, on the threads of mdrun (by default
  all `-ntomp` of them). `GMX_QMMM_DFTB_NTHREADS=N` sets another number. DFTB+ must not be
  built with MPI (non-MPI DFTB+ uses OpenMP threads by default, `Parallel { UseOmpThreads }`).
- **OpenBLAS.** Use the OpenMP build of OpenBLAS (Debian/Ubuntu: `libopenblas0-openmp`), so
  that LAPACK runs on the same threads. mdrun reports at startup which one it found; with the
  pthreads build it sets the number of OpenBLAS threads to that of DFTB+ and prints a NOTE,
  since that pool is not pinned by mdrun.
- **Pinning.** libgomp lets the surplus threads of its pool exit when a smaller team is
  started (OpenBLAS does that for small matrices) and creates new threads later, which
  inherit the affinity of the main thread. To keep them on the cores of mdrun, the main
  thread of a QM/MM run gets the union of the cores of its pinned threads (printed as
  `QM/MM threads: the main thread may run on the N cores ...`). Without that, the new threads
  of a test run shared one core with the main thread and the step took 10× longer.
- **PLUMED.** Runs on the same single rank. The native PLUMED interface of GROMACS 2026 does
  not pass a thread count to PLUMED (`PLUMED_NUM_THREADS`, default 1); the PLUMED patch for
  GROMACS 2026 does.
- **Timing.** `GMX_QMMM_TIMING=N` prints every N steps the wall time per step of: the
  potential of the MM atoms on the QM atoms, the DFTB+ calls (and the PME potential of the
  QM images within them), the QM/MM gradients. The `QM` row of the cycle accounting in
  `md.log` is the DFTB+ part.

Measured on 12165 atoms (TIP3P box, flexible, PME), cores of one Threadripper 3990X,
ms/step:

| QM zone | 1 thread | 2 | 4 | 8 | of which at 8 threads: DFTB+ / QM–MM terms |
|---|---|---|---|---|---|
| 600 atoms (200 H₂O) | 10064 | 6309 | 4784 | 3987 | 3963 / 58 |
| 150 atoms (50 H₂O) | 395 | 322 | 320 | 308 | 291 / 19 |

The QM–MM terms scale 5–5.6× on 8 threads; the rest is DFTB+, ~87 % of it the
diagonalisation (LAPACK), which scales 2.3–2.6× on 8 threads in a standalone DFTB+ run as
well. A faster threaded LAPACK (e.g. MKL with the GNU OpenMP layer) speeds up exactly that part.
