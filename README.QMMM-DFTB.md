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
