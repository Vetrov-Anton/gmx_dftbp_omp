/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "qmmm.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#if defined(__linux__)
#    include <sched.h>
#endif
#include <map>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_internal.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdlib/qm_gamess.h"
#include "gromacs/mdlib/qm_gaussian.h"
#include "gromacs/mdlib/qm_mopac.h"
#include "gromacs/mdlib/qm_orca.h"
#include "gromacs/mdlib/qmmm_threading.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
//#include "gromacs/mdtypes/nblist.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunreachable-code"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

void put_cluster_in_MMlist_verlet(int                            ck, // cluster number
                                  int                            na_ck, // # of atoms in cluster
                                  int                            nrQMatoms,
                                  const int*                     indexQM,
                                  const gmx::ArrayRef<const int> atomIndices,
			                      int*                           shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                         pbc,
				                  const rvec*                    x,
                                  const gmx::ArrayRef<const int> globalToLocalAtomMap,
                                  const gmx::ArrayRef<const int> localToGlobalAtomMap);
/*
std::unique_ptr<QMMM_rec>
void init_QMMM_rec(const t_commrec  *cr,
              const gmx_mtop_t *mtop,
              const t_inputrec *ir,
              const t_forcerec *fr,
              const gmx_wallcycle_t gmx_unused wcycle)
{
    return std::make_unique<QMMM_rec>(cr, mtop, ir, fr, wcycle);
}
*/

static real call_QMroutine(//const t_commrec*  cr,
                           QMMM_rec*         qr,
                           QMMM_QMrec*       qm,
                           QMMM_MMrec*       mm,
                           rvec              f[],
                           rvec              fshift[],
                           t_nrnb*           nrnb,
                           gmx_wallcycle*    wcycle)
{
    // Makes a call to the requested QM routine (qm->QMmethod).
    // Note that f is actually the gradient, i.e. -f

    if (GMX_QMMM_MOPAC)
    {
        return call_mopac(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_GAMESS)
    {
        return call_gamess(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_GAUSSIAN)
    {
        return qm->gaussian.call_gaussian(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_ORCA)
    {
        return call_orca(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_DFTBPLUS)
    {
        return call_dftbplus(qr, qm, *mm, f, fshift, nrnb, wcycle); // cr
    }
    else
    {
        gmx_fatal(FARGS, "Unknown QM software -- should never land here :-/");
    }
}

// Update QM and MM coordinates in the QM/MM data structures.
// New version of the function:
// Update the coordinates of the MM atoms on the short-range neighborlist!
// The NBlist needs to have been created previously by either group or Verlet scheme.
void QMMM_rec::update_QMMM_coord(const t_commrec*  cr,
                                 rvec*             shift_vec,
                                 const rvec        x[],
                                 // gmx::ArrayRef<const gmx::RVec> x,
                                 const t_mdatoms*  md,
                                 const matrix      box)
{
    // Shifts the QM and MM atoms into the central box and
    //   stores the shifted coordinates in the coordinate arrays of QMMMrec.
    // These coordinates are passed on the QM subroutines.
    //
    // Only MM atoms up to the distance fr->rcoulomb from the respective
    //   nearest QM atoms are considered;
    // in case fr->rcoulomb == 0. is detected,
    //   all of the MM atoms are considered.

    QMMM_QMrec& qm_ = qm[0];
    QMMM_MMrec& mm_ = mm[0];
    real rcut = qm_.rcoulomb > 0.1 ? qm_.rcoulomb : 999999.; // infinity
    // char, not bool: the elements are written by several threads
    std::vector<char> isCurrentMMatom(mm_.nrMMatoms_nbl, 0);

 // printf("Original Gromacs coordinates\n");
 // for (int i = 0; i < qm_.nrQMatoms; i++)
 // {
 //     printf("QM atom %d: %8.5f %8.5f %8.5f\n", qm_.indexQM[i]+1,
 //                                                 x[globalToLocalAtomMap[qm_.indexQM[i]]][XX],
 //                                                 x[globalToLocalAtomMap[qm_.indexQM[i]]][YY],
 //                                                 x[globalToLocalAtomMap[qm_.indexQM[i]]][ZZ]);
 // }

    // shift the QM atoms into the central box
    for (int i = 0; i < qm_.nrQMatoms; i++)
    {
        rvec_sub(x[globalToLocalAtomMap[qm_.indexQM[i]]], shift_vec[qm_.shiftQM[i]], qm_.xQM[i]);
    }

 // printf("QM coordinates updated\n");
 // for (int i = 0; i < qm_.nrQMatoms; i++)
 // {
 //     printf("QM atom %d: %8.5f %8.5f %8.5f\n", qm_.indexQM[i]+1, qm_.xQM[i][XX], qm_.xQM[i][YY], qm_.xQM[i][ZZ]);
 // }

    // copy box size
    copy_mat(box, qm_.box);

    // initialize PBC for MM coordinate manipulation
    t_pbc pbc;
    gmx::IVec null_ivec;
    clear_ivec(null_ivec);
    set_pbc_dd(&pbc, pbcType, haveDDAtomOrdering(*cr) ? &cr->dd->numCells : &null_ivec, false, box);

 // for (int s = 0; s < pbc.ntric_vec; s++)
 // {
 //     printf("SHIFT[%2d] = %d %d %d\n", s, pbc.tric_shift[s][0], pbc.tric_shift[s][1], pbc.tric_shift[s][2]);
 //  // printf("SHIFT[%2d] = %8.5f %8.5f %8.5f\n", s, pbc.tric_vec[s][0], pbc.tric_vec[s][1], pbc.tric_vec[s][2]);
 // }

    // DECIDE IF WE WANT TO APPLY A CUTOFF ON THE ATOMS FROM THE SR NEIGHBORLIST!

    // DO WE NEED TO RE-ALLOCATE THE ARRAYS TO BE FILLED?
    //   YES!
    //
    // FIRST, IDENTIFY THE MM ATOMS UP TO CUTOFF AT THIS STEP AND COUNT THEM:

    // Among the atoms found as candidates for being MM atoms in neighborsearching,
    // find those that are within electrostatics cut-off.
    // For the cutoff, use the value "rcut"
    int nrMMatoms = 0;
#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_.nrMMatoms_nbl, 256)) schedule(static) reduction(+ : nrMMatoms)
    for (int i = 0; i < mm_.nrMMatoms_nbl; i++)
    {
	    isCurrentMMatom[i] = false;
	 // printf("DEBUG_MM TEST %5d %5d %2d", i, mm_.indexMM_nbl[i], mm_.shiftMM_nbl[i]);
	    // loop over all QM atoms here
	    for (int q=0; q<qm_.nrQMatoms; q++)
	    {
            rvec bond;
            pbc_dx_aiuc(&pbc, x[globalToLocalAtomMap[qm_.indexQM[q]]], x[globalToLocalAtomMap[mm_.indexMM_nbl[i]]], bond);
	     // printf(" %8.5f\n", norm(bond));
	        if (norm(bond) < rcut)
            {
	            isCurrentMMatom[i] = true;
	            nrMMatoms++;
	         // printf("DEBUG_MM %5d %5d %2d %6.4f\n", i, mm_.indexMM_nbl[i], mm_.shiftMM_nbl[i], distance);
                break;
            }
	    }
	 // printf("\n");
    }
 // printf("Number of actual    MM atoms in the current MD step               : %d\n", nrMMatoms);

    // ALLOCATION
    mm_.nrMMatoms = nrMMatoms;
    mm_.indexMM.resize(nrMMatoms);
    mm_.MMcharges.resize(nrMMatoms);
    mm_.shiftMM.resize(nrMMatoms);
    mm_.xMM.resizeWithPadding(nrMMatoms);

    int index = 0; // runs over the identified MM atoms
    for (int i = 0; i < mm_.nrMMatoms_nbl; i++)
    {
	    if (isCurrentMMatom[i])
	    {
	        // Add to list!
	        mm_.indexMM[index] = mm_.indexMM_nbl[i];

	        // Also add charge
	        mm_.MMcharges[index] = md->chargeA[globalToLocalAtomMap[mm_.indexMM[index]]] * mm_.scalefactor;

            // Having obtained the shift at NS time (update_qmmmrec),
            //   merely copy it here to shiftMM[]
	        mm_.shiftMM[index] = mm_.shiftMM_nbl[i];

	        // one MM atom found => increment index */
	        index++;
	    }
    }

    // The short-range MM list has just been rebuilt: map the charges of the QM--MM
    //   electrostatics (boundary charge scheme) onto it.
    update_QMMM_boundary_SR();

    // also shift the MM atoms into the central box

 //   for (int a=0; a<45; a++)
 //     printf("SHIFT %2d: %7.3f %7.3f %7.3f\n", a,
 //     fr->shift_vec[a][XX], fr->shift_vec[a][YY], fr->shift_vec[a][ZZ]);

#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_.nrMMatoms, 1024)) schedule(static)
    for (int ind = 0; ind < mm_.nrMMatoms; ind++)
    {
        rvec_sub(x[globalToLocalAtomMap[mm_.indexMM[ind]]], shift_vec[mm_.shiftMM[ind]], mm_.xMM[ind]);
 //     printf("COORD MM %4d %2d\n", mm_.indexMM[ind], mm_.shiftMM[ind]);
    }

    // For DFTB, also update the coordinates of *all* of the MM atoms,
    //   not only those on the short-range neighborlist.
    // Do not shift the MM atoms into the central box!
    //   It might break the calculation of the surface correction in the Ewald sum.
    if (GMX_QMMM_DFTBPLUS)
    {
#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_.nrMMatoms_full, 1024)) schedule(static)
        for (int i = 0; i < mm_.nrMMatoms_full; i++)
        {
            copy_rvec(x[globalToLocalAtomMap[mm_.indexMM_full[i]]], mm_.xMM_full[i]);
        }
    }
} // update_QMMM_coord

void QMMM_QMrec::init_QMrec(int               grpnr,
                            int               nr,
                            const int*        atomarray,
                            const gmx_mtop_t* mtop,
                            const t_inputrec* ir)
{
    nrQMatoms = nr;
    snew(xQM, nr);
    snew(indexQM, nr);
    snew(shiftQM, nr);
    for (int i = 0; i < nr; i++)
    {
        indexQM[i] = atomarray[i];
    }

    snew(atomicnumberQM, nr);
    int molb = 0;
    for (int i = 0; i < nrQMatoms; i++)
    {
        const t_atom &atom = mtopGetAtomParameters(*mtop, indexQM[i], &molb);
        nelectrons        += atom.atomnumber; // mtop->atomtypes.atomnumber[atom.type];
        atomicnumberQM[i]  = atom.atomnumber; // mtop->atomtypes.atomnumber[atom.type];
    }

    QMcharge      = ir->opts.QMcharge[grpnr];
    multiplicity  = ir->opts.QMmult[grpnr];
    nelectrons   -= ir->opts.QMcharge[grpnr];

    QMmethod      = ir->opts.QMmethod[grpnr];
    QMbasis       = ir->opts.QMbasis[grpnr];

    // hack to prevent gaussian from reinitializing all the time
    gaussian.nQMcpus = 0; // number of CPU's to be used by g01, is set
                          // upon initializing gaussian with init_gaussian()

    rcoulomb      = ir->rcoulomb;
    ewaldcoeff_q  = calc_ewaldcoeff_q(ir->rcoulomb, ir->ewald_rtol);
    epsilon_r     = ir->epsilon_r;

    snew(pot_qmmm, nr);
    snew(pot_qmqm, nr);

} // init_QMrec

int QMMM_QMrec::nrQMatoms_get() const
{
    return nrQMatoms;
}

int QMMM_QMrec::qmmm_variant_get() const
{
    return qmmm_variant;
}

double QMMM_QMrec::xQM_get(const int atom, const int coordinate) const
{
    return xQM[atom][coordinate];
}

real QMMM_QMrec::QMcharges_get(const int atom) const
{
    return QMcharges[atom];
}

void QMMM_QMrec::QMcharges_set(const int atom, const real value)
{
    QMcharges[atom] = value;
}

real QMMM_QMrec::QMatomicShifts_get(const int atom) const
{
    return QMatomicShifts[atom];
}

void QMMM_QMrec::QMatomicShifts_set(const int atom, const real value)
{
    QMatomicShifts[atom] = value;
}

double QMMM_QMrec::pot_qmmm_get(const int atom) const
{
    return pot_qmmm[atom];
}

double QMMM_QMrec::pot_qmqm_get(const int atom) const
{
    return pot_qmqm[atom];
}

void QMMM_QMrec::pot_qmmm_set(const int atom, const double value)
{
    pot_qmmm[atom] = value;
}

void QMMM_QMrec::pot_qmqm_set(const int atom, const double value)
{
    pot_qmqm[atom] = value;
}

int QMMM_QMrec::atomicnumberQM_get(const int atom)const
{
    return atomicnumberQM[atom];
}

int QMMM_QMrec::QMcharge_get()const
{
    return QMcharge;
}

int QMMM_QMrec::multiplicity_get()const
{
    return multiplicity;
}

QMmethodType QMMM_QMrec::QMmethod_get()const
{
    return QMmethod;
}

QMbasisType QMMM_QMrec::QMbasis_get()const
{
    return QMbasis;
}

int QMMM_QMrec::nelectrons_get()const
{
    return nelectrons;
}

int QMMM_QMrec::CASelectrons_get()const
{
    return CASelectrons;
}

int QMMM_QMrec::CASorbitals_get()const
{
    return CASorbitals;
}

real QMMM_QMrec::box_xx_get() const
{
    return box[0][0];
}

real QMMM_QMrec::box_yy_get() const
{
    return box[1][1];
}

real QMMM_QMrec::box_zz_get() const
{
    return box[2][2];
}

void QMMM_MMrec::init_MMrec(real scalefactor_in,
                            int  nrMMatoms_full_in,
                            int  natoms,
                            int  nrQMatoms,
                            const int* indexQM,
                            int* found_mm_atoms)
{
    scalefactor    = scalefactor_in;
    nrMMatoms_full = nrMMatoms_full_in;
    indexMM_full.resize(nrMMatoms_full); // ???
    xMM_full.resizeWithPadding(nrMMatoms_full); // ???
    MMcharges_full.resize(nrMMatoms_full); // ???
    shiftMM_full.resize(nrMMatoms_full); // ???

    // fill the indexMM_full array
    *found_mm_atoms = 0;
    for (int i=0; i<natoms; i++)
    {
        bool is_mm_atom = true;
        for (int j=0; j<nrQMatoms; j++)
        {
            if (i == indexQM[j])
            {
                 is_mm_atom = false;
            }
        }
        if (is_mm_atom)
        {
            indexMM_full[*found_mm_atoms] = i;
	        (*found_mm_atoms)++;
        }
    }
}

QMMM_rec::QMMM_rec(const t_commrec*                 cr,
                   const gmx_mtop_t*                mtop,
                   const t_inputrec*                ir,
                   const t_forcerec*                fr)
 //                const gmx_wallcycle*  gmx_unused wcycle)
    : nAtoms(mtop->natoms)
{
#if GMX_QMMM
    // Put the atom numbers of atoms that belong to the QMMM group
    // into an array that will be copied later to QMMMrec->indexQM[..].
    // Also, it will be used to create an index array QMMMrec->bQMMM[],
    // which contains true/false for QM and MM (the other) atoms.

    if (!GMX_QMMM)
    {
        gmx_incons("Compiled without QMMM");
    }

    // issue a fatal if the user wants to run with more than one node
    if (cr->commMyGroup.isParallel())
    {
        gmx_fatal(FARGS, "QM/MM may not work in parallel due to neighborsearching issues, \
              use a single processor instead!\n");
    }

    // The array bQMMM[] contains true/false for atoms that are QM/not QM.
    // We first set all elements at false.
    // Afterwards we use qm_arr (= MMrec->indexQM) to change
    // the elements corresponding to the QM atoms at true.

    // We take the possibility into account
    // that a user has defined more than one QM group:
    // HOW SHOULD WE PROCEED IN THAT CASE?
    // IT WOULD BE COOL TO BE ABLE TO DO IT!

    // An ugly work-around in case there is only one group.
    // In this case, the whole system is treated as QM.
    // Otherwise, the second group is always the rest of the total system
    //   and is treated as MM.

    // Small problem if there is only QM... so no MM. */

    pbcType = fr->pbcType;

    int numQmmmGroups = ir->opts.ngQM;

    if (numQmmmGroups > 1) {
        fprintf(stderr, "\nQM/MM cannot calculate more than 1 group of atoms at the moment\nExiting!\n\n");
        exit(-1);
    }

    // There are numQmmmGroups groups of QM atoms.
    // Previously, multiple QM groups typically meant
    // that the user wanted to do ONIOM.
    // However, maybe it should also be possible to define
    // more than one QM subsystem with independent neighbourlists.
    // Gerrit Groenhof said he would have to think about that...
    // (11-11-2003)

    std::vector<int> qmmmAtoms = qmmmAtomIndices(*ir, *mtop);

    qm.resize(numQmmmGroups);

    // Standard QMMM (no ONIOM).
    // All layers are merged together, so there is one QM subsystem and one MM subsystem.
    // Also, we set the charges to zero in mtop
    //   to prevent the innerloops from doubly counting the electrostatic QM--MM interaction.
    // TODO: Consider doing this in grompp instead.

    // store QM atoms in the QMrec and initialise
    qm[0].init_QMrec(0, qmmmAtoms.size(), qmmmAtoms.data(), mtop, ir);

    // print the current layer to allow users to check their input
    fprintf(stderr, "Layer %d\nnr of QM atoms %d\n", 0, qm[0].nrQMatoms);
    fprintf(stderr, "QMlevel: %s/%s\n\n",
            enumValueToString(qm[0].QMmethod),
            enumValueToString(qm[0].QMbasis));

    // MM rec creation
    int nrMMatoms_full_in = (mtop->natoms)-(qm[0].nrQMatoms); // rest of the atoms
    int found_mm_atoms = 0;
    mm.resize(1);
    QMMM_MMrec& mm_ = mm[0];
    mm_.init_MMrec(ir->scalefactor, nrMMatoms_full_in, mtop->natoms, qm[0].nrQMatoms, qm[0].indexQM, &found_mm_atoms); 

    printf ("(mtop->natoms) = %d\n(qr->qm[0]->nrQMatoms) = %d\nmm->nrMMatoms_full = %d\n",
            (mtop->natoms), (qm[0].nrQMatoms), mm_.nrMMatoms_full);
    printf ("(found_mm_atoms) = %d\n", found_mm_atoms);

    // these variables get updated in the update QMMMrec // ???

    // OLD COMMENT but maybe useful in the future:
    //   With only one layer there is only one initialization needed.
    //   Multilayer is a bit more complicated as it requires
    //   a re-initialization at every step of the simulation.
    //   This is due to the use of COMMON blocks in Fortran QM subroutines.

    if (GMX_QMMM_MOPAC)
    {
        init_mopac(qm[0]);
    }
    else if (GMX_QMMM_GAMESS)
    {
        init_gamess(cr, qm[0], mm_);
    }
    else if (GMX_QMMM_GAUSSIAN)
    {
        qm[0].gaussian.init_gaussian();
    }
    else if (GMX_QMMM_ORCA)
    {
        init_orca(&(qm[0]));
    }
    else if (GMX_QMMM_DFTBPLUS)
    {
        // The bonded and LJ interactions at the QM/MM boundary are treated according to
        //   the schemes selected in grompp, and the outcome is stored in the tpr file.
        if (getenv("GMX_QMMM_BONDED_SCHEME") != nullptr || getenv("GMX_QMMM_LJ_SCHEME") != nullptr)
        {
            fprintf(stdout,
                    "NOTE: GMX_QMMM_BONDED_SCHEME and GMX_QMMM_LJ_SCHEME are evaluated by grompp, "
                    "not by mdrun.\n"
                    "      The treatment of the bonded and LJ interactions at the QM/MM boundary "
                    "is fixed in the tpr file.\n");
        }

        // Look how the QM/MM electrostatics shall be treated.
        // In the future, this could be performed for QM/MM in general,
        //   not only with DFTB+.
        char *env1 = getenv("GMX_QMMM_VARIANT");
        char *env2 = getenv("GMX_QMMM_PME_DIPCOR");
        if (env1 == nullptr)
        {
            qm[0].qmmm_variant = eqmmmVACUO;
		    fprintf(stdout, "No electrostatic QM/MM interaction.\nTo change, set environment variable GMX_QMMM_VARIANT.\n");
        }
        else
        {
            sscanf(env1, "%d", &(qm[0].qmmm_variant));
            switch (qm[0].qmmm_variant) {
		    case eqmmmVACUO: // 0
		                    fprintf(stdout, "No electrostatic QM/MM interaction.\n");
		                    break;
		    case eqmmmPME: // 1
                {
		               if (pbcType != PbcType::Xyz)
                       {
		                   fprintf(stderr, "PME treatment of QM/MM electrostatics only possible with triclinic periodic system!\n");
		                   exit(-1);
		               }
		               fprintf(stdout, "Electrostatic QM/MM interaction calculated with full PME treatment.\n");

                       pme.resize(2);
                       pmedata = nullptr; // This will be initialized later in the PME routines
                       QMMM_PME& pme_full   = pme[0];
                       QMMM_PME& pme_qmonly = pme[1];

                       // PME data structure for the entire system
                       pme_full.x.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.q.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.f.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       snew(pme_full.pot, qm[0].nrQMatoms);
                       
                       // PME data structure for the QM-only system
                       pme_qmonly.x.resizeWithPadding(qm[0].nrQMatoms);
                       pme_qmonly.q.resizeWithPadding(qm[0].nrQMatoms);
                       pme_qmonly.f.resizeWithPadding(qm[0].nrQMatoms);
                       snew(pme_qmonly.pot, qm[0].nrQMatoms);
                       
                       if (env2 != nullptr)
                       {
                           pme_full.surf_corr_pme   = true;
                           pme_full.epsilon_r       = qm[0].epsilon_r;
                           pme_qmonly.surf_corr_pme = true;
                           pme_qmonly.epsilon_r     = qm[0].epsilon_r;
					       fprintf(stdout, "Dipole (surface) correction for QM/MM PME applied ");
					       fprintf(stdout, "with a permittivity of %5.1f.\n", pme_qmonly.epsilon_r);
					       fprintf(stdout, "\nCurrently disabled due to solvent molecules broken across box boundary!\nExiting!\n\n");
                           exit(-1);
                       }
                       else
                       {
                           pme_full.surf_corr_pme   = false;
                           pme_qmonly.surf_corr_pme = false;
					       fprintf(stdout, "No dipole (surface) correction for QM/MM PME, i.e. tin-foil boundary conditions.\n");
                       }
		               break;
                }
			case eqmmmSWITCH: // 2
		                 fprintf(stdout, "Electrostatic QM/MM interaction calculated with a switched cut-off.\n");
		                 break;
		    case eqmmmRFIELD: // 3
		                 fprintf(stdout, "Electrostatic QM/MM interaction calculated with a reaction-field cut-off.\n");
		                 break;
		    case eqmmmSHIFT: // 4
		                 fprintf(stdout, "Electrostatic QM/MM interaction calculated with a shifted cut-off.\n");
		                 break;
		    default:
		            fprintf(stderr, "Unrecognized choice for treatment of QM/MM electrostatics.\n");
		            fprintf(stderr, "Set environment variable GMX_QMMM_VARIANT to either 0, 1, 2, 3, or 4.\n");
	                exit(-1);
		    }
        }
        snew(qm[0].QMcharges, qm[0].nrQMatoms);
        snew(qm[0].QMatomicShifts, qm[0].nrQMatoms);

        // Boundary charge scheme of the QM--MM electrostatics
        init_QMMM_boundary(mtop);

        init_dftbplus(&(qm[0]), this, ir, cr); //, wcycle);
    }
    else
    {
        gmx_fatal(FARGS, "Unknown QM software -- should never land here :-/");
    }
#else // GMX_QMMM
    gmx_incons("Compiled without QMMM");
    (void) cr;
    (void) mtop;
    (void) ir;
    (void) fr;
#endif
} // init_QMMMrec

QMMM_rec::~QMMM_rec() = default;

int qmmmShareThreadAffinity(int numThreads)
{
#if defined(__linux__)
    if (numThreads <= 1)
    {
        return 0;
    }
    std::vector<cpu_set_t> masks(numThreads);
    std::vector<int>       ok(numThreads, 0);
#    pragma omp parallel num_threads(numThreads)
    {
        const int t = gmx_omp_get_thread_num();
        CPU_ZERO(&masks[t]);
        ok[t] = (sched_getaffinity(0, sizeof(cpu_set_t), &masks[t]) == 0);
    }
    cpu_set_t all;
    CPU_ZERO(&all);
    for (int t = 0; t < numThreads; t++)
    {
        if (!ok[t])
        {
            return 0;
        }
        CPU_OR(&all, &all, &masks[t]);
    }
    const int numCores = CPU_COUNT(&all);
    if (numCores <= CPU_COUNT(&masks[0]) || sched_setaffinity(0, sizeof(cpu_set_t), &all) != 0)
    {
        return 0;
    }
    std::string list;
    for (int c = 0; c < CPU_SETSIZE && gmx::ssize(list) < 200; c++)
    {
        if (CPU_ISSET(c, &all))
        {
            list += (list.empty() ? "" : ",") + std::to_string(c);
        }
    }
    printf("QM/MM threads: the main thread may run on the %d cores of the %d pinned threads (%s),\n"
           "  so that the OpenMP threads that DFTB+ and OpenBLAS create later stay on them.\n",
           numCores, numThreads, list.c_str());
    fflush(stdout);
    return numCores;
#else
    (void) numThreads;
    return 0;
#endif
}

std::vector<int> qmmmAtomIndices(const t_inputrec& ir, const gmx_mtop_t& mtop)
{
    const int               numQmmmGroups = ir.opts.ngQM;
    const SimulationGroups& groups        = mtop.groups;
    std::vector<int>        qmmmAtoms;
    for (int i = 0; i < numQmmmGroups; i++)
    {
        for (const AtomProxy atomP : AtomRange(mtop))
        {
            int index = atomP.globalAtomNumber();
            if (getGroupType(groups, SimulationAtomGroupType::QuantumMechanics, index) == i)
            {
                qmmmAtoms.push_back(index);
            }
        }
    }
    return qmmmAtoms;
}

void removeQmmmAtomCharges(gmx_mtop_t* mtop, gmx::ArrayRef<const int> qmmmAtoms)
{
    int molb = 0;
    for (gmx::Index i = 0; i < qmmmAtoms.ssize(); i++)
    {
        int indexInMolecule;
        mtopGetMolblockIndex(*mtop, qmmmAtoms[i], &molb, nullptr, &indexInMolecule);
        t_atom* atom = &mtop->moltype[mtop->molblock[molb].type].atoms.atom[indexInMolecule];
        atom->q      = 0.0;
        atom->qB     = 0.0;
    }
}

// Updates the shift and charges of *all of the* MM atoms in QMMMrec.
//   Only with DFTB.
//   (Not nice, should be done in a more elegant way...)
void QMMM_rec::update_QMMMrec_dftb(const t_commrec*  cr,
                                   rvec*             shift_vec,
                                   const rvec        x[],
                                   const t_mdatoms*  md,
                                   const matrix      box)
{
    // INHERITED NOTE: is NOT yet working if there are no PBC.
    // Also in ns.c, simple NS needs to be fixed!
    //   As of 2019, ns.c does not exist any longer.

    // copy pointers
    QMMM_QMrec& qm_ = qm[0]; // in case of normal QMMM, there is only one group
    QMMM_MMrec& mm_ = mm[0];

    // init_pbc(box); needs to be called first, see pbc.h
    gmx::IVec null_ivec;
    clear_ivec(null_ivec);
    t_pbc pbc;
//  set_pbc_dd(&pbc, pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : null_ivec, false, box);
    set_pbc_dd(&pbc, pbcType, haveDDAtomOrdering(*cr) ? &cr->dd->numCells : &null_ivec, false, box);

 // printf("There are %d QM atoms, namely:", qm_.nrQMatoms);
 // for (int i=0; i<qm_.nrQMatoms; i++)
 //   printf(" %d", qm_.indexQM[i]);
 // printf("\n");

    // Compute the shift for the MM atoms with respect to QM atom [0].
    // TODO: This looks like a viable first guess, but is that correct?
    // Related to the problem of contributions to virial pressure
    //   in a system treated with particle--mesh Ewald.
    rvec crd;
    rvec_sub(x[globalToLocalAtomMap[qm_.indexQM[0]]], shift_vec[qm_.shiftQM[0]], crd);
#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_.nrMMatoms_full, 1024)) schedule(static)
    for (int i=0; i<mm_.nrMMatoms_full; i++) {
        rvec dx;
        mm_.shiftMM_full[i] = pbc_dx_aiuc(&pbc, crd, x[globalToLocalAtomMap[mm_.indexMM_full[i]]], dx);
    }

 // // previous version of the loop
 // for (i=0; i<mm_.nrMMatoms; i++) {
 //     ivec dx;
 //     current_shift = pbc_dx_aiuc(&pbc, x[globalToLocalAtomMap[qm_.indexQM[0]]], x[globalToLocalAtomMap[mm_.indexMM[i]]], dx);
 //     crd[0] = IS2X(QMMMlist->shift[i]) + IS2X(qm_i_particles[i].shift);
 //     crd[1] = IS2Y(QMMMlist->shift[i]) + IS2Y(qm_i_particles[i].shift);
 //     crd[2] = IS2Z(QMMMlist->shift[i]) + IS2Z(qm_i_particles[i].shift);
 //     is     = static_cast<int>(XYZ2IS(crd[0], crd[1], crd[2]));
 //     mm_.shiftMM[i] = is;
 // }

#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_.nrMMatoms_full, 1024)) schedule(static)
    for (int i = 0; i < mm_.nrMMatoms_full; i++) // no free energy yet
    {
        mm_.MMcharges_full[i] = md->chargeA[globalToLocalAtomMap[mm_.indexMM_full[i]]] * mm_.scalefactor;
    }
    update_QMMM_boundary_full();
} // update_QMMMrec_dftb

// ADD THE NON-QM ATOMS IN THE VERLET CLUSTER ck TO THE LIST OF MM ATOMS
void put_cluster_in_MMlist_verlet(int                            ck, // cluster number
                                  int                            na_ck, // # of atoms in cluster
                                  int                            nrQMatoms,
                                  const int*                     indexQM,
                                  const gmx::ArrayRef<const int> atomIndices,
			                      int*                           shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                         pbc,
				                  const rvec*                    x,
                                  const gmx::ArrayRef<const int> globalToLocalAtomMap,
                                  const gmx::ArrayRef<const int> localToGlobalAtomMap)
{
 //  * This calculation of shift would be desirable,
 //  * but it does not seem to work properly!
 // ivec crd;
 // crd[XX] = (is_j_cluster ? 1 : -1) * IS2X(shift) + IS2X(qm->shiftQM[qm_atom]);
 // crd[YY] = (is_j_cluster ? 1 : -1) * IS2Y(shift) + IS2Y(qm->shiftQM[qm_atom]);
 // crd[ZZ] = (is_j_cluster ? 1 : -1) * IS2Z(shift) + IS2Z(qm->shiftQM[qm_atom]);
 // int is = IVEC2IS(crd);

    // Loop over the atoms in the cluster ck.
    for (int k=0; k<na_ck; k++)  // NA_CK IS USUALLY 4 (SIMD RELATED)
    {
	    const int localAtom = atomIndices[na_ck * ck + k];
	    if (localAtom < 0)
	    {
	        // The value of -1 in the Verlet list is for padding purpose only.
	        // It does not correspond to any atom.
	        // Therefore, ignore!
	        continue;
	    }
	    const int globalAtom = localToGlobalAtomMap[localAtom];
	    if (globalAtom < 0)
	    {
	        continue;
	    }
	    // In the following loop, determine 2 things:
	    // 1: the shift to put the k-th atom (a putative MM atom)
	    //    shortest-distance with respect to the nearest QM atom
	    // 2: whether the k-th atom is a QM atom
	    real dist = 1000.;
	    int sh = -1;
	    bool is_qmatom = false;
	    for (int q=0; q<nrQMatoms; q++)
	    {
	        // 1: the shift -- this calculation looks OK!
	        rvec bond;
	        int sh_t = pbc_dx_aiuc(pbc,
                                   x[globalToLocalAtomMap[indexQM[q]]],
                                   x[globalToLocalAtomMap[globalAtom]],
                                   bond);
	        if (norm(bond) < dist)
	        {
	            dist = norm(bond);
		        sh = sh_t;
	        }
	        // 2: a QM atom?
	        if (globalAtom == indexQM[q])
	        {
	            is_qmatom = true;
            }
	    }
	    // If it is not a QM atom, then put it in the list and store the shift.
	    if (!is_qmatom)
	    {
	        // printf("FOUND_MM_ATOM %5d in cluster %4d\n", globalAtom, ck);
	        shiftMMatom[globalAtom] = sh; // true;
	    }
    }
}

// Construct the maps for finding the local atom index of a global atom index and vice versa.
void QMMM_rec::update_QMMMrec_map(const t_commrec* cr) //, int nAtoms)
{
    globalToLocalAtomMap.resize(nAtoms, -1);
    localToGlobalAtomMap.resize(nAtoms, -1);

    if (haveDDAtomOrdering(*cr))
    {
        for (int globalAtom = 0; globalAtom < nAtoms; ++globalAtom)
        {
            if (const int* home = cr->dd->ga2la->findHome(globalAtom))
            {
                globalToLocalAtomMap[globalAtom] = *home;
                localToGlobalAtomMap[*home] = globalAtom;
            }
            else
            {
                gmx_fatal(FARGS, "Global atom %d has no home in domain decomposition",
                    globalAtom);
            }
        }
    }
    else
    {
        // no DD, identity mapping
        for (int iAtom = 0; iAtom < nAtoms; ++iAtom)
        {
            globalToLocalAtomMap[iAtom] = iAtom;
            localToGlobalAtomMap[iAtom] = iAtom;
        }
    }
}

// create the SR MM list using the Verlet neighborlist
void QMMM_rec::update_QMMMrec_verlet_ns(const t_commrec*    cr,
                                        const gmx::nonbonded_verlet_t* nbv,
                                        const rvec          x[],
                                        const t_mdatoms*    md,
                                        const matrix        box)
{
    (void) md;
 //  * COMMENTS TO THE FORMER GROUP-SCHEME BASED VERSION OF THIS FUNCTION:
 //  *********************************************************************
 //  * Create/update a number of QMMMrec entries:
 //  * 1) shiftQM -- shifts of the QM atoms
 //  * 2) indexMM -- indices of the MM atoms
 //  * 3) shiftMM -- shifts of the MM atoms
 //  * 4) shifted coordinates of the MM atoms
 //  *       --- NOT THIS ONE, BECAUSE A SEPARATE ROUTINE IS USED!
 //  * (The shifts are used to compute the virial of the QM/MM particles.)
 //  *

 //  * if atom i shall be considered as MM,
 //  *   isMMatom[i] = true
 //  * UPDATE: store the shift in this array,
 //  * and change 'bool' to 'int':
 //  *   shiftMMatom[i] = the value of shift

    std::vector<int> shiftMMatom(nAtoms, -1); // ALL ATOMS IN SIMULATION - IS THAT NECESSARY???

    // init PBC
    gmx::IVec null_ivec;
    clear_ivec(null_ivec);
    t_pbc pbc;
//  set_pbc_dd(&pbc, pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : null_ivec, false, box);
    set_pbc_dd(&pbc, pbcType, haveDDAtomOrdering(*cr) ? &cr->dd->numCells : &null_ivec, false, box);

    // copy pointers
    QMMM_QMrec&                           qm_  = qm[0];
    QMMM_MMrec&                           mm_  = mm[0];
    gmx::ArrayRef<const gmx::NbnxnPairlistCpu> nbl = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local).cpuLists();
    int                                   nnbl = nbl.ssize();
//  const gmx::ArrayRef<const int>        atomIndices = nbv->pairSearch_->gridSet().atomIndices();
    const gmx::ArrayRef<const int>        atomIndices = nbv->getLocalAtomOrder();

    // QM shift array
    // !!! CHECK THIS !!!
    rvec dx;
    qm_.shiftQM[0] = gmx::xyzToShiftIndex(0, 0, 0);
    for (int i = 1; i < qm_.nrQMatoms; i++)
    {
        qm_.shiftQM[i] = pbc_dx_aiuc(&pbc,
                                     x[globalToLocalAtomMap[qm_.indexQM[0]]],
                                     x[globalToLocalAtomMap[qm_.indexQM[i]]],
                                     dx);
    }
 // for (int i = 0; i < qm->nrQMatoms; i++)
 // {
 //     printf("VERLET QM SHIFT [%d] = %d\n", i, qm->shiftQM[i]);
 // }

    /* The MM candidates are the non-QM atoms of the clusters that share a pair list entry
     * with a cluster holding a QM atom. Each atom found gets the shift of its nearest QM atom.
     *
     * The search runs in two passes, both on the OpenMP threads of mdrun:
     *   1. the pair lists (one per thread of the pair search) are scanned in parallel, and
     *      the non-QM atoms of the clusters next to a QM cluster are collected per list;
     *   2. the shifts of the atoms found are computed in parallel, atom by atom.
     * The QM test is a table lookup, not a loop over the QM atoms. The result (atoms in
     * the order of their global index, and their shifts) does not depend on the threads.
     */
    std::vector<char> isQMglobal(nAtoms, 0);
    for (int q = 0; q < qm_.nrQMatoms; q++)
    {
        isQMglobal[qm_.indexQM[q]] = 1;
    }
    const auto globalOf = [&](int localAtom) {
        return localAtom < 0 ? -1 : localToGlobalAtomMap[localAtom];
    };
    const auto clusterHasQM = [&](int cluster, int na) {
        for (int i = 0; i < na; i++)
        {
            const int g = globalOf(atomIndices[na * cluster + i]);
            if (g >= 0 && isQMglobal[g])
            {
                return true;
            }
        }
        return false;
    };
    std::vector<std::vector<int>> candidates(nnbl);
#pragma omp parallel for num_threads(std::min(nnbl, qmmm_omp::maxThreads())) schedule(dynamic)
    for (int inbl = 0; inbl < nnbl; inbl++)
    {
        try
        {
            const gmx::NbnxnPairlistCpu& list  = nbl[inbl];
            std::vector<int>&            found = candidates[inbl];
            const auto addCluster = [&](int cluster, int na) {
                for (int i = 0; i < na; i++)
                {
                    const int g = globalOf(atomIndices[na * cluster + i]);
                    if (g >= 0 && !isQMglobal[g])
                    {
                        found.push_back(g);
                    }
                }
            };
            // loop over CI clusters
            for (const auto& ciEntry : list.ci)
            {
                const bool qmInCi = clusterHasQM(ciEntry.ci, list.na_ci);
                // loop over the corresponding CJ clusters
                for (int cj = ciEntry.cj_ind_start; cj < ciEntry.cj_ind_end; cj++)
                {
                    const int  cjCluster = list.cj.cj(cj);
                    const bool qmInCj    = clusterHasQM(cjCluster, list.na_cj);
                    // a QM atom in cluster CI: the non-QM atoms of cluster CJ are MM candidates
                    if (qmInCi)
                    {
                        addCluster(cjCluster, list.na_cj);
                    }
                    // and vice versa
                    if (qmInCj)
                    {
                        addCluster(ciEntry.ci, list.na_ci);
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    std::vector<int> mmCandidates;
    for (const std::vector<int>& found : candidates)
    {
        for (int g : found)
        {
            if (shiftMMatom[g] == -1)
            {
                shiftMMatom[g] = 0; // marked, the shift follows
                mmCandidates.push_back(g);
            }
        }
    }
    // the shift that puts the MM atom next to its nearest QM atom
    const int numCandidates = gmx::ssize(mmCandidates);
#pragma omp parallel for num_threads(qmmm_omp::numThreads(numCandidates, 64)) schedule(static)
    for (int c = 0; c < numCandidates; c++)
    {
        const int globalAtom = mmCandidates[c];
        real      dist       = 1000.;
        int       sh         = -1;
        for (int q = 0; q < qm_.nrQMatoms; q++)
        {
            rvec bond;
            const int sh_t = pbc_dx_aiuc(&pbc,
                                         x[globalToLocalAtomMap[qm_.indexQM[q]]],
                                         x[globalToLocalAtomMap[globalAtom]],
                                         bond);
            if (norm(bond) < dist)
            {
                dist = norm(bond);
                sh   = sh_t;
            }
        }
        shiftMMatom[globalAtom] = sh;
    }

    // count the MM atoms found in the above search
    int nrMMatoms = 0;
    for (int i=0; i<nAtoms; i++) {
        // criterium for MM atom found
        if (shiftMMatom[i] != -1)
	    {
         // printf("VERLET MM ATOM %5d with shift %d\n", i, shiftMMatom[i]);
	        nrMMatoms++;
	    }
    }

 // printf("Number of potential MM atoms as found in the Verlet neighbor lists: %d\n", nrMMatoms);

    // allocate space and fill the array with atom numbers
    mm_.nrMMatoms_nbl = nrMMatoms;
    mm_.indexMM_nbl.resize(nrMMatoms);
    mm_.shiftMM_nbl.resize(nrMMatoms);
    // index i runs along isMMatom / shiftMMatom,
    //       j runs along the new indexMM_nbl array

    int count=0;
    for (int atom=0; atom<nAtoms; atom++) {
        // criterium for MM atom found
        if (shiftMMatom[atom] != -1)
	    {
	        mm_.indexMM_nbl[count] = atom;
	        mm_.shiftMM_nbl[count] = shiftMMatom[atom];
	        count++;
	    }
    }
    // check
    if (count != mm_.nrMMatoms_nbl) {
        printf("ERROR IN MM ATOM SEARCH -- VERLET BASED SCHEME\n");
	    exit(-1);
    }
} // update_QMMMrec_verlet_ns

namespace
{

/*! \brief Whether the per-atom QM/MM report files are written.
 *
 * On by default; GMX_QMMM_REPORTS set to 0, no, off or false switches off the
 * reports of both grompp and mdrun, together with the lines that point to them.
 */
bool qmmmReportsEnabled()
{
    const char* env = std::getenv("GMX_QMMM_REPORTS");
    if (env == nullptr)
    {
        return true;
    }
    for (const char* off : { "0", "no", "off", "false" })
    {
        if (gmx_strcasecmp(env, off) == 0)
        {
            return false;
        }
    }
    return true;
}

//! Name of a boundary charge scheme, as set with GMX_QMMM_POT_SCHEME
const char* potSchemeName(QMMM_rec::PotScheme scheme)
{
    switch (scheme)
    {
        case QMMM_rec::PotScheme::RC: return "RC";
        case QMMM_rec::PotScheme::RCD: return "RCD";
        case QMMM_rec::PotScheme::CS: return "CS";
        case QMMM_rec::PotScheme::Amber: return "AMBER";
        default: return "none";
    }
}

//! "RESnr NAME" label of a global atom
std::string qmmmGlobalAtomLabel(const gmx_mtop_t& mtop, int globalIndex)
{
    int         molb    = 0;
    int         resnr   = 0;
    const char* name    = nullptr;
    const char* resname = nullptr;
    mtopGetAtomAndResidueName(mtop, globalIndex, &molb, &name, &resnr, &resname, nullptr);
    return gmx::formatString("%s%d %s", resname, resnr, name);
}

//! Writes, atom by atom, the boundary treatment of the QM--MM electrostatics
void writeQmmmChargeReport(const gmx_mtop_t& mtop, const QMMM_rec& qr, const std::vector<bool>& bQM, const char* fileName)
{
    FILE* fp = std::fopen(fileName, "w");
    if (fp == nullptr)
    {
        fprintf(stderr, "WARNING: could not open %s for the QM/MM charge report\n", fileName);
        return;
    }
    std::vector<bool> isLA(bQM.size(), false);
    for (const QMMM_rec::LinkAtom& l : qr.linkAtoms)
    {
        isLA[l.la] = true;
    }
    const auto atomText = [&](int a) {
        return gmx::formatString("%7d %-14s %s", a + 1, qmmmGlobalAtomLabel(mtop, a).c_str(),
                                 isLA[a] ? "LA" : (bQM[a] ? "QM" : "MM"));
    };
    const auto chargeOf = [&mtop](int a) {
        int molb = 0;
        return mtopGetAtomParameters(mtop, a, &molb).q;
    };
    const QMMM_QMrec& qm_ = qr.qm[0];

    fprintf(fp, "; QM/MM electrostatics at the QM/MM boundary, written by gmx mdrun\n");
    fprintf(fp, "; atom numbers are global and 1-based, i.e. the numbering of the input .gro file\n");
    fprintf(fp, "; labels are RESIDUEnumber ATOMNAME from the topology; LA = link atom\n");
    fprintf(fp, "; boundary charge scheme: GMX_QMMM_POT_SCHEME = %s\n", potSchemeName(qr.potScheme));
    fprintf(fp, "; the external potential passed to DFTB+ and the QM/MM gradient use the same charges;\n");
    fprintf(fp, ";   the MM--MM interactions keep the charges of the topology\n\n");

    fprintf(fp, "[ qm_atoms ]\n; %d atoms, in the order of the QM group\n", qm_.nrQMatoms_get());
    for (int a = 0; a < static_cast<int>(bQM.size()); a++)
    {
        if (bQM[a])
        {
            fprintf(fp, "%s\n", atomText(a).c_str());
        }
    }

    fprintf(fp, "\n[ link_atoms ]\n; virtual sites of the QM group constructed from a QM1 and an MM1 atom: %zu\n",
            qr.linkAtoms.size());
    fprintf(fp, "; %-26s %-26s %s\n", "link atom", "QM1", "MM1");
    for (const QMMM_rec::LinkAtom& l : qr.linkAtoms)
    {
        fprintf(fp, "%s %s %s\n", atomText(l.la).c_str(), atomText(l.qm1).c_str(), atomText(l.mm1).c_str());
    }

    fprintf(fp, "\n[ removed_mm1_charges ]\n");
    if (qr.potScheme == QMMM_rec::PotScheme::None)
    {
        fprintf(fp, "; none: without a boundary charge scheme every MM charge enters the QM--MM electrostatics in full\n");
    }
    else
    {
        fprintf(fp, "; MM1 atoms whose charge is removed from the QM--MM electrostatics of every QM atom\n");
        fprintf(fp, "; %-26s %10s\n", "MM1 atom", "q_MM1");
        for (const QMMM_rec::LinkAtom& l : qr.linkAtoms)
        {
            fprintf(fp, "%s %+10.5f\n", atomText(l.mm1).c_str(), chargeOf(l.mm1));
        }
    }

    fprintf(fp, "\n[ charge_shift ]\n");
    if (qr.potScheme == QMMM_rec::PotScheme::Amber)
    {
        fprintf(fp, "; AMBER: the MM1 charges are spread over the other MM atoms of their molecule,\n");
        fprintf(fp, ";   in the QM--MM electrostatics only; every atom that receives a share:\n");
        fprintf(fp, "; %-26s %10s %14s\n", "MM atom", "q_MM", "added");
        for (size_t a = 0; a < qr.potChargeShift.size(); a++)
        {
            if (qr.potChargeShift[a] != real(0.0))
            {
                fprintf(fp, "%s %+10.5f %+14.8e\n", atomText(static_cast<int>(a)).c_str(),
                        chargeOf(static_cast<int>(a)), qr.potChargeShift[a]);
            }
        }
    }
    else
    {
        fprintf(fp, "; none\n");
    }

    fprintf(fp, "\n[ point_charges ]\n");
    if (qr.potPoints.empty())
    {
        fprintf(fp, "; none\n");
    }
    else
    {
        fprintf(fp, "; fictitious charges at x(MM1) + f * (x(MM2) - x(MM1)), seen by the QM atoms only;\n");
        fprintf(fp, ";   f = 1 is a change of the charge of MM2 itself. Their forces are passed to\n");
        fprintf(fp, ";   MM1 and MM2 with the weights (1 - f) and f, as for a two-atom virtual site\n");
        fprintf(fp, "; %-26s %-26s %6s %10s\n", "MM1 atom", "MM2 atom", "f", "charge");
        for (const QMMM_rec::PotPoint& pt : qr.potPoints)
        {
            fprintf(fp, "%s %s %6.3f %+10.5f\n", atomText(pt.a).c_str(), atomText(pt.b).c_str(), pt.f, pt.q);
        }
    }
    std::fclose(fp);
}

} // namespace

// Set up the boundary charge scheme of the QM--MM electrostatics (GMX_QMMM_POT_SCHEME).
//
// With a scheme, the charge of every MM1 atom (the MM atom from which a link atom is
// constructed) is removed from the QM--MM electrostatics of every QM atom, and replaced by
// fictitious point charges near the MM1--MM2 bonds (Lin and Truhlar, J. Phys. Chem. A 109,
// 3991 (2005) for RC and RCD; Sherwood et al., THEOCHEM 632, 1 (2003) for CS), or spread
// evenly over the other MM atoms of the molecule (AMBER). The external potential of DFTB+ and
// the QM/MM gradient use these same charges, so the forces stay the gradient of the energy.
// The topology and the MM--MM interactions keep the charges of the force field.
void QMMM_rec::init_QMMM_boundary(const gmx_mtop_t* mtop)
{
    QMMM_QMrec& qm_ = qm[0];

    const char* env = getenv("GMX_QMMM_POT_SCHEME");
    potScheme       = PotScheme::None;
    if (env != nullptr)
    {
        if (gmx_strcasecmp(env, "none") == 0)
        {
            potScheme = PotScheme::None;
        }
        else if (gmx_strcasecmp(env, "RC") == 0)
        {
            potScheme = PotScheme::RC;
        }
        else if (gmx_strcasecmp(env, "RCD") == 0)
        {
            potScheme = PotScheme::RCD;
        }
        else if (gmx_strcasecmp(env, "CS") == 0)
        {
            potScheme = PotScheme::CS;
        }
        else if (gmx_strcasecmp(env, "AMBER") == 0)
        {
            potScheme = PotScheme::Amber;
        }
        else
        {
            gmx_fatal(FARGS, "GMX_QMMM_POT_SCHEME must be none, RC, RCD, CS or AMBER, but it is '%s'.", env);
        }
    }

    const int         natoms = mtop->natoms;
    std::vector<bool> bQM(natoms, false);
    for (int j = 0; j < qm_.nrQMatoms; j++)
    {
        bQM[qm_.indexQM[j]] = true;
    }

    // Link atoms: QM virtual sites constructed from one QM and one MM atom.
    linkAtoms.clear();
    std::vector<bool> isLA(natoms, false);
    int               atomOffset = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& molt = mtop->moltype[molb.type];
        for (int mol = 0; mol < molb.nmol; mol++, atomOffset += molt.atoms.nr)
        {
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                if (!IS_VSITE(ftype) || ftype == InteractionFunction::VirtualSiteN)
                {
                    continue;
                }
                const int              nral = NRAL(ftype);
                const InteractionList& il   = molt.ilist[ftype];
                for (int i = 0; i < il.size(); i += 1 + nral)
                {
                    const int site = atomOffset + il.iatoms[i + 1];
                    if (!bQM[site])
                    {
                        continue;
                    }
                    std::vector<int> qmBuild, mmBuild;
                    for (int c = 2; c <= nral; c++)
                    {
                        const int a = atomOffset + il.iatoms[i + c];
                        (bQM[a] ? qmBuild : mmBuild).push_back(a);
                    }
                    if (mmBuild.empty())
                    {
                        continue; // a virtual site inside the QM region
                    }
                    if (qmBuild.size() != 1 || mmBuild.size() != 1)
                    {
                        gmx_fatal(FARGS,
                                  "QM virtual site %d is constructed from MM atoms, but it is not a "
                                  "link atom constructed from one QM and one MM atom.",
                                  site + 1);
                    }
                    linkAtoms.push_back({ site, qmBuild[0], mmBuild[0] });
                    isLA[site] = true;
                }
            }
        }
    }

    // Chemical bonds of the entire system, in global atom numbering, without the link
    //   atoms. The bonds inside the QM region are connections by now, which are chemical
    //   bonds as well.
    std::vector<std::vector<int>> bonds(natoms);
    atomOffset = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& molt = mtop->moltype[molb.type];
        for (int mol = 0; mol < molb.nmol; mol++, atomOffset += molt.atoms.nr)
        {
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                if (!IS_CHEMBOND(ftype))
                {
                    continue;
                }
                const InteractionList& il = molt.ilist[ftype];
                for (int i = 0; i < il.size(); i += 3)
                {
                    const int a = atomOffset + il.iatoms[i + 1];
                    const int b = atomOffset + il.iatoms[i + 2];
                    if (isLA[a] || isLA[b])
                    {
                        continue;
                    }
                    bonds[a].push_back(b);
                    bonds[b].push_back(a);
                }
            }
        }
    }

    const auto chargeOf = [mtop](int a) {
        int molb = 0;
        return mtopGetAtomParameters(*mtop, a, &molb).q;
    };

    potPoints.clear();
    potChargeShift.clear();
    isRemovedMM1.clear();
    if (potScheme != PotScheme::None)
    {
        if (linkAtoms.empty())
        {
            fprintf(stdout,
                    "NOTE: GMX_QMMM_POT_SCHEME=%s, but the QM region has no link atoms: there is no "
                    "MM1 charge to redistribute.\n",
                    potSchemeName(potScheme));
        }
        isRemovedMM1.assign(natoms, false);
        for (const LinkAtom& l : linkAtoms)
        {
            if (isRemovedMM1[l.mm1])
            {
                gmx_fatal(FARGS,
                          "GMX_QMMM_POT_SCHEME=%s: MM atom %d is the MM1 atom of more than one link "
                          "atom. This is not supported.",
                          potSchemeName(potScheme), l.mm1 + 1);
            }
            isRemovedMM1[l.mm1] = true;
        }
    }

    if (potScheme == PotScheme::Amber)
    {
        // As in AMBER: the MM1 charges are removed, and their sum is added evenly to the
        //   other MM atoms of the same molecule.
        std::vector<int> molStart(natoms), molSize(natoms);
        int              start = 0;
        for (const gmx_molblock_t& molb : mtop->molblock)
        {
            const int nat = mtop->moltype[molb.type].atoms.nr;
            for (int mol = 0; mol < molb.nmol; mol++, start += nat)
            {
                std::fill(molStart.begin() + start, molStart.begin() + start + nat, start);
                std::fill(molSize.begin() + start, molSize.begin() + start + nat, nat);
            }
        }
        std::map<int, double> sumOfMolecule; // first atom of the molecule -> sum of its MM1 charges
        for (const LinkAtom& l : linkAtoms)
        {
            sumOfMolecule[molStart[l.mm1]] += chargeOf(l.mm1);
        }
        potChargeShift.assign(natoms, real(0.0));
        for (const auto& mol : sumOfMolecule)
        {
            std::vector<int> receivers;
            for (int a = mol.first; a < mol.first + molSize[mol.first]; a++)
            {
                if (!bQM[a] && !isRemovedMM1[a])
                {
                    receivers.push_back(a);
                }
            }
            if (receivers.empty())
            {
                gmx_fatal(FARGS,
                          "GMX_QMMM_POT_SCHEME=AMBER: the molecule of atoms %d-%d has no MM atom to "
                          "spread the MM1 charges to.",
                          mol.first + 1, mol.first + molSize[mol.first]);
            }
            const real shift = static_cast<real>(mol.second / receivers.size());
            for (int a : receivers)
            {
                potChargeShift[a] = shift;
            }
            fprintf(stdout,
                    "QM/MM electrostatics with the boundary charge scheme AMBER: molecule of atoms "
                    "%d-%d, MM1 charges (sum %+.5f) removed and spread over its %zu other MM atoms, "
                    "%+.5e each.\n",
                    mol.first + 1, mol.first + molSize[mol.first], mol.second, receivers.size(), shift);
        }
    }
    else if (potScheme != PotScheme::None)
    {
        for (const LinkAtom& l : linkAtoms)
        {
            std::vector<int> mm2;
            for (int b : bonds[l.mm1])
            {
                if (!bQM[b])
                {
                    mm2.push_back(b);
                }
            }
            if (mm2.empty())
            {
                gmx_fatal(FARGS,
                          "GMX_QMMM_POT_SCHEME=%s: the MM1 atom %d of link atom %d has no MM2 atom to "
                          "redistribute its charge to.",
                          potSchemeName(potScheme), l.mm1 + 1, l.la + 1);
            }
            for (int b : mm2)
            {
                if (isRemovedMM1[b])
                {
                    gmx_fatal(FARGS,
                              "GMX_QMMM_POT_SCHEME=%s: atom %d is an MM2 atom of MM1 atom %d and an MM1 "
                              "atom itself. This is not supported.",
                              potSchemeName(potScheme), b + 1, l.mm1 + 1);
                }
                for (int c : bonds[b])
                {
                    if (bQM[c])
                    {
                        gmx_fatal(FARGS,
                                  "GMX_QMMM_POT_SCHEME=%s: the MM2 atom %d of MM1 atom %d is bonded to "
                                  "the QM atom %d. This is not supported.",
                                  potSchemeName(potScheme), b + 1, l.mm1 + 1, c + 1);
                    }
                }
            }

            const real q0 = chargeOf(l.mm1) / mm2.size();
            for (int b : mm2)
            {
                switch (potScheme)
                {
                    case PotScheme::RC: potPoints.push_back({ l.mm1, b, real(0.5), q0 }); break;
                    case PotScheme::RCD:
                        potPoints.push_back({ l.mm1, b, real(0.5), real(2.0) * q0 });
                        potPoints.push_back({ l.mm1, b, real(1.0), -q0 });
                        break;
                    case PotScheme::CS:
                    {
                        /* Charge shift: q0 moves from MM1 onto MM2, which changes the dipole of
                         * the MM1--MM2 bond by +q0*b (b = MM1->MM2). The compensating pair sits
                         * at the fractions fMinus and fPlus of the bond and has the dipole
                         * qPair*(fMinus - fPlus)*b, so qPair = q0 / (fPlus - fMinus) cancels it.
                         */
                        constexpr real fMinus = 0.94;
                        constexpr real fPlus  = 1.06;
                        const real     qPair  = q0 / (fPlus - fMinus);
                        potPoints.push_back({ l.mm1, b, real(1.0), q0 });
                        potPoints.push_back({ l.mm1, b, fMinus, qPair });
                        potPoints.push_back({ l.mm1, b, fPlus, -qPair });
                        break;
                    }
                    default: break;
                }
            }
        }
    }

    if (potScheme == PotScheme::None)
    {
        fprintf(stdout,
                "QM/MM electrostatics without a boundary charge scheme -- every MM atom within the "
                "cut-off interacts with the QM atoms with its full charge (%zu link atoms).\n"
                "To change, set environment variable GMX_QMMM_POT_SCHEME to RC, RCD, CS or AMBER.\n",
                linkAtoms.size());
    }
    else
    {
        fprintf(stdout,
                "QM/MM electrostatics with the boundary charge scheme %s: %zu MM1 charges removed, "
                "%zu fictitious point charges added%s.\n"
                "  The external potential passed to DFTB+ and the QM/MM forces are built from the "
                "same charges;\n  the forces on the point charges are passed to their MM1 and MM2 "
                "atoms as for two-atom virtual sites.\n  The MM--MM interactions keep the charges of "
                "the topology.\n",
                potSchemeName(potScheme), linkAtoms.size(), potPoints.size(),
                potScheme == PotScheme::Amber ? " (charges spread over the molecule instead)" : "");
    }

    if (qmmmReportsEnabled())
    {
        const char* reportFile = getenv("GMX_QMMM_EXCLUSION_REPORT");
        if (reportFile == nullptr)
        {
            reportFile = "qmmm_exclusion_report.txt";
        }
        writeQmmmChargeReport(*mtop, *this, bQM, reportFile);
        fprintf(stdout,
                "The QM/MM boundary charge scheme is listed atom by atom in %s\n"
                "  (file name set with GMX_QMMM_EXCLUSION_REPORT, switched off with GMX_QMMM_REPORTS=off).\n",
                reportFile);
    }
}

// The charges of the QM--MM electrostatics on the current short-range list. Has to be
//   redone whenever that list changes, i.e. in every step.
void QMMM_rec::update_QMMM_boundary_SR()
{
    QMMM_MMrec& mm_ = mm[0];

    if (static_cast<int>(localIndexOfAtom.size()) != nAtoms)
    {
        localIndexOfAtom.assign(nAtoms, -1);
    }
    else
    {
        std::fill(localIndexOfAtom.begin(), localIndexOfAtom.end(), -1);
    }
    qmmmChargesSR.resize(mm_.nrMMatoms);
    qmmmScaleSR.resize(mm_.nrMMatoms);
    for (int k = 0; k < mm_.nrMMatoms; k++)
    {
        const int a         = mm_.indexMM[k];
        localIndexOfAtom[a] = k;
        qmmmChargesSR[k]    = mm_.MMcharges[k];
        if (!potChargeShift.empty())
        {
            qmmmChargesSR[k] += potChargeShift[a] * mm_.scalefactor;
        }
        qmmmScaleSR[k] = (!isRemovedMM1.empty() && isRemovedMM1[a]) ? real(0.0) : real(1.0);
    }

    // With PME, the reciprocal-space part of a removed MM1 charge is on the grid and can only
    //   be taken out through the short-range list, see calculate_SR_QM_MM().
    if (!isRemovedMM1.empty())
    {
        for (const LinkAtom& l : linkAtoms)
        {
            if (localIndexOfAtom[l.mm1] < 0)
            {
                gmx_fatal(FARGS,
                          "QM/MM boundary charge scheme: the MM1 atom %d is not on the QM/MM "
                          "short-range list. Increase rcoulomb.",
                          l.mm1 + 1);
            }
        }
    }
}

// The charges of the QM--MM electrostatics on the full MM list (PME).
void QMMM_rec::update_QMMM_boundary_full()
{
    QMMM_MMrec& mm_ = mm[0];
    qmmmChargesFull.resize(mm_.nrMMatoms_full);
    for (int k = 0; k < mm_.nrMMatoms_full; k++)
    {
        qmmmChargesFull[k] = mm_.MMcharges_full[k];
        if (!potChargeShift.empty())
        {
            qmmmChargesFull[k] += potChargeShift[mm_.indexMM_full[k]] * mm_.scalefactor;
        }
    }
}

real QMMM_rec::calculate_QMMM(// const t_commrec*      cr,
                              gmx::ForceWithVirial* forceWithVirial,
                              t_nrnb*               nrnb,
                              gmx_wallcycle*        wcycle)
{
    if (!GMX_QMMM)
    {
        gmx_incons("Compiled without QMMM");
    }

    real QMener = 0.0;
    // A selection for the QM package depending on which is requested
    // (Gaussian, GAMESS-UK, MOPAC or ORCA) needs to be implemented here.
    // Now it works through defines.
    //   ... Not so nice yet

    QMMM_QMrec* qm_ = &(qm[0]);
    QMMM_MMrec* mm_ = &(mm[0]);

    rvec *forces = nullptr,
         *fshift = nullptr;
 
 // gmx::ArrayRef<gmx::RVec> fMM      = forceWithVirial->force_.data();
                  gmx::RVec *fMM      = forceWithVirial->force_.data();
 // gmx::ArrayRef<gmx::RVec> fshiftMM = forceWithShiftForces->shiftForces();

    if (GMX_QMMM_DFTBPLUS)
    {
        snew(forces, (qm_->nrQMatoms + mm_->nrMMatoms + mm_->nrMMatoms_full));
     // snew(fshift, (qm_->nrQMatoms + mm_->nrMMatoms + mm_->nrMMatoms_full));
    }
    else
    {
        snew(forces, (qm_->nrQMatoms + mm_->nrMMatoms));
     // snew(fshift, (qm_.nrQMatoms + mm_.nrMMatoms));
    }

    computeVirial = forceWithVirial->computeVirial_;
    clear_mat(recipVirialCorrection);
    QMener = call_QMroutine(this, qm_, mm_, forces, fshift, nrnb, wcycle); // (cr,)

    if (GMX_QMMM_DFTBPLUS)
    {
        for (int i = 0; i < qm_->nrQMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[globalToLocalAtomMap[qm_->indexQM[i]]][j]        -= forces[i][j];
             // fshiftMM[globalToLocalAtomMap[qm_->shiftQM[i]]][j]   += fshift[i][j];
            }
         // printf("F[%5d] = %8.2f %8.2f %8.2f\n", qm_->indexQM[i], forces[i][0], forces[i][1], forces[i][2]);
        }
        // the MM atoms are distinct, so the threads write distinct forces
#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_->nrMMatoms, 1024)) schedule(static)
        for (int i = 0; i < mm_->nrMMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[globalToLocalAtomMap[mm_->indexMM[i]]][j]        -= forces[qm_->nrQMatoms+i][j];
             // fshiftMM[globalToLocalAtomMap[mm_->shiftMM[i]]][j]   += fshift[qm_->nrQMatoms+i][j];
            }
         // if (i<30) if (norm(forces[qm_->nrQMatoms+i]) > 10.)
         //   printf("F_MM[%5d] = %8.2f %8.2f %8.2f\n", mm_->indexMM[i],
         //     forces[qm_->nrQMatoms+i][0], forces[qm_->nrQMatoms+i][1], forces[qm_->nrQMatoms+i][2]);
        }
#pragma omp parallel for num_threads(qmmm_omp::numThreads(mm_->nrMMatoms_full, 1024)) schedule(static)
        for (int i = 0; i < mm_->nrMMatoms_full; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[globalToLocalAtomMap[mm_->indexMM_full[i]]][j]        -= forces[qm_->nrQMatoms+mm_->nrMMatoms+i][j];
             // fshiftMM[globalToLocalAtomMap[mm_->shiftMM_full[i]]][j]   += fshift[qm_->nrQMatoms+mm_->nrMMatoms+i][j];
            }
         // if (i<100) if (norm(forces[qm_->nrQMatoms+mm_->nrMMatoms+i]) > 10.)
         //   printf("F_MM_F[%5d] = %8.2f %8.2f %8.2f\n", mm_->indexMM_full[i],
         //   forces[qm_->nrQMatoms+mm_->nrMMatoms+i][0], forces[qm_->nrQMatoms+mm_->nrMMatoms+i][1], forces[qm_->nrQMatoms+mm_->nrMMatoms+i][2]);
        }
    }
    else
    {
        for (int i = 0; i < qm_->nrQMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[globalToLocalAtomMap[qm_->indexQM[i]]][j]          -= forces[i][j];
             // fshiftMM[globalToLocalAtomMap[qm_->shiftQM[i]]][j]     += fshift[i][j];
            }
        }
        for (int i = 0; i < mm_->nrMMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[globalToLocalAtomMap[mm_->indexMM[i]]][j]      -= forces[qm_->nrQMatoms+i][j];
             // fshiftMM[globalToLocalAtomMap[mm_->shiftMM[i]]][j] += fshift[qm_->nrQMatoms+i][j];
            }
        }
    }

    /* Virial of the QM/MM forces. They are collected in a buffer of their own (see
     * init_forcerec(), which marks a QM/MM run as having direct virial contributions), so the
     * single sum over the shift forces does not see them and the virial is supplied here.
     * Every force is paired with the position it was computed from: the QM and MM atoms are
     * taken as the periodic images nearest to the first QM atom, which is what the
     * minimum-image QM--MM terms and the contiguous QM cluster of DFTB+ use. That is exact for
     * everything computed in real space, the forces of the fictitious charges of the boundary
     * scheme included (they are linear in the positions of MM1 and MM2); the reciprocal-space
     * part of PME is replaced by its exact virial, which gradient_QM_MM() has prepared.
     * Rectangular boxes only, like the rest of this interface.
     */
    if (GMX_QMMM_DFTBPLUS && forceWithVirial->computeVirial_)
    {
        const auto nearestImage = [qm_](const rvec x, rvec image) {
            for (int d = 0; d < DIM; d++)
            {
                const real L  = qm_->box[d][d];
                real       dx = x[d] - qm_->xQM[0][d];
                if (L > 0)
                {
                    dx -= L * std::round(dx / L);
                }
                image[d] = qm_->xQM[0][d] + dx;
            }
        };
        const auto addTerm = [](matrix w, const rvec x, const rvec gradient) {
            // the force is minus the stored gradient
            for (int a = 0; a < DIM; a++)
            {
                for (int b = 0; b < DIM; b++)
                {
                    w[a][b] -= x[a] * gradient[b];
                }
            }
        };
        matrix w;
        clear_mat(w);
        rvec   image;
        for (int i = 0; i < qm_->nrQMatoms; i++)
        {
            nearestImage(qm_->xQM[i], image);
            addTerm(w, image, forces[i]);
        }
        for (int i = 0; i < mm_->nrMMatoms; i++)
        {
            nearestImage(mm_->xMM[i].as_vec(), image);
            addTerm(w, image, forces[qm_->nrQMatoms + i]);
        }
        {
            // the long list of all MM atoms by thread, added in the order of the threads
            struct VirialPart
            {
                matrix m;
            };
            const int               nrFull = mm_->nrMMatoms_full;
            const int               nth    = qmmm_omp::numThreads(nrFull, 1024);
            std::vector<VirialPart> wThread(nth);
#pragma omp parallel for num_threads(nth) schedule(static)
            for (int t = 0; t < nth; t++)
            {
                clear_mat(wThread[t].m);
                rvec imageT;
                for (int i = qmmm_omp::chunkStart(nrFull, t, nth); i < qmmm_omp::chunkStart(nrFull, t + 1, nth); i++)
                {
                    nearestImage(mm_->xMM_full[i].as_vec(), imageT);
                    addTerm(wThread[t].m, imageT, forces[qm_->nrQMatoms + mm_->nrMMatoms + i]);
                }
            }
            for (int t = 0; t < nth; t++)
            {
                m_add(w, wThread[t].m, w);
            }
        }
        matrix virial;
        msmul(w, -0.5, virial);
        m_add(virial, recipVirialCorrection, virial);
        forceWithVirial->addVirialContribution(virial);
    }

    sfree(forces);
 // sfree(fshift);

    return QMener;
} // calculate_QMMM

#pragma GCC diagnostic pop
