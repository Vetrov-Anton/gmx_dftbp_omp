/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
//#include "gromacs/mdlib/qmmm.h"

#ifndef GMX_MDLIB_QM_DFTBPLUS_H
#define GMX_MDLIB_QM_DFTBPLUS_H

#include "config.h"

#define HARTREE_TO_EV     (27.211396132)

void
init_dftbplus(QMMM_QMrec*       qm,
              QMMM_rec*         qr,
           // const t_forcerec* fr,
              const t_inputrec* ir,
              const t_commrec*  cr);
           // gmx_wallcycle*    wcycle);

real
call_dftbplus(QMMM_rec*         qr,
          //  const t_commrec*  cr,
              QMMM_QMrec*       qm,
              const QMMM_MMrec& mm,
              rvec              f[],
              rvec              fshift[],
              t_nrnb*           nrnb,
              gmx_wallcycle*    wcycle);

struct Context;

void initialize_context(Context*          cont,
                        int               nrQMatoms,
                        int               qmmm_variant,
                        QMMM_rec*         qr_in,
                     // const t_forcerec* fr_in,
                        const t_inputrec* ir_in,
                        const t_commrec*  cr_in);
                     // gmx_wallcycle*    wcycle_in);

void calcQMextPotPME(Context *cont, double *q, double *extpot);

extern "C" void calcqmextpot(void *refptr, double *q, double *extpot);

extern "C" void calcqmextpotgrad(void *refptr, gmx_unused double *q, double *extpotgrad);


// The C API of DFTB+, from the installed header of the DFTB+ that GROMACS is linked
// against (see cmake/gmxManageDftbPlus.cmake). Its content differs between the
// versions of DFTB+; GMX_QMMM_DFTBPLUS_ATOM_LIST and GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS
// tell which variant is in use.
#if GMX_QMMM_DFTBPLUS
// Some versions of dftbplus.h use the C99 keyword _Bool, which is not valid C++.
#    define _Bool bool
#    include <dftbplus.h>
#    undef _Bool
#endif

#endif
