/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "config.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

#if GMX_QMMM_DFTBPLUS

#include <dlfcn.h>
#include <chrono>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/ewald/pme.h"
//#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/qmmm_threading.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

//#include "dftbplus_gromacs.h"
#include "gromacs/mdlib/qm_dftbplus.h"

// TODO -- remove t_commrec* cr from all functions, since it is not needed in gmx_pme_do, any longer
typedef struct Context {
  bool              pme;
  int               n;
  const t_commrec*  cr;
  QMMM_rec*         qr;
//const t_forcerec* fr;
  t_nrnb*           nrnb;
  gmx_wallcycle*    wcycle;
  real              rcoul;
  real              ewaldcoeff_q;
} Context;

// extern "C" {void readdftbplusinput();}
// extern "C" {void performdftbpluscalculation(real *e, int *n, real x[], real q[], real extshift[], real f[]);}

/* Fill up the context (status structure) with all the relevant data
 */
void initialize_context(Context*          cont,
                        int               nrQMatoms,
                        int               qmmm_variant,
                        QMMM_rec*         qr_in,
                     // const t_forcerec* fr_in,
                        const t_inputrec* ir_in,
                        const t_commrec*  cr_in)
                     // gmx_wallcycle*    wcycle_in)
//                      const real        rcoul_in,
//                      const real        ewaldcoeff_q_in)
{
  /* The "cr" and "qr" structures will be initialized
   *   at the start of every MD step,
   *   and will be remembered for all of the SCC iterations.
   * That way, these data need not pass through the DFTB+ program.
   *
   * NOTE: In the new "C++" implementation,
   *       this is called from the constructor QMMM_rec::QMMM_rec(),
   *       so the object qr is not available yet!
   */

  if (qmmm_variant == eqmmmPME) {
      cont->pme = true;
  } else {
      cont->pme = false;
  }
  printf("qmmm_variant = %d\n", qmmm_variant);
  printf("cont->pme = %s\n", cont->pme ? "true" : "false");
//cont->n = fr_in->qr->qm[0].nrQMatoms; // object qr not available yet!
  cont->n = nrQMatoms;
  printf("cont->n = %d\n", cont->n);
  if (cont->pme)
  {
      cont->cr           = cr_in;
      cont->qr           = qr_in;
   // cont->wcycle       = wcycle_in;
      cont->rcoul        = ir_in->rcoulomb;
      cont->ewaldcoeff_q = calc_ewaldcoeff_q(ir_in->rcoulomb, ir_in->ewald_rtol);
      printf("cont->cr = %p\n", static_cast<const void*>(cont->cr));
      printf("cont->qr = %p\n", static_cast<const void*>(cont->qr));
   // printf("&(cont->qr->pmedata) = %p\n", static_cast<void*>(&(cont->qr->pmedata)));
   // printf("cont->qr->pmedata = %p\n", static_cast<void*>(cont->qr->pmedata));
   // printf("cont->qr->pmedata = %p\n", static_cast<void*>(cont->qr->pmedata.get())); UNCOMMENT!
      printf("cont->rcoul = %f\n", cont->rcoul);
      printf("cont->ewaldcoeff_q = %f\n", cont->ewaldcoeff_q);
  }

  return;
}

/* Wall-time breakdown of the QM/MM step, printed every GMX_QMMM_TIMING steps (in ms/step):
 *   the potential of the MM atoms on the QM atoms (real space and PME), the DFTB+ calls with
 *   the PME potential of the periodic QM images computed in the SCC iterations (callback),
 *   and the QM/MM gradients.
 */
struct QmmmTiming
{
    using Clock = std::chrono::steady_clock;
    double mmPotential = 0, dftb = 0, imageCallback = 0, gradient = 0;
    int    numCallbacks = 0, numSteps = 0;
    static double since(Clock::time_point t0)
    {
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }
};
static QmmmTiming qmmmTiming;

/* Calculate the external potential due to periodic images of QM atoms with PME.
 */
void calcQMextPotPME(Context *cont, double *q, double *extpot)
{
//int n = cont->fr->qr->qm[0]->nrQMatoms;
  int n = cont->n;

//printf("calcQMextPotPME: %s\n", setup ? "setup" : "calculation");

  if (cont->pme)
  {
      const auto t0 = QmmmTiming::Clock::now();
      /* PERFORM THE REAL CALCULATION */
      real *extpot_real;
      snew(extpot_real, n);
      for (int i=0; i<n; i++)
      {
          cont->qr->qm[0].QMcharges_set(i, (real) -q[i]); // check sign TODO
      }
      cont->qr->calculate_complete_QM_QM(cont->nrnb, cont->wcycle, extpot_real); // cont->cr ... *cont->qr->pmedata, extpot_real);
      for (int i=0; i<n; i++)
      {
          extpot[i] = (double) - extpot_real[i]; // sign OK
      }
      sfree(extpot_real);
      qmmmTiming.imageCallback += QmmmTiming::since(t0);
      qmmmTiming.numCallbacks++;
  }
  else
  {
      for (int i=0; i<n; i++)
      {
          extpot[i] = 0.;
      }
  }

  return;
}

/* The wrapper for the calculation of external potential
 *   due to the periodic images of QM atoms with PME.
 */
extern "C" void calcqmextpot(void *refptr, double *q, double *extpot)
{
  Context *cont = (Context *) refptr;
  calcQMextPotPME(cont, q, extpot);
  return;
}

/* The wrapper for the calculation of grdient of external potential
 *   due to the periodic images of QM atoms with PME.
 * This is function does not calculate anything,
 *   because it is not needed in this DFTB+Gromacs implementation.
 */
extern "C" void calcqmextpotgrad(void *refptr, gmx_unused double *q, double *extpotgrad)
{
  Context *cont = (Context *) refptr;
  for (int i=0; i<3*cont->n; i++)
      extpotgrad[i] = 0.;
  return;
}

/* The OpenMP threads of DFTB+.
 *   DFTB+ runs in the process of mdrun, in its OpenMP thread pool: by default on all of the
 *   threads of mdrun (-ntomp), the threads that mdrun has pinned. GMX_QMMM_DFTB_NTHREADS
 *   sets another number, e.g. fewer threads for a small QM zone.
 */
static int dftbNumThreads()
{
    static int numThreads = -1;
    if (numThreads < 0)
    {
        numThreads = qmmm_omp::maxThreads();
        if (const char* env = getenv("GMX_QMMM_DFTB_NTHREADS"))
        {
            const int value = atoi(env);
            if (value > 0)
            {
                numThreads = value;
            }
            else
            {
                printf("NOTE: GMX_QMMM_DFTB_NTHREADS=%s is not a positive number; DFTB+ uses %d threads.\n",
                       env, numThreads);
            }
        }
    }
    return numThreads;
}

/* The threading of the BLAS/LAPACK library of DFTB+, if it is OpenBLAS (looked up at run time,
 *   so any BLAS works). OpenBLAS built on OpenMP shares the thread pool of mdrun and DFTB+;
 *   the pthreads build has a pool of its own, whose threads are not pinned by mdrun. Its
 *   number of threads is then set to that of DFTB+.
 */
static void setUpDftbBlasThreads(int numThreads)
{
    using GetInt    = int (*)();
    using GetString = char* (*)();
    using SetInt    = void (*)(int);
    const auto getParallel = reinterpret_cast<GetInt>(dlsym(RTLD_DEFAULT, "openblas_get_parallel"));
    const auto getConfig   = reinterpret_cast<GetString>(dlsym(RTLD_DEFAULT, "openblas_get_config"));
    const auto setThreads  = reinterpret_cast<SetInt>(dlsym(RTLD_DEFAULT, "openblas_set_num_threads"));
    if (getParallel == nullptr)
    {
        printf("QM/MM threads: the BLAS/LAPACK of DFTB+ is not OpenBLAS; its threading is not controlled by mdrun.\n");
        return;
    }
    const int   mode     = getParallel(); // 0 sequential, 1 pthreads, 2 OpenMP
    const char* modeName = mode == 0 ? "sequential" : (mode == 1 ? "pthreads" : "OpenMP");
    printf("QM/MM threads: BLAS/LAPACK of DFTB+ is OpenBLAS (%s), %s.\n",
           getConfig != nullptr ? getConfig() : "unknown version", modeName);
    if (mode == 1)
    {
        if (setThreads != nullptr)
        {
            setThreads(numThreads);
        }
        printf("NOTE: OpenBLAS uses a pthreads pool of its own (%d threads), which mdrun does not pin.\n"
               "      The OpenMP build of OpenBLAS (Debian/Ubuntu: libopenblas0-openmp) runs on the\n"
               "      threads of mdrun instead and is recommended.\n",
               numThreads);
    }
}

/* Runs the DFTB+ calls of one MD step on numThreads OpenMP threads, and restores the number of
 *   threads of mdrun afterwards (gmx_omp_nthreads_init() has set it as the default).
 */
class DftbThreadScope
{
public:
    explicit DftbThreadScope(int numThreads) : previous_(gmx_omp_get_max_threads())
    {
        if (numThreads != previous_)
        {
            gmx_omp_set_num_threads(numThreads);
        }
    }
    ~DftbThreadScope()
    {
        if (gmx_omp_get_max_threads() != previous_)
        {
            gmx_omp_set_num_threads(previous_);
        }
    }

private:
    int previous_;
};

/* DFTBPLUS interface routines */

void init_dftbplus(QMMM_QMrec*       qm,
                   QMMM_rec*         qr,
                // const t_forcerec* fr,
                   const t_inputrec* ir,
                   const t_commrec*  cr)
                // gmx_wallcycle*    wcycle)
//void init_dftbplus(t_forcerec *fr)
{
    /* perhaps check the geometry first, to see which elements we have? */

    // variables
    int nAtom;
    int nSpecies;
    char ptrElement[20][3]; // element names may have up to 2 characters (+1 null termination)
    int *ptrSpecies;
    int *atomicNumber;
    int atomicNumberBySpecies[20];
    char chemSymbol[55][3] = { "XX",
          "H",                                      "He",
          "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
          "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
          "K",  "Ca",
          "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                      "Ga", "Ge", "As", "Se", "Br", "Kr",
          "Rb", "Sr",
          "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                      "In", "Sn", "Sb", "Te", "I",  "Xe" };

    static DftbPlus  calculator;
    DftbPlusInput    input;
#if GMX_QMMM_DFTBPLUS_ATOM_LIST
    DftbPlusAtomList atomList;
#endif

    /* This structure will be passed through DFTB+
     *   into the Gromacs calculator calcQMextPot
     */
  //static Context *cont;
  //snew(cont, 1);

  //initialize_context(cont, qm->nrQMatoms, qm->qmmm_variant, fr, ir, cr, nrnb, wcycle);
    snew(qm->dftbContext, 1);
    initialize_context(qm->dftbContext, qm->nrQMatoms_get(), qm->qmmm_variant_get(), qr, ir, cr); //, wcycle);
  //qm->dftbContext = cont;

    snew(qm->dpcalc, 1);

    /* Fill up the context with all the relevant data! */

    /* The layout of the parallel run: one MPI rank, OpenMP threads for everything. */
    {
        const int numThreadsDftb = dftbNumThreads();
        printf("QM/MM threads: %d MPI rank(s); MM forces and the QM/MM electrostatics on %d OpenMP "
               "thread(s), PME on %d, DFTB+ on %d.\n",
               cr->commMySim.size(), qmmm_omp::maxThreads(),
               gmx_omp_nthreads_get(ModuleMultiThread::Pme), numThreadsDftb);
        if (numThreadsDftb > qmmm_omp::maxThreads())
        {
            printf("NOTE: DFTB+ uses more threads (%d) than mdrun (%d); the extra threads are not pinned.\n",
                   numThreadsDftb, qmmm_omp::maxThreads());
        }
        setUpDftbBlasThreads(numThreadsDftb);
        fflush(stdout);
    }

    /* Initialize the DFTB+ calculator, on the threads it will run on */
    DftbThreadScope threadScope(dftbNumThreads());
    dftbp_init(&calculator, "dftb_in.out");
    printf("DFTB+ calculator has been created!\n");
    {
        int apiMajor = 0, apiMinor = 0, apiPatch = 0;
        dftbp_api(&apiMajor, &apiMinor, &apiPatch);
        printf("Linked DFTB+ C API version %d.%d.%d; QM atoms and species are taken %s.\n",
               apiMajor, apiMinor, apiPatch,
               GMX_QMMM_DFTBPLUS_ATOM_LIST ? "from the GROMACS topology"
                                           : "from the Geometry block of dftb_in.hsd");
    }

    /* Parse the input file and store the input-tree */
    dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
    printf("DFTB+ input has been read!\n");

    /* Pass the list of QM atoms to DFTB+ */
    nAtom = qm->nrQMatoms_get();
    snew(ptrSpecies, nAtom);
    snew(atomicNumber, nAtom);
    // read the atomic numbers of QM atoms, and set up the lists
    nSpecies = 0;
    for (int i=0; i<nAtom; i++)
    {
        atomicNumber[i] = qm->atomicnumberQM_get(i);
        // is this a new chemical element?
        bool newElement = true;
        for (int j=0; j<i; j++)
            if (atomicNumber[i] == atomicNumber[j])
                newElement = false;
        // if it is a new element, introduce it in the list
        if (newElement)
        {
            atomicNumberBySpecies[nSpecies] = atomicNumber[i];
            nSpecies++;
            ptrSpecies[i] = nSpecies; // this numbering will start at 1 (and not 0)
        }
        else
        {
            for (int k=0; k<nSpecies; k++)
                if (atomicNumber[i] == atomicNumberBySpecies[k])
                    ptrSpecies[i] = k+1; // because numbering starts at 1
        }
    }
    // assemble the list of chemical species
    for (int k=0; k<nSpecies; k++)
        strcpy(ptrElement[k], chemSymbol[atomicNumberBySpecies[k]]);
    // check what is being passed to DFTB+
    printf("This is being passed from Gromacs to DFTB+:\n");
    printf("No. of QM atoms = %d\n", nAtom);
    printf("No. of chem. species = %d:", nSpecies);
    for (int k=0; k<nSpecies; k++)
        printf(" %d=%s", k+1, ptrElement[k]);
    printf("\n");
    printf("Species by atom:\n");
    for (int i=0; i<nAtom; i++)
        printf("Atom %d is species %d\n", i+1, ptrSpecies[i]);

#if GMX_QMMM_DFTBPLUS_ATOM_LIST
    // finally, call the DFTB+ routine
    dftbp_get_atom_list(&atomList, &nAtom, &nSpecies, (char *) ptrElement, ptrSpecies);
    printf("DFTB+ has obtained the list of QM atoms!\n");
#endif

    sfree(atomicNumber);
    sfree(ptrSpecies);

    /* Set up the calculator by processing the input tree.
     * Without the atom list (the C API of the DFTB+ releases), DFTB+ takes the atoms
     *   from the Geometry block of dftb_in.hsd, which must therefore list the QM atoms
     *   in the order of the QM group, with the species printed above.
     */
#if GMX_QMMM_DFTBPLUS_ATOM_LIST
    dftbp_process_input(&calculator, &input, &atomList);
#else
    dftbp_process_input(&calculator, &input);
#endif
    printf("DFTB+ input has been processed!\n");

    const int nAtomDftb = dftbp_get_nr_atoms(&calculator);
    if (nAtomDftb != nAtom)
    {
        gmx_fatal(FARGS,
                  "DFTB+ has %d atoms in the Geometry block of dftb_in.hsd, but the QM group of "
                  "the run has %d atoms. The Geometry block must list the QM atoms in the order "
                  "of the QM group.",
                  nAtomDftb, nAtom);
    }

    qm->dpcalc = &calculator;

    /* Register the callback functions which calculate
     * the external potential and its gradient
     */
    dftbp_register_ext_pot_generator(qm->dpcalc,
                                   //cont,
                                     qm->dftbContext,
                                     calcqmextpot,
                                     calcqmextpotgrad);

    /* Alternatively - modify the geometry/elements information here,
     * based on what we have in the Gromacs topology!
     */

    return;
} /* init_dftbplus */

/* The stride in steps of a DFTB output file, from its environment variable: a positive
 * integer, or 0 (no output) for anything else, which is reported. */
static int dftbOutputStride(const char* name, const char* env)
{
    const int stride = atoi(env);
    if (stride <= 0)
    {
        printf("NOTE: %s=%s is not a positive number of steps; no output is written.\n", name, env);
        return 0;
    }
    return stride;
}

real call_dftbplus(QMMM_rec*         qr,
                // const t_commrec*  cr,
                   QMMM_QMrec*       qm,
                   const QMMM_MMrec& mm,
                   rvec              f[],
                   rvec              fshift[],
		           t_nrnb*           nrnb,
                   gmx_wallcycle*    wcycle)
{
    static int step = 0;
    static FILE *f_q = nullptr;
    static FILE *f_sh = nullptr;
    static FILE *f_p = nullptr;
    static FILE *f_x_qm = nullptr;
    static FILE *f_x_mm = nullptr;
    static FILE *f_x_mm_full = nullptr;
    static FILE *f_p_split = nullptr;
    static FILE *f_grad = nullptr;
    static FILE *f_grad_full = nullptr;
    static FILE *f_energy_corr = nullptr;
    static int output_freq_p_split;
    static int output_freq_grad;
    static int output_freq_grad_full;
    static int output_freq_energy_corr;
    static int output_freq_q;
    static int output_freq_sh;
    static int output_freq_p;
    static int output_freq_x_qm;
    static int output_freq_x_mm;
    static int output_freq_x_mm_full;

    double QMener;
 // bool lPme = (qm->qmmm_variant == eqmmmPME);
    (void) fshift;

    qm->dftbContext->nrnb = nrnb;
    int n = qm->nrQMatoms_get();

    double *x, *grad, *pot, *potgrad, *q, *atomicShifts; // real instead of rvec, to help pass data to fortran
    real *pot_sr = nullptr, *pot_lr = nullptr;
    rvec *QMgrad = nullptr, *MMgrad = nullptr, *MMgrad_full = nullptr;

    snew(x, 3*n);
    snew(grad, 3*n);
    for (int i=0; i<3*n; i++)
        grad[i] = 0.;
    snew(pot, n);
    snew(potgrad, 3*n); // dummy parameter; not used at this moment
    for (int i=0; i<3*n; i++)
        potgrad[i] = 0.;
    snew(q, n);
    for (int i=0; i<n; i++)
        q[i] = 0.;
    snew(atomicShifts, n);
    for (int i=0; i<n; i++)
        atomicShifts[i] = 0.;
    snew(pot_sr, n);
    snew(pot_lr, n);

    snew(QMgrad, qm->nrQMatoms_get());

    if (step == 0) {
        char *env;
        const bool withMM = (qm->qmmm_variant_get() != eqmmmVACUO);
        const bool withPME = (qm->qmmm_variant_get() == eqmmmPME);

        if ((env = getenv("GMX_DFTB_CHARGES")) != nullptr
            && (output_freq_q = dftbOutputStride("GMX_DFTB_CHARGES", env)) > 0)
        {
            f_q = fopen("qm_dftb_charges.xvg", "a");
            printf("The QM charges will be saved in file qm_dftb_charges.xvg every %d steps.\n", output_freq_q);
        }

        if ((env = getenv("GMX_DFTB_ATOMIC_SHIFTS")) != nullptr)
        {
            if (GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS)
            {
                if ((output_freq_sh = dftbOutputStride("GMX_DFTB_ATOMIC_SHIFTS", env)) > 0)
                {
                    f_sh = fopen("qm_dftb_atomic_shifts.xvg", "a");
                    printf("The QM atomic shifts will be saved in file qm_dftb_atomic_shifts.xvg every %d steps.\n", output_freq_sh);
                }
            }
            else
            {
                printf("NOTE: GMX_DFTB_ATOMIC_SHIFTS is ignored: the linked DFTB+ does not provide "
                       "dftbp_get_atomic_shifts() in its C API.\n");
            }
        }

        if (withMM && (env = getenv("GMX_DFTB_ESP")) != nullptr
            && (output_freq_p = dftbOutputStride("GMX_DFTB_ESP", env)) > 0)
        {
            f_p = fopen("qm_dftb_esp.xvg", "a");
            printf("The MM potential induced on QM atoms will be saved in file qm_dftb_esp.xvg every %d steps.\n", output_freq_p);
        }

        // The same potential, but with the two contributions kept apart:
        //   the MM atoms (with the boundary charge scheme) and the periodic images of the QM charges.
        if (withMM && (env = getenv("GMX_DFTB_ESP_SPLIT")) != nullptr
            && (output_freq_p_split = dftbOutputStride("GMX_DFTB_ESP_SPLIT", env)) > 0)
        {
            f_p_split = fopen("qm_dftb_esp_split.xvg", "a");
            printf("The MM and the QM-image contributions to the potential on QM atoms will be saved separately in file qm_dftb_esp_split.xvg every %d steps.\n", output_freq_p_split);
        }

        if ((env = getenv("GMX_DFTB_QM_COORD")) != nullptr
            && (output_freq_x_qm = dftbOutputStride("GMX_DFTB_QM_COORD", env)) > 0)
        {
            f_x_qm = fopen("qm_dftb_qm.qxyz", "a");
            printf("The QM coordinates (XYZQ) will be saved in file qm_dftb_qm.qxyz every %d steps.\n", output_freq_x_qm);
        }

        if (withMM && (env = getenv("GMX_DFTB_MM_COORD")) != nullptr
            && (output_freq_x_mm = dftbOutputStride("GMX_DFTB_MM_COORD", env)) > 0)
        {
            f_x_mm = fopen("qm_dftb_mm.qxyz", "a");
            printf("The MM coordinates (XYZQ) will be saved in file qm_dftb_mm.qxyz every %d steps.\n", output_freq_x_mm);
        }

        if (withMM && (env = getenv("GMX_DFTB_MM_COORD_FULL")) != nullptr
            && (output_freq_x_mm_full = dftbOutputStride("GMX_DFTB_MM_COORD_FULL", env)) > 0)
        {
            f_x_mm_full = fopen("qm_dftb_mm_full.qxyz", "a");
            printf("The full MM coordinates (XYZQ) will be saved in file qm_dftb_mm_full.qxyz every %d steps.\n", output_freq_x_mm_full);
        }

        // Gradients on the QM atoms and on the MM atoms of the short-range list.
        if (withMM && (env = getenv("GMX_DFTB_QMMM_GRAD")) != nullptr
            && (output_freq_grad = dftbOutputStride("GMX_DFTB_QMMM_GRAD", env)) > 0)
        {
            f_grad = fopen("qm_dftb_grad.xvg", "a");
            printf("The gradients on the QM atoms and on the short-range MM atoms will be saved in file qm_dftb_grad.xvg every %d steps.\n", output_freq_grad);
        }

        // Gradients on all of the MM atoms -- only meaningful with PME,
        //   where the long-range contribution is evaluated for the entire MM subsystem.
        if (withPME && (env = getenv("GMX_DFTB_QMMM_GRAD_FULL")) != nullptr
            && (output_freq_grad_full = dftbOutputStride("GMX_DFTB_QMMM_GRAD_FULL", env)) > 0)
        {
            f_grad_full = fopen("qm_dftb_grad_full.xvg", "a");
            printf("The gradients on all of the MM atoms will be saved in file qm_dftb_grad_full.xvg every %d steps.\n", output_freq_grad_full);
        }

        // The correction of the energy returned by DFTB+ (the periodic QM images, PME).
        if ((env = getenv("GMX_DFTB_ENERGY_CORR")) != nullptr
            && (output_freq_energy_corr = dftbOutputStride("GMX_DFTB_ENERGY_CORR", env)) > 0)
        {
            f_energy_corr = fopen("qm_dftb_energy_corr.xvg", "a");
            printf("The correction of the DFTB+ energy and the corrected QM energy will be saved in file qm_dftb_energy_corr.xvg every %d steps.\n", output_freq_energy_corr);
        }
    }

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<DIM; j++)
        {
            x[3*i+j] = qm->xQM_get(i,j) / gmx::c_bohr2Nm; // to bohr units for DFTB+
        }
    }

    /* calculate the QM/MM electrostatics beforehand with Gromacs!
     * with cut-off treatment, this will be the entire QM/MM
     * with PME, this will only include MM atoms,
     *    and the contribution from QM atoms will be added in every SCC iteration
     */
    auto tStart = QmmmTiming::Clock::now();
    qr->calculate_SR_QM_MM(qm->qmmm_variant_get(), pot_sr);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
     // gmx_pme_init_qmmm(&(qr->pme->pmedata), true, fr->pmedata);
        qr->calculate_LR_QM_MM(nrnb, wcycle, pot_lr); // cr ... *qr->pmedata, pot_lr);
        for (int i=0; i<n; i++)
        {
            pot[i] = (double) - (pot_sr[i] + pot_lr[i]);
        }
    }
    else
    {
        for (int i=0; i<n; i++)
        {
            pot_lr[i] = 0.;
            pot[i] = (double) - pot_sr[i];
        }
    }
    // save the potential in the QMMM_QMrec structure
    for (int j=0; j<n; j++)
    {
        qm->pot_qmmm_set(j, (double) - pot[j] * HARTREE_TO_EV); // in volt units
    }
 // // DEBUG
 // for (int i=0; i<n; i++)
 //     printf("pot_sr[%d] = %9.5f pot_lr[%d] = %9.5f\n", i+1, pot_sr[i], i+1, pot_lr[i]);

    /* Set up the data structures needed for the PME calculation
	 *   of QM--imageQM electrostatics.
	 * During the iterative SCC calculation, this routine will be called
	 *   directly from DFTB+, and will use those data structures.
	 */
 //  calcQMextPotPME(nullptr, nullptr, true, lPme, fr, cr, wcycle, fr->rcoulomb, fr->ewaldcoeff_q);
    /* This was already done in the initialization procedure! */

    for (int i=0; i<n; i++)
    {
        q[i] = 0.;
    }

    qmmmTiming.mmPotential += QmmmTiming::since(tStart);
    tStart = QmmmTiming::Clock::now();
    /* DFTB+ calculation itself, on its OpenMP threads (dftbNumThreads()).
     * The callback for the periodic QM images (PME) runs PME on the PME threads of mdrun. */
    wallcycle_start(wcycle, WallCycleCounter::QM);
    DftbThreadScope threadScope(dftbNumThreads());
    dftbp_set_coords(qm->dpcalc, x); // unit OK
    dftbp_set_external_potential(qm->dpcalc, pot, potgrad); // unit and sign OK
    dftbp_get_energy(qm->dpcalc, &QMener); // unit OK
    dftbp_get_gross_charges(qm->dpcalc, q);
 // for (int i=0; i<n; i++)
 //     printf("%d %6.3f\n", i+1, q[i]);
#if GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS
    dftbp_get_atomic_shifts(qm->dpcalc, atomicShifts);
#endif
    dftbp_get_gradients(qm->dpcalc, grad);
    wallcycle_stop(wcycle, WallCycleCounter::QM);
    qmmmTiming.dftb += QmmmTiming::since(tStart);

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) grad[3*i+j]; // negative of force -- sign OK
        }
    }

 // /* Print the QM pure gradient */
 // for (int i=0; i<n; i++)
 // {
 //     printf("GRAD QM %d: %8.2f %8.2f %8.2f\n", i+1,
 //         QMgrad[i][XX] * HARTREE_BOHR2MD, QMgrad[i][YY] * HARTREE_BOHR2MD, QMgrad[i][ZZ] * HARTREE_BOHR2MD);
 // }

    /* Save the QM charges */
    for (int i=0; i<n; i++)
    {
        qm->QMcharges_set(i, (real) q[i]); // sign OK
     // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qm->QMcharges[i]);
        qm->QMatomicShifts_set(i, (real) atomicShifts[i] * HARTREE_TO_EV); // in volt units
    }

    /* Calculate the QM/MM forces
     *   (these are not covered by DTFB+,
     *    because it does not know the position and charges of individual atoms,
     *    and instead, it only obtains the external potentials induced at QM atoms.
     */
    snew(MMgrad, mm.nrMMatoms);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(MMgrad_full, mm.nrMMatoms_full);
    }

    /* The interaction of the QM charges with their own periodic images (PME).
     * The callback hands DFTB+ the image potential V_img in every SCC iteration, and DFTB+
     *   counts q.V_img in its energy in full. The Ewald energy of a charge distribution with
     *   its own images is 1/2 q.V_img -- the other half would count every image pair twice --
     *   and the image forces of gradient_QM_MM() (and the virial) are those of the halved
     *   term. The energy is brought in line with them here.
     */
    const double QMenerDftb = QMener;
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        double eImage = 0.;
        for (int i=0; i<n; i++)
        {
            eImage += q[i] * qm->pot_qmqm_get(i) / HARTREE_TO_EV;
        }
        QMener -= 0.5 * eImage;
    }
    if (f_energy_corr && step % output_freq_energy_corr == 0)
    {
        if (step == 0)
        {
            fprintf(f_energy_corr, "# step, energy of DFTB+, correction for the periodic QM images, corrected QM energy (kJ/mol)\n");
        }
        const double toKJ = gmx::c_hartree2Kj * gmx::c_avogadro;
        fprintf(f_energy_corr, "%10d %20.10f %20.10f %20.10f\n", step, QMenerDftb * toKJ,
                (QMener - QMenerDftb) * toKJ, QMener * toKJ);
        fflush(f_energy_corr);
    }

    rvec *partgrad;
    snew(partgrad, qm->nrQMatoms_get());
    tStart = QmmmTiming::Clock::now();
    qr->gradient_QM_MM(nrnb, wcycle, // cr ... (qm->qmmm_variant_get() == eqmmmPME ? *qr->pmedata : nullptr),
                   qm->qmmm_variant_get(), partgrad, MMgrad, MMgrad_full);
    qmmmTiming.gradient += QmmmTiming::since(tStart);
    qmmmTiming.numSteps++;
    {
        static const int timingStride = getenv("GMX_QMMM_TIMING") ? atoi(getenv("GMX_QMMM_TIMING")) : 0;
        if (timingStride > 0 && qmmmTiming.numSteps % timingStride == 0)
        {
            const double ns = qmmmTiming.numSteps;
            printf("QM/MM timing (ms/step, average of %d steps): MM potential %.2f, DFTB+ %.2f "
                   "(of which QM-image PME %.2f in %.1f calls), gradients %.2f, %d MM atoms in the short-range list\n",
                   qmmmTiming.numSteps, qmmmTiming.mmPotential / ns, qmmmTiming.dftb / ns,
                   qmmmTiming.imageCallback / ns, qmmmTiming.numCallbacks / ns, qmmmTiming.gradient / ns,
                   mm.nrMMatoms);
            fflush(stdout);
        }
    }

    /* Optionally, write out the gradients while they are still separated.
     * At this point, and in atomic units (hartree/bohr):
     *   QMgrad[]   is the gradient obtained from DFTB+, i.e. the QM subsystem itself;
     *   partgrad[] is the electrostatic gradient due to the environment -- with PME, this
     *     includes the periodic images of the QM charges;
     *   MMgrad[]   is the electrostatic gradient on the MM atoms of the short-range list,
     *     the forces of the fictitious charges of the boundary scheme included.
     * The sum of the first two is the total gradient on the QM atom.
     */
    if (f_grad && step % output_freq_grad == 0)
    {
        fprintf(f_grad, "\nQM gradients: DFTB+, electrostatic, total (hartree/bohr) step %d\n", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_grad, "QM %5d %12.7f%12.7f%12.7f %12.7f%12.7f%12.7f %12.7f%12.7f%12.7f\n", i+1,
                QMgrad[i][XX], QMgrad[i][YY], QMgrad[i][ZZ],
                partgrad[i][XX], partgrad[i][YY], partgrad[i][ZZ],
                QMgrad[i][XX] + partgrad[i][XX],
                QMgrad[i][YY] + partgrad[i][YY],
                QMgrad[i][ZZ] + partgrad[i][ZZ]);
        }
        fprintf(f_grad, "MM gradients on the short-range list (hartree/bohr) step %d\n", step);
        for (int i=0; i<mm.nrMMatoms; i++)
        {
            fprintf(f_grad, "MM %5d %8d %12.7f%12.7f%12.7f\n", i+1, mm.indexMM[i] + 1,
                MMgrad[i][XX], MMgrad[i][YY], MMgrad[i][ZZ]);
        }
        fflush(f_grad);
    }

    if (f_grad_full && step % output_freq_grad_full == 0)
    {
        fprintf(f_grad_full, "\nMM gradients on all MM atoms, reciprocal space (hartree/bohr) step %d\n", step);
        for (int i=0; i<mm.nrMMatoms_full; i++)
        {
            fprintf(f_grad_full, "%8d %12.7f%12.7f%12.7f\n", mm.indexMM_full[i] + 1,
                MMgrad_full[i][XX], MMgrad_full[i][YY], MMgrad_full[i][ZZ]);
        }
        fflush(f_grad_full);
    }
    for (int i=0; i<n; i++)
    {
        rvec_inc(QMgrad[i], partgrad[i]); // sign OK
     // printf("GRAD QM FULL %d: %8.2f %8.2f %8.2f\n", i+1,
     //     QMgrad[i][XX] * HARTREE_BOHR2MD, QMgrad[i][YY] * HARTREE_BOHR2MD, QMgrad[i][ZZ] * HARTREE_BOHR2MD);
    }
    sfree(partgrad);

    /* Put the QMMM forces in the force array and to the fshift.
     * Convert to MD units.
     */
    for (int i = 0; i < qm->nrQMatoms_get(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i][j]      = gmx::c_hartreeBohr2Md*QMgrad[i][j];
         // fshift[i][j] = gmx::c_hartreeBohr2Md*QMgrad[i][j];
        }
    }
    for (int i = 0; i < mm.nrMMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms_get()][j]      = gmx::c_hartreeBohr2Md*MMgrad[i][j];
         // fshift[i+qm->nrQMatoms_get()][j] = gmx::c_hartreeBohr2Md*MMgrad[i][j];
        }
    }
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        for (int i = 0; i < mm.nrMMatoms_full; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                f[i+qm->nrQMatoms_get()+mm.nrMMatoms][j]      = gmx::c_hartreeBohr2Md*MMgrad_full[i][j];
             // fshift[i+qm->nrQMatoms_get()+mm.nrMMatoms][j] = gmx::c_hartreeBohr2Md*MMgrad_full[i][j];
            }
        }
    }

    if (f_q && step % output_freq_q == 0)
    {
        fprintf(f_q, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_q, " %8.5f", qm->QMcharges_get(i));
        }
        fprintf(f_q, "\n");
    }

    if (f_sh && step % output_freq_sh == 0)
    {
        fprintf(f_sh, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_sh, " %8.5f", qm->QMatomicShifts_get(i));
        }
        fprintf(f_sh, "\n");
    }

    if (f_p && step % output_freq_p == 0)
    {
        fprintf(f_p, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_p, " %8.5f", qm->pot_qmmm_get(i) + qm->pot_qmqm_get(i));
         // if (qm.qmmm_variant == eqmmmPME)
         // {
         //     fprintf(f_p, " %8.5f %8.5f %8.5f", qm.pot_qmmm[i], qm.pot_qmqm[i], qm.pot_qmmm[i] + qm.pot_qmqm[i]);
         // }
         // else
         // {
         //     fprintf(f_p, " %8.5f", qm.pot_qmmm[i]);
         // }
        }
        fprintf(f_p, "\n");
    }

    /* The same potential as above, with the two contributions written out separately:
     *   the one induced by the MM atoms (with the boundary charge scheme), and the one induced
     *   by the periodic images of the QM charges (identically zero unless PME is used).
     */
    if (f_p_split && step % output_freq_p_split == 0)
    {
        fprintf(f_p_split, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_p_split, " %8.5f %8.5f %8.5f",
                qm->pot_qmmm_get(i), qm->pot_qmqm_get(i),
                qm->pot_qmmm_get(i) + qm->pot_qmqm_get(i));
        }
        fprintf(f_p_split, "\n");
    }

    if (f_x_qm && step % output_freq_x_qm == 0)
    {
        const char periodic_system[37][3]={"XX",
            "h",                               "he",
            "li","be","b", "c", "n", "o", "f", "ne",
            "na","mg","al","si","p", "s", "cl","ar",
            "k", "ca","sc","ti","v", "cr","mn","fe","co",
            "ni","cu","zn","ga","ge","as","se","br","kr"};
        fprintf(f_x_qm, "\nQM coordinates and charges step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_x_qm, "%-2s %10.5f%10.5f%10.5f %10.7f\n",
                periodic_system[qm->atomicnumberQM_get(i)],
                qm->xQM_get(i, 0) * 10., qm->xQM_get(i, 1) * 10., qm->xQM_get(i, 2) * 10.,
                qm->QMcharges_get(i));
        }
    }

    if (f_x_mm && step % output_freq_x_mm == 0)
    {
        /* The charges of the QM--MM electrostatics: with a boundary charge scheme a removed MM1
         * atom has zero charge here, the AMBER shares are included, and the fictitious charges
         * follow the MM atoms, so that the potential on the QM atoms can be recomputed from
         * this file alone (with the cut-off variants; with PME the reciprocal space adds to it).
         */
        fprintf(f_x_mm, "\nMM coordinates and charges step %d\n", step);
        for (int i=0; i<mm.nrMMatoms; i++) {
            fprintf(f_x_mm, "%10.7f%10.5f%10.5f%10.5f\n",
                qr->qmmmChargesSR[i] * qr->qmmmScaleSR[i],
                mm.xMM[i][0] * 10., mm.xMM[i][1] * 10., mm.xMM[i][2] * 10.);
        }
        if (!qr->potPoints.empty())
        {
            fprintf(f_x_mm, "fictitious charges of the boundary scheme step %d\n", step);
            for (const QMMM_rec::PotPoint& pt : qr->potPoints)
            {
                rvec xp;
                qr->boundary_point_position(pt, xp);
                fprintf(f_x_mm, "%10.7f%10.5f%10.5f%10.5f\n", pt.q * mm.scalefactor,
                        xp[0] * 10., xp[1] * 10., xp[2] * 10.);
            }
        }
    }

    if (f_x_mm_full && step % output_freq_x_mm_full == 0)
    {
        fprintf(f_x_mm_full, "\nbox step %d : %10.5f%10.5f%10.5f\n", step, qm->box_xx_get(), qm->box_yy_get(), qm->box_zz_get());
        fprintf(f_x_mm_full, "\nfull MM coordinates and charges step %d\n", step);
        for (int i=0; i<mm.nrMMatoms_full; i++) {
            fprintf(f_x_mm_full, "%10.7f%10.5f%10.5f%10.5f\n",
                qr->qmmmChargesFull[i], // the charges on the PME grid (incl. the AMBER shares)
                mm.xMM_full[i][0] * 10., mm.xMM_full[i][1] * 10., mm.xMM_full[i][2] * 10.);
        }
    }

    sfree(QMgrad);
    sfree(MMgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(MMgrad_full);
    }

    sfree(x);
    sfree(grad);
    sfree(pot);
    sfree(potgrad);
    sfree(q);
    sfree(atomicShifts);
    sfree(pot_sr);
    sfree(pot_lr);

    step++;

    return (real) QMener * gmx::c_hartree2Kj * gmx::c_avogadro;
} /* call_dftbplus */

/* end of dftbplus sub routines */

#endif

#pragma GCC diagnostic pop

