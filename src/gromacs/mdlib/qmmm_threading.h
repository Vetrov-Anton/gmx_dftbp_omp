/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */
/*! \internal \file
 * \brief OpenMP helpers of the QM/MM interface (one MPI rank, several OpenMP threads).
 *
 * The QM/MM interface runs on a single rank. Its loops over the MM atoms use the OpenMP
 * threads of mdrun (-ntomp), the same pool that the MM force kernels, PME and DFTB+ use.
 * The MM atoms are split into contiguous chunks, one per thread, and the per-QM-atom sums
 * are collected in buffers of their own and added up in the order of the threads, so a run
 * is reproducible for a given number of threads, and a run on one thread adds the terms in
 * the order of the serial code.
 */
#ifndef GMX_MDLIB_QMMM_THREADING_H
#define GMX_MDLIB_QMMM_THREADING_H

#include <algorithm>

#include "gromacs/mdlib/gmx_omp_nthreads.h"

namespace qmmm_omp
{

//! The OpenMP threads of mdrun for the non-PME work.
inline int maxThreads()
{
    return std::max(1, gmx_omp_nthreads_get(ModuleMultiThread::Default));
}

/*! \brief The number of threads for a loop over \p numItems items
 *
 * At least \p minItemsPerThread items per thread, so short loops (or a small system)
 * do not pay for the start of a parallel region.
 */
inline int numThreads(int numItems, int minItemsPerThread = 64)
{
    /* All of the threads or one: a team of another size makes libgomp let the surplus threads
     * of its pool exit, and the threads created anew later do not have the pinning of mdrun
     * (see qmmmShareThreadAffinity()). */
    return numItems >= 2 * std::max(1, minItemsPerThread) ? maxThreads() : 1;
}

//! The first item of chunk \p thread of \p numThreads over \p numItems items
inline int chunkStart(int numItems, int thread, int numThreads)
{
    return static_cast<int>((static_cast<long long>(numItems) * thread) / numThreads);
}

} // namespace qmmm_omp

/*! \brief Lets OpenMP threads created later run on the cores of mdrun
 *
 * mdrun pins every thread of its OpenMP pool to one core. DFTB+ and OpenBLAS start teams of
 * other sizes, and libgomp then lets the surplus threads of its pool exit and creates new ones
 * later, which inherit the affinity of the main thread, i.e. a single core. To keep them on
 * the cores of mdrun, the main thread gets the union of the affinity masks of the team.
 * Call after the threads have been pinned (gmx_set_thread_affinity()). Linux only.
 *
 * \returns the number of cores in the union, 0 if nothing was changed
 */
int qmmmShareThreadAffinity(int numThreads);

#endif
