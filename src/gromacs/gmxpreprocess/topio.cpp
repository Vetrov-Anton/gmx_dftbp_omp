/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "topio.h"

#include <cassert>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>

#include <sys/types.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gmxcpp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/readir.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/gmxpreprocess/topshake.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/gmxpreprocess/vsite_parm.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct t_nbparam;

#define OPENDIR '['  /* starting sign for directive */
#define CLOSEDIR ']' /* ending sign for directive   */

static void gen_pairs(const InteractionsOfType& nbs, InteractionsOfType* pairs, real fudge, CombinationRule comb)
{
    real scaling;
    int  ntp = nbs.size();
    int  nnn = static_cast<int>(std::sqrt(static_cast<double>(ntp)));
    GMX_ASSERT(nnn * nnn == ntp,
               "Number of pairs of generated non-bonded parameters should be a perfect square");
    int nrfp  = NRFP(InteractionFunction::LennardJonesShortRange);
    int nrfpA = interaction_function[InteractionFunction::LennardJones14].nrfpA;
    int nrfpB = interaction_function[InteractionFunction::LennardJones14].nrfpB;

    if ((nrfp != nrfpA) || (nrfpA != nrfpB))
    {
        gmx_incons("Number of force parameters in gen_pairs wrong");
    }

    fprintf(stderr, "Generating 1-4 interactions: fudge = %g\n", fudge);
    pairs->interactionTypes.clear();
    int                             i = 0;
    std::array<int, 2>              atomNumbers;
    std::array<real, MAXFORCEPARAM> forceParam = { NOTSET };
    for (const auto& type : nbs.interactionTypes)
    {
        /* Copy type.atoms */
        atomNumbers = { i / nnn, i % nnn };
        /* Copy normal and FEP parameters and multiply by fudge factor */
        gmx::ArrayRef<const real> existingParam = type.forceParam();
        GMX_RELEASE_ASSERT(2 * nrfp <= MAXFORCEPARAM,
                           "Can't have more parameters than half of maximum parameter number");
        for (int j = 0; j < nrfp; j++)
        {
            /* If we are using sigma/epsilon values, only the epsilon values
             * should be scaled, but not sigma.
             * The sigma values have even indices 0,2, etc.
             */
            if ((comb == CombinationRule::Arithmetic || comb == CombinationRule::GeomSigEps)
                && (j % 2 == 0))
            {
                scaling = 1.0;
            }
            else
            {
                scaling = fudge;
            }

            forceParam[j]        = scaling * existingParam[j];
            forceParam[nrfp + j] = scaling * existingParam[j];
        }
        pairs->interactionTypes.emplace_back(atomNumbers, forceParam);
        i++;
    }
}

double check_mol(const gmx_mtop_t* mtop, WarningHandler* wi)
{
    char   buf[256];
    int    i, ri;
    double q;
    real   m, mB;

    /* Check mass and charge */
    q = 0.0;

    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const t_atoms* atoms = &mtop->moltype[molb.type].atoms;
        for (i = 0; (i < atoms->nr); i++)
        {
            q += molb.nmol * atoms->atom[i].q;
            m               = atoms->atom[i].m;
            mB              = atoms->atom[i].mB;
            ParticleType pt = atoms->atom[i].ptype;
            /* If the particle is an atom or a nucleus it must have a mass,
             * else, if it is a shell, a vsite or a bondshell it can have mass zero
             */
            if (((m <= 0.0) || (mB <= 0.0)) && ((pt == ParticleType::Atom) || (pt == ParticleType::Nucleus)))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf,
                        "atom %s (Res %s-%d) has mass %g (state A) / %g (state B)\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m,
                        mB);
                wi->addError(buf);
            }
            else if (((m != 0) || (mB != 0)) && (pt == ParticleType::VSite))
            {
                ri = atoms->atom[i].resind;
                sprintf(buf,
                        "virtual site %s (Res %s-%d) has non-zero mass %g (state A) / %g (state "
                        "B)\n"
                        "     Check your topology.\n",
                        *(atoms->atomname[i]),
                        *(atoms->resinfo[ri].name),
                        atoms->resinfo[ri].nr,
                        m,
                        mB);
                wi->addError(buf);
                /* The following statements make LINCS break! */
                /* atoms->atom[i].m=0; */
            }
        }
    }
    return q;
}

/*! \brief Describe molecule and involved atoms for a dihedral interaction of a given type
 *
 * Searches in the dihedrals of type 3 (Ryckaert-Bellemans or Fourier) for an interaction
 * type matching the interactionType input parameters, returning for the first match
 * a string with the name of the  molecule and the 4 involved atoms.
 *
 * Precondition: the interaction should exist in the topology
 *
 */
static std::string describeAtomsForRBDihedralOfGivenType(const gmx_mtop_t& mtop, int interactionType)
{
    for (const auto& molt : mtop.moltype)
    {
        const int* ia = molt.ilist[InteractionFunction::RyckaertBellemansDihedrals].iatoms.data();
        for (int i = 0; (i < molt.ilist[InteractionFunction::RyckaertBellemansDihedrals].size());)
        {
            const int                 type  = ia[0];
            const InteractionFunction ftype = mtop.ffparams.functype[type];
            const int                 nra   = interaction_function[ftype].nratoms;
            if (type == interactionType)
            {
                return gmx::formatString(
                        "First such dihedral in molecule %s, involving atoms %d %d %d %d",
                        *(molt.name),
                        ia[1],
                        ia[2],
                        ia[3],
                        ia[4]);
            }
            ia += nra + 1;
            i += nra + 1;
        }
    }
    gmx_fatal(FARGS, "Precondition violation: could not find RB interaction of given type %d", interactionType);
}

void checkRBDihedralSum(const gmx_mtop_t& mtop, const t_inputrec& ir, WarningHandler* wi)
{
    /*
     * The sum of the RB dihedral coefficient being zero is relevant when:
     *     - Free energy computation is being performed (dHdl)
     *     - Comparing energies between force field ports and/or other MD codes
     *  because this affect the value of the potential energy.
     *
     *  We can clearly detect problems in the first situation: free energy is enabled, stateA
     *  and stateB have a different sum (=> potential energy "offset" at 0 degree), so we emit a
     *  warning. The second case is more subtle, since formally the potential is up to a constant,
     *  which does not affect the dynamics. We therefore only emit a note.
     */

    // Mistakes here are typically due to using the wrong formula to port dihedrals, not numerical
    // issues, so we use a relatively large tolerance
    const real         absoluteTolerance        = 0.01;
    int                numSumCoefficientNotZero = 0;
    std::optional<int> indexOfFirstSumCoefficientNotZero;
    int                numSumCoefficientDifferentInStateAStateB = 0;
    std::optional<int> indexOfSumCoefficientDifferentInStateAStateB;

    for (int j = 0; j < gmx::ssize(mtop.ffparams.functype); ++j)
    {
        InteractionFunction ftype = mtop.ffparams.functype[j];
        if (ftype == InteractionFunction::RyckaertBellemansDihedrals)
        {
            const t_iparams& params = mtop.ffparams.iparams[j];
            const real       sum_a =
                    std::accumulate(std::begin(params.rbdihs.rbcA), std::end(params.rbdihs.rbcA), 0.0);
            const real sum_b =
                    std::accumulate(std::begin(params.rbdihs.rbcB), std::end(params.rbdihs.rbcB), 0.0);

            if (std::abs(sum_a - sum_b) > absoluteTolerance)
            {
                numSumCoefficientDifferentInStateAStateB++;
                if (!indexOfSumCoefficientDifferentInStateAStateB.has_value())
                {
                    indexOfSumCoefficientDifferentInStateAStateB = j;
                }
            }

            if (std::abs(sum_a) > absoluteTolerance || std::abs(sum_b) > absoluteTolerance)
            {
                numSumCoefficientNotZero++;
                if (!indexOfFirstSumCoefficientNotZero.has_value())
                {
                    indexOfFirstSumCoefficientNotZero = j;
                }
            }
        }
    }

    // At this stage of grompp, we only have the interactions that are going to be used in the
    // simulation - this eliminates warning unrelated to the user's system, but also unfortunately
    // means that we cannot easily identify the file/line where the offending parameters were
    // originally defined. We can however go through the molecule types and print one instance of
    // the offending dihedral.

    auto generateMessage = [](int numDihedrals, const std::string& note, const std::string& involvedAtoms)
    {
        return gmx::formatString(
                "%d dihedrals with function type 3 (Ryckaert-Bellemans or Fourier) have "
                "coefficients %s"
                "\n%s",
                numDihedrals,
                note.c_str(),
                involvedAtoms.c_str());
    };

    if (numSumCoefficientNotZero > 0)
    {
        std::string involvedAtoms =
                describeAtomsForRBDihedralOfGivenType(mtop, indexOfFirstSumCoefficientNotZero.value());
        std::string note =
                "that do not sum to zero. This does not affect the simulation and can "
                "be ignored, unless you are comparing potential energy values with other force "
                "field ports and/or MD software.";
        std::string message = generateMessage(numSumCoefficientNotZero, note, involvedAtoms);
        wi->addNote(message);
    }

    if (numSumCoefficientDifferentInStateAStateB > 0 && ir.efep != FreeEnergyPerturbationType::No)
    {
        std::string involvedAtoms = describeAtomsForRBDihedralOfGivenType(
                mtop, indexOfSumCoefficientDifferentInStateAStateB.value());
        std::string note =
                "whose sums do not match in state A and B. This could introduce an "
                "undesired offset in dHdl values.";
        std::string message =
                generateMessage(numSumCoefficientDifferentInStateAStateB, note, involvedAtoms);
        wi->addWarning(message);
    }
}

/*! \brief Returns the rounded charge of a molecule, when close to integer, otherwise returns the original charge.
 *
 * The results of this routine are only used for checking and for
 * printing warning messages. Thus we can assume that charges of molecules
 * should be integer. If the user wanted non-integer molecular charge,
 * an undesired warning is printed and the user should use grompp -maxwarn 1.
 *
 * \param qMol     The total, unrounded, charge of the molecule
 * \param sumAbsQ  The sum of absolute values of the charges, used for determining the tolerance for the rounding.
 */
static double roundedMoleculeCharge(double qMol, double sumAbsQ)
{
    /* We use a tolerance of 1e-6 for inaccuracies beyond the 6th decimal
     * of the charges for ascii float truncation in the topology files.
     * Although the summation here uses double precision, the charges
     * are read and stored in single precision when real=float. This can
     * lead to rounding errors of half the least significant bit.
     * Note that, unfortunately, we can not assume addition of random
     * rounding errors. It is not entirely unlikely that many charges
     * have a near half-bit rounding error with the same sign.
     */
    double tolAbs = 1e-6;
    double tol    = std::max(tolAbs, 0.5 * GMX_REAL_EPS * sumAbsQ);
    double qRound = std::round(qMol);
    if (std::abs(qMol - qRound) <= tol)
    {
        return qRound;
    }
    else
    {
        return qMol;
    }
}

static void sum_q(const t_atoms* atoms, int numMols, double* qTotA, double* qTotB)
{
    /* sum charge */
    double qmolA    = 0;
    double qmolB    = 0;
    double sumAbsQA = 0;
    double sumAbsQB = 0;
    for (int i = 0; i < atoms->nr; i++)
    {
        qmolA += atoms->atom[i].q;
        qmolB += atoms->atom[i].qB;
        sumAbsQA += std::abs(atoms->atom[i].q);
        sumAbsQB += std::abs(atoms->atom[i].qB);
    }

    *qTotA += numMols * roundedMoleculeCharge(qmolA, sumAbsQA);
    *qTotB += numMols * roundedMoleculeCharge(qmolB, sumAbsQB);
}

static void get_nbparm(char* nb_str, char* comb_str, VanDerWaalsPotential* nb, CombinationRule* comb, WarningHandler* wi)
{
    *nb = VanDerWaalsPotential::Count;
    for (auto i : gmx::EnumerationArray<VanDerWaalsPotential, bool>::keys())
    {
        if (gmx_strcasecmp(nb_str, enumValueToString(i)) == 0)
        {
            *nb = i;
        }
    }
    if (*nb == VanDerWaalsPotential::Count)
    {
        int integerValue = std::strtol(nb_str, nullptr, 10);
        if ((integerValue < 1) || (integerValue >= static_cast<int>(VanDerWaalsPotential::Count)))
        {
            std::string message =
                    gmx::formatString("Invalid nonbond function selector '%s' using %s",
                                      nb_str,
                                      enumValueToString(VanDerWaalsPotential::LJ));
            wi->addError(message);
            *nb = VanDerWaalsPotential::LJ;
        }
        else
        {
            *nb = static_cast<VanDerWaalsPotential>(integerValue);
        }
    }
    *comb = CombinationRule::Count;
    for (auto i : gmx::EnumerationArray<CombinationRule, bool>::keys())
    {
        if (gmx_strcasecmp(comb_str, enumValueToString(i)) == 0)
        {
            *comb = i;
        }
    }
    if (*comb == CombinationRule::Count)
    {
        int integerValue = std::strtol(comb_str, nullptr, 10);
        if ((integerValue < 1) || (integerValue >= static_cast<int>(CombinationRule::Count)))
        {
            std::string message =
                    gmx::formatString("Invalid combination rule selector '%s' using %s",
                                      comb_str,
                                      enumValueToString(CombinationRule::Geometric));
            wi->addError(message);
            *comb = CombinationRule::Geometric;
        }
        else
        {
            *comb = static_cast<CombinationRule>(integerValue);
        }
    }
}

/*! \brief Parses define and include flags.
 *
 * Returns a vector of parsed include/define flags, with an extra nullptr entry at the back
 * for consumers that expect null-terminated char** structures.
 */
static std::vector<char*> cpp_opts(const char* define, const char* include, WarningHandler* wi)
{
    int         n, len;
    const char* cppadds[2];
    const char* option[2] = { "-D", "-I" };
    const char* nopt[2]   = { "define", "include" };
    const char* ptr;
    const char* rptr;
    char*       buf;
    char        warn_buf[STRLEN];

    cppadds[0] = define;
    cppadds[1] = include;
    std::vector<char*> cppOptions;
    for (n = 0; (n < 2); n++)
    {
        if (cppadds[n])
        {
            ptr = cppadds[n];
            while (*ptr != '\0')
            {
                while ((*ptr != '\0') && std::isspace(*ptr))
                {
                    ptr++;
                }
                rptr = ptr;
                while ((*rptr != '\0') && !std::isspace(*rptr))
                {
                    rptr++;
                }
                len = (rptr - ptr);
                if (len > 2)
                {
                    snew(buf, (len + 1));
                    std::strncpy(buf, ptr, len);
                    if (std::strstr(ptr, option[n]) != ptr)
                    {
                        wi->setFileAndLineNumber("mdp file", -1);
                        sprintf(warn_buf, "Malformed %s option %s", nopt[n], buf);
                        wi->addWarning(warn_buf);
                    }
                    else
                    {
                        cppOptions.emplace_back(gmx_strdup(buf));
                    }
                    sfree(buf);
                    ptr = rptr;
                }
            }
        }
    }
    // Users of cppOptions expect a null last element.
    cppOptions.emplace_back(nullptr);
    return cppOptions;
}


static void make_atoms_sys(gmx::ArrayRef<const gmx_molblock_t>      molblock,
                           gmx::ArrayRef<const MoleculeInformation> molinfo,
                           t_atoms*                                 atoms)
{
    atoms->nr   = 0;
    atoms->atom = nullptr;

    for (const gmx_molblock_t& molb : molblock)
    {
        const t_atoms& mol_atoms = molinfo[molb.type].atoms;

        srenew(atoms->atom, atoms->nr + molb.nmol * mol_atoms.nr);

        for (int m = 0; m < molb.nmol; m++)
        {
            for (int a = 0; a < mol_atoms.nr; a++)
            {
                atoms->atom[atoms->nr++] = mol_atoms.atom[a];
            }
        }
    }
}


static char** read_topol(const char*                                 infile,
                         const std::optional<std::filesystem::path>& outfile,
                         const char*                                 define,
                         const char*                                 include,
                         t_symtab*                                   symtab,
                         PreprocessingAtomTypes*                     atypes,
                         std::vector<MoleculeInformation>*           molinfo,
                         std::unique_ptr<MoleculeInformation>*       intermolecular_interactions,
                         gmx::EnumerationArray<InteractionFunction, InteractionsOfType>& interactions,
                         CombinationRule*             combination_rule,
                         double*                      reppow,
                         t_gromppopts*                opts,
                         real*                        fudgeQQ,
                         std::vector<gmx_molblock_t>* molblock,
                         bool*                        ffParametrizedWithHBondConstraints,
                         bool                         bFEP,
                         bool                         bZero,
                         bool                         usingFullRangeElectrostatics,
                         WarningHandler*              wi,
                         const gmx::MDLogger&         logger)
{
    FILE*                out;
    int                  sl;
    char *               pline = nullptr, **title = nullptr;
    char                 line[STRLEN], errbuf[256], comb_str[256], nb_str[256];
    char                 genpairs[32];
    char *               dirstr, *dummy2;
    int                  nrcopies, nscan, ncombs, ncopy;
    double               fLJ, fQQ, fPOW;
    MoleculeInformation* mi0 = nullptr;
    DirStack*            DS;
    Directive            d, newd;
    t_nbparam **         nbparam, **pair;
    real                 fudgeLJ = -1; /* Multiplication factor to generate 1-4 from LJ */
    bool                 bReadDefaults, bReadMolType, bGenPairs, bWarn_copy_A_B;
    double               qt = 0, qBt = 0; /* total charge */
    int                  dcatt = -1, nmol_couple;
    /* File handling variables */
    int         status;
    bool        done;
    gmx_cpp_t   handle;
    char*       tmp_line = nullptr;
    char        warn_buf[STRLEN];
    const char* floating_point_arithmetic_tip =
            "Total charge should normally be an integer. See\n"
            "https://manual.gromacs.org/current/user-guide/floating-point.html\n"
            "for discussion on how close it should be to an integer.\n";
    /* We need to open the output file before opening the input file,
     * because cpp_open_file can change the current working directory.
     */
    if (outfile)
    {
        out = gmx_fio_fopen(outfile.value(), "w");
    }
    else
    {
        out = nullptr;
    }

    /* open input file */
    auto cpp_opts_return = cpp_opts(define, include, wi);
    status               = cpp_open_file(infile, &handle, cpp_opts_return.data());
    if (status != 0)
    {
        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
    }

    /* some local variables */
    DS_Init(&DS);                   /* directive stack	*/
    d       = Directive::d_invalid; /* first thing should be a directive */
    nbparam = nullptr;              /* The temporary non-bonded matrix */
    pair    = nullptr;              /* The temporary pair interaction matrix */
    std::vector<std::vector<gmx::ExclusionBlock>> exclusionBlocks;
    VanDerWaalsPotential                          nb_funct = VanDerWaalsPotential::LJ;

    *reppow = 12.0; /* Default value for repulsion power     */

    /* Init the number of CMAP torsion angles */
    interactions[InteractionFunction::DihedralEnergyCorrectionMap].numCmaps_ = 0;

    bWarn_copy_A_B = bFEP;

    PreprocessingBondAtomType bondAtomType;
    /* parse the actual file */
    bReadDefaults = FALSE;
    bGenPairs     = FALSE;
    bReadMolType  = FALSE;
    nmol_couple   = 0;

    do
    {
        status = cpp_read_line(&handle, STRLEN, line);
        done   = (status == eCPP_EOF);
        if (!done)
        {
            if (status != eCPP_OK)
            {
                gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
            }
            else if (out)
            {
                fprintf(out, "%s\n", line);
            }

            wi->setFileAndLineNumber(cpp_cur_file(&handle), cpp_cur_linenr(&handle));

            pline = gmx_strdup(line);

            /* Strip trailing '\' from pline, if it exists */
            sl = std::strlen(pline);
            if ((sl > 0) && (pline[sl - 1] == CONTINUE))
            {
                pline[sl - 1] = ' ';
            }

            /* build one long line from several fragments - necessary for CMAP */
            while (continuing(line))
            {
                status = cpp_read_line(&handle, STRLEN, line);
                wi->setFileAndLineNumber(cpp_cur_file(&handle), cpp_cur_linenr(&handle));

                /* Since we depend on the '\' being present to continue to read, we copy line
                 * to a tmp string, strip the '\' from that string, and cat it to pline
                 */
                tmp_line = gmx_strdup(line);

                sl = std::strlen(tmp_line);
                if ((sl > 0) && (tmp_line[sl - 1] == CONTINUE))
                {
                    tmp_line[sl - 1] = ' ';
                }

                done = (status == eCPP_EOF);
                if (!done)
                {
                    if (status != eCPP_OK)
                    {
                        gmx_fatal(FARGS, "%s", cpp_error(&handle, status));
                    }
                    else if (out)
                    {
                        fprintf(out, "%s\n", line);
                    }
                }

                srenew(pline, std::strlen(pline) + std::strlen(tmp_line) + 1);
                std::strcat(pline, tmp_line);
                sfree(tmp_line);
            }

            /* skip trailing and leading spaces and comment text */
            strip_comment(pline);
            trim(pline);

            /* if there is something left... */
            if (static_cast<int>(std::strlen(pline)) > 0)
            {
                if (pline[0] == OPENDIR)
                {
                    /* A directive on this line: copy the directive
                     * without the brackets into dirstr, then
                     * skip spaces and tabs on either side of directive
                     */
                    dirstr = gmx_strdup((pline + 1));
                    if ((dummy2 = std::strchr(dirstr, CLOSEDIR)) != nullptr)
                    {
                        (*dummy2) = 0;
                    }
                    trim(dirstr);

                    if ((newd = str2dir(dirstr)) == Directive::d_invalid)
                    {
                        sprintf(errbuf, "Invalid directive %s", dirstr);
                        wi->addError(errbuf);
                    }
                    else
                    {
                        /* Directive found */
                        if (DS_Check_Order(DS, newd))
                        {
                            DS_Push(&DS, newd);
                            d = newd;
                        }
                        else
                        {
                            /* we should print here which directives should have
                               been present, and which actually are */
                            gmx_fatal(FARGS,
                                      "%s\nInvalid order for directive %s",
                                      cpp_error(&handle, eCPP_SYNTAX),
                                      enumValueToString(newd));
                            /* d = Directive::d_invalid; */
                        }

                        if (d == Directive::d_intermolecular_interactions)
                        {
                            if (*intermolecular_interactions == nullptr)
                            {
                                /* We (mis)use the moleculetype processing
                                 * to process the intermolecular interactions
                                 * by making a "molecule" of the size of the system.
                                 */
                                *intermolecular_interactions = std::make_unique<MoleculeInformation>();
                                mi0 = intermolecular_interactions->get();
                                mi0->initMolInfo();
                                make_atoms_sys(*molblock, *molinfo, &mi0->atoms);
                            }
                        }
                    }
                    sfree(dirstr);
                }
                else if (d != Directive::d_invalid)
                {
                    /* Not a directive, just a plain string
                     * use a gigantic switch to decode,
                     * if there is a valid directive!
                     */
                    switch (d)
                    {
                        case Directive::d_defaults:
                            if (bReadDefaults)
                            {
                                gmx_fatal(FARGS,
                                          "%s\nFound a second defaults directive.\n",
                                          cpp_error(&handle, eCPP_SYNTAX));
                            }
                            bReadDefaults = TRUE;
                            nscan         = sscanf(
                                    pline, "%s%s%s%lf%lf%lf", nb_str, comb_str, genpairs, &fLJ, &fQQ, &fPOW);
                            if (nscan < 2)
                            {
                                too_few(wi);
                            }
                            else
                            {
                                bGenPairs = FALSE;
                                fudgeLJ   = 1.0;
                                *fudgeQQ  = 1.0;

                                get_nbparm(nb_str, comb_str, &nb_funct, combination_rule, wi);
                                if (nscan >= 3)
                                {
                                    bGenPairs = (gmx::equalCaseInsensitive(genpairs, "Y", 1));
                                    if (nb_funct != VanDerWaalsPotential::LJ && bGenPairs)
                                    {
                                        gmx_fatal(FARGS,
                                                  "Generating pair parameters is only supported "
                                                  "with LJ non-bonded interactions");
                                    }
                                }
                                if (nscan >= 4)
                                {
                                    fudgeLJ = fLJ;
                                }
                                if (nscan >= 5)
                                {
                                    *fudgeQQ = fQQ;
                                }
                                if (nscan >= 6)
                                {
                                    *reppow = fPOW;
                                }
                            }
                            nb_funct = static_cast<VanDerWaalsPotential>(ifunc_index(
                                    Directive::d_nonbond_params, static_cast<int>(nb_funct)));

                            break;
                        case Directive::d_atomtypes:
                            push_at(atypes,
                                    &bondAtomType,
                                    pline,
                                    static_cast<int>(nb_funct),
                                    &nbparam,
                                    bGenPairs ? &pair : nullptr,
                                    wi);
                            break;

                        case Directive::d_bondtypes: // Intended to fall through
                        case Directive::d_constrainttypes:
                            push_bt(d, interactions, 2, nullptr, &bondAtomType, pline, wi);
                            break;
                        case Directive::d_pairtypes:
                            if (bGenPairs)
                            {
                                push_nbt(d,
                                         pair,
                                         atypes,
                                         pline,
                                         static_cast<int>(InteractionFunction::LennardJones14),
                                         wi);
                            }
                            else
                            {
                                push_bt(d, interactions, 2, atypes, nullptr, pline, wi);
                            }
                            break;
                        case Directive::d_angletypes:
                            push_bt(d, interactions, 3, nullptr, &bondAtomType, pline, wi);
                            break;
                        case Directive::d_dihedraltypes:
                            /* Special routine that can read both 2 and 4 atom dihedral definitions. */
                            push_dihedraltype(d, interactions, &bondAtomType, pline, wi);
                            break;

                        case Directive::d_nonbond_params:
                            push_nbt(d, nbparam, atypes, pline, static_cast<int>(nb_funct), wi);
                            break;

                        case Directive::d_implicit_genborn_params: // NOLINT bugprone-branch-clone
                            // Skip this line, so old topologies with
                            // GB parameters can be read.
                            break;

                        case Directive::d_implicit_surface_params:
                            // Skip this line, so that any topologies
                            // with surface parameters can be read
                            // (even though these were never formally
                            // supported).
                            break;

                        case Directive::d_cmaptypes:
                            push_cmaptype(d,
                                          interactions,
                                          NRAL(InteractionFunction::DihedralEnergyCorrectionMap),
                                          atypes,
                                          &bondAtomType,
                                          pline,
                                          wi);
                            break;

                        case Directive::d_moleculetype:
                        {
                            if (!bReadMolType)
                            {
                                int ntype;
                                if (opts->couple_moltype != nullptr
                                    && (opts->couple_lam0 == ecouplamNONE || opts->couple_lam0 == ecouplamQ
                                        || opts->couple_lam1 == ecouplamNONE
                                        || opts->couple_lam1 == ecouplamQ))
                                {
                                    dcatt = add_atomtype_decoupled(
                                            atypes, &nbparam, bGenPairs ? &pair : nullptr);
                                }
                                ntype  = atypes->size();
                                ncombs = (ntype * (ntype + 1)) / 2;
                                generate_nbparams(*combination_rule,
                                                  static_cast<int>(nb_funct),
                                                  &(interactions[static_cast<int>(nb_funct)]),
                                                  atypes,
                                                  wi);
                                ncopy = copy_nbparams(nbparam,
                                                      static_cast<int>(nb_funct),
                                                      &(interactions[static_cast<int>(nb_funct)]),
                                                      ntype);
                                GMX_LOG(logger.info)
                                        .asParagraph()
                                        .appendTextFormatted(
                                                "Generated %d of the %d non-bonded parameter "
                                                "combinations",
                                                ncombs - ncopy,
                                                ncombs);
                                free_nbparam(nbparam, ntype);
                                if (bGenPairs)
                                {
                                    gen_pairs((interactions[static_cast<int>(nb_funct)]),
                                              &(interactions[InteractionFunction::LennardJones14]),
                                              fudgeLJ,
                                              *combination_rule);
                                    ncopy = copy_nbparams(
                                            pair,
                                            static_cast<int>(nb_funct),
                                            &(interactions[InteractionFunction::LennardJones14]),
                                            ntype);
                                    GMX_LOG(logger.info)
                                            .asParagraph()
                                            .appendTextFormatted(
                                                    "Generated %d of the %d 1-4 parameter "
                                                    "combinations",
                                                    ncombs - ncopy,
                                                    ncombs);
                                    free_nbparam(pair, ntype);
                                }
                                /* Copy GBSA parameters to atomtype array? */

                                bReadMolType = TRUE;
                            }

                            push_molt(symtab, molinfo, pline, wi);
                            exclusionBlocks.emplace_back();
                            mi0                    = &molinfo->back();
                            mi0->atoms.haveMass    = TRUE;
                            mi0->atoms.haveCharge  = TRUE;
                            mi0->atoms.haveType    = TRUE;
                            mi0->atoms.haveBState  = TRUE;
                            mi0->atoms.havePdbInfo = FALSE;
                            break;
                        }
                        case Directive::d_atoms:
                            push_atom(symtab, &(mi0->atoms), atypes, pline, wi);
                            break;

                        case Directive::d_pairs:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      FALSE,
                                      bGenPairs,
                                      *fudgeQQ,
                                      bZero,
                                      false,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;
                        case Directive::d_pairs_nb:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      FALSE,
                                      FALSE,
                                      1.0,
                                      bZero,
                                      false,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;

                        case Directive::d_vsites1:
                        case Directive::d_vsites2:
                        case Directive::d_vsites3:
                        case Directive::d_vsites4:
                        case Directive::d_bonds:
                        case Directive::d_angles:
                        case Directive::d_constraints:
                        case Directive::d_settles:
                        case Directive::d_position_restraints:
                        case Directive::d_angle_restraints:
                        case Directive::d_angle_restraints_z:
                        case Directive::d_distance_restraints:
                        case Directive::d_orientation_restraints:
                        case Directive::d_dihedral_restraints:
                        case Directive::d_dihedrals:
                        case Directive::d_polarization:
                        case Directive::d_water_polarization:
                        case Directive::d_thole_polarization:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_bond(d,
                                      interactions,
                                      mi0->interactions,
                                      &(mi0->atoms),
                                      atypes,
                                      pline,
                                      TRUE,
                                      bGenPairs,
                                      *fudgeQQ,
                                      bZero,
                                      cpp_find_define(&handle, "_FF_AMBER_LEAP_ATOM_REORDERING") != nullptr,
                                      &bWarn_copy_A_B,
                                      wi);
                            break;
                        case Directive::d_cmap:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_cmap(d, interactions, mi0->interactions, &(mi0->atoms), atypes, pline, wi);
                            break;

                        case Directive::d_vsitesn:
                            GMX_RELEASE_ASSERT(
                                    mi0,
                                    "Need to have a valid MoleculeInformation object to work on");
                            push_vsitesn(d, mi0->interactions, &(mi0->atoms), pline, wi);
                            break;
                        case Directive::d_exclusions:
                            GMX_ASSERT(!exclusionBlocks.empty(),
                                       "exclusionBlocks must always be allocated so exclusions can "
                                       "be processed");
                            if (exclusionBlocks.back().empty())
                            {
                                GMX_RELEASE_ASSERT(mi0,
                                                   "Need to have a valid MoleculeInformation "
                                                   "object to work on");
                                exclusionBlocks.back().resize(mi0->atoms.nr);
                            }
                            push_excl(pline, exclusionBlocks.back(), wi);
                            break;
                        case Directive::d_system:
                            trim(pline);
                            title = put_symtab(symtab, pline);
                            break;
                        case Directive::d_molecules:
                        {
                            int  whichmol;
                            bool bCouple;

                            push_mol(*molinfo, pline, &whichmol, &nrcopies, wi);
                            mi0 = &((*molinfo)[whichmol]);
                            molblock->resize(molblock->size() + 1);
                            molblock->back().type = whichmol;
                            molblock->back().nmol = nrcopies;

                            bCouple = (opts->couple_moltype != nullptr
                                       && (gmx_strcasecmp("system", opts->couple_moltype) == 0
                                           || std::strcmp(*(mi0->name), opts->couple_moltype) == 0));
                            if (bCouple)
                            {
                                nmol_couple += nrcopies;
                            }

                            if (mi0->atoms.nr == 0)
                            {
                                gmx_fatal(FARGS, "Molecule type '%s' contains no atoms", *mi0->name);
                            }
                            GMX_LOG(logger.info)
                                    .asParagraph()
                                    .appendTextFormatted(
                                            "Excluding %d bonded neighbours molecule type '%s'",
                                            mi0->nrexcl,
                                            *mi0->name);

                            if (!mi0->bProcessed)
                            {
                                generate_excl(mi0->nrexcl, mi0->atoms.nr, mi0->interactions, &(mi0->excls));
                                gmx::mergeExclusions(&(mi0->excls), exclusionBlocks[whichmol]);
                                make_shake(mi0->interactions, &mi0->atoms, opts->nshake, logger);

                                if (bCouple)
                                {
                                    convert_moltype_couple(mi0,
                                                           dcatt,
                                                           *fudgeQQ,
                                                           opts->couple_lam0,
                                                           opts->couple_lam1,
                                                           opts->bCoupleIntra,
                                                           static_cast<int>(nb_funct),
                                                           &(interactions[static_cast<int>(nb_funct)]),
                                                           wi);
                                }
                                stupid_fill_block(&mi0->mols, mi0->atoms.nr, TRUE);
                                mi0->bProcessed = TRUE;
                            }

                            // After, potentially, applying decoupling we can accumulate the charge sum
                            sum_q(&mi0->atoms, nrcopies, &qt, &qBt);

                            break;
                        }
                        default:
                            if (d == Directive::d_intermolecular_interactions)
                            {
                                gmx_fatal(FARGS,
                                          "Expected a directive after directive '%s', not a line "
                                          "with: '%s'",
                                          enumValueToString(d),
                                          line);
                            }
                            else
                            {
                                GMX_RELEASE_ASSERT(
                                        false, "Unhandled combination of a line after a directive");
                            }
                    }
                }
            }
            sfree(pline);
            pline = nullptr;
        }
    } while (!done);

    // Check that all strings defined with -D were used when processing topology
    std::string unusedDefineWarning = checkAndWarnForUnusedDefines(*handle);
    if (!unusedDefineWarning.empty())
    {
        wi->addWarning(unusedDefineWarning);
    }

    for (char* element : cpp_opts_return)
    {
        sfree(element);
    }

    if (out)
    {
        gmx_fio_fclose(out);
    }

    /* List of GROMACS define names for force fields that have been
     * parametrized using constraints involving hydrogens only.
     *
     * We should avoid hardcoded names, but this is hopefully only
     * needed temparorily for discouraging use of constraints=all-bonds.
     */
    const std::array<std::string, 3> ffDefines = { "_FF_AMBER", "_FF_CHARMM", "_FF_OPLSAA" };
    *ffParametrizedWithHBondConstraints        = false;
    for (const std::string& ffDefine : ffDefines)
    {
        if (cpp_find_define(&handle, ffDefine))
        {
            *ffParametrizedWithHBondConstraints = true;
        }
    }

    if (cpp_find_define(&handle, "_FF_AMBER_LEAP_ATOM_REORDERING") != nullptr)
    {
        for (auto ftype : { InteractionFunction::ProperDihedrals,
                            InteractionFunction::RyckaertBellemansDihedrals,
                            InteractionFunction::ImproperDihedrals,
                            InteractionFunction::PeriodicImproperDihedrals })
        {
            GMX_ASSERT(interactions[ftype].leapDihedralTypes_.size()
                               == interactions[ftype].leapDihedralIndices_.size(),
                       "Numbers of AMBER LEaP dihedral types and their first occurrences should "
                       "match");
            if (!interactions[ftype].leapDihedralTypes_.empty())
            {
                GMX_LOG(logger.info)
                        .asParagraph()
                        .appendTextFormatted(
                                "To match AMBER LEaP ordering, reordered %zu and kept %zu %s "
                                "ordered as encountered (%zu %s types total)\n",
                                interactions[ftype].numLeapReorderingPerformed,
                                interactions[ftype].numLeapReorderingNotNecessary,
                                interaction_function[ftype].longname,
                                interactions[ftype].leapDihedralTypes_.size(),
                                interaction_function[ftype].longname);
            }
        }
    }

    if (cpp_find_define(&handle, "_FF_GROMOS96") != nullptr)
    {
        wi->addWarning(
                "The GROMOS force fields have been parametrized with a physically incorrect "
                "multiple-time-stepping scheme for a twin-range cut-off. When used with "
                "a single-range cut-off (or a correct Trotter multiple-time-stepping scheme), "
                "physical properties, such as the density, might differ from the intended values. "
                "Since there are researchers actively working on validating GROMOS with modern "
                "integrators we have not yet removed the GROMOS force fields, but you should be "
                "aware of these issues and check if molecules in your system are affected before "
                "proceeding. "
                "Further information is available at "
                "https://gitlab.com/gromacs/gromacs/-/issues/2884, "
                "and a longer explanation of our decision to remove physically incorrect "
                "algorithms "
                "can be found at https://doi.org/10.26434/chemrxiv.11474583.v1 .");
    }
    // TODO: Update URL for Issue #2884 in conjunction with updating grompp.warn in regressiontests.

    cpp_done(handle);

    if (opts->couple_moltype)
    {
        if (nmol_couple == 0)
        {
            gmx_fatal(FARGS, "Did not find any molecules of type '%s' for coupling", opts->couple_moltype);
        }
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted(
                        "Coupling %d copies of molecule type '%s'", nmol_couple, opts->couple_moltype);
    }

    /* this is not very clean, but fixes core dump on empty system name */
    if (!title)
    {
        title = put_symtab(symtab, "");
    }

    if (std::fabs(qt) > 1e-4)
    {
        sprintf(warn_buf, "System has non-zero total charge: %.6f\n%s\n", qt, floating_point_arithmetic_tip);
        wi->addNote(warn_buf);
    }
    if (std::fabs(qBt) > 1e-4 && !gmx_within_tol(qBt, qt, 1e-6))
    {
        sprintf(warn_buf, "State B has non-zero total charge: %.6f\n%s\n", qBt, floating_point_arithmetic_tip);
        wi->addNote(warn_buf);
    }
    if (usingFullRangeElectrostatics && (std::fabs(qt) > 1e-4 || std::fabs(qBt) > 1e-4))
    {
        wi->addWarning(
                "You are using Ewald electrostatics in a system with net charge. This can lead to "
                "severe artifacts, such as ions moving into regions with low dielectric, due to "
                "the uniform background charge. We suggest to neutralize your system with counter "
                "ions, possibly in combination with a physiological salt concentration.");
        please_cite(stdout, "Hub2014a");
    }

    DS_Done(&DS);

    if (*intermolecular_interactions != nullptr)
    {
        sfree(intermolecular_interactions->get()->atoms.atom);
    }

    return title;
}

char** do_top(bool                                                            bVerbose,
              const char*                                                     topfile,
              const std::optional<std::filesystem::path>&                     topppfile,
              t_gromppopts*                                                   opts,
              bool                                                            bZero,
              t_symtab*                                                       symtab,
              gmx::EnumerationArray<InteractionFunction, InteractionsOfType>& interactions,
              CombinationRule*                                                combination_rule,
              double*                                                         repulsion_power,
              real*                                                           fudgeQQ,
              PreprocessingAtomTypes*                                         atypes,
              std::vector<MoleculeInformation>*                               molinfo,
              std::unique_ptr<MoleculeInformation>* intermolecular_interactions,
              const t_inputrec*                     ir,
              std::vector<gmx_molblock_t>*          molblock,
              bool*                                 ffParametrizedWithHBondConstraints,
              WarningHandler*                       wi,
              const gmx::MDLogger&                  logger)
{
    char** title;

    if (bVerbose)
    {
        GMX_LOG(logger.info).asParagraph().appendTextFormatted("processing topology...");
    }
    title = read_topol(topfile,
                       topppfile,
                       opts->define,
                       opts->include,
                       symtab,
                       atypes,
                       molinfo,
                       intermolecular_interactions,
                       interactions,
                       combination_rule,
                       repulsion_power,
                       opts,
                       fudgeQQ,
                       molblock,
                       ffParametrizedWithHBondConstraints,
                       ir->efep != FreeEnergyPerturbationType::No,
                       bZero,
                       usingFullElectrostatics(ir->coulombtype),
                       wi,
                       logger);

    if ((*combination_rule != CombinationRule::Geometric) && (ir->vdwtype == VanDerWaalsType::User))
    {
        wi->addWarning(
                "Using sigma/epsilon based combination rules with"
                " user supplied potential function may produce unwanted"
                " results");
    }

    return title;
}

/*! \brief Bookkeeping of the force-field terms dropped by generate_qmexcl_moltype().
 *
 * Counted per interaction type and split by how the atoms of the term are
 * distributed over the two regions, so that the effect of the chosen scheme
 * on the QM/MM boundary can be read off the grompp output.
 */
struct QmmmRemovedInteractions
{
    //! Terms all of whose atoms are QM -- described by the QM calculation itself.
    gmx::EnumerationArray<InteractionFunction, int> allQm = {};
    //! Terms with atoms in both regions -- the QM/MM boundary terms.
    gmx::EnumerationArray<InteractionFunction, int> boundary = {};
    //! Terms without a single QM atom.
    gmx::EnumerationArray<InteractionFunction, int> mmOnly = {};
};

/*! \brief Per-atom record of what generate_qmexcl_moltype() changed.
 *
 * Atom numbers are global and 0-based here; they are written 1-based, i.e. in the
 * numbering of the input coordinate file. The labels are stored together with the
 * numbers, because the molecule types are modified while the report is filled.
 */
struct QmmmTopologyReport
{
    //! One atom of a reported term: global index, "RESnr NAME" label, QM or not
    struct Atom
    {
        int         index;
        std::string label;
        bool        isQm;
    };
    //! A force-field term (bonded interaction or pair) with its atoms
    struct Term
    {
        InteractionFunction ftype;
        std::vector<Atom>   atoms;
    };
    //! A bond with exactly one QM atom, and whether it makes its MM atom a boundary atom
    struct BoundaryBond
    {
        InteractionFunction ftype;
        Atom                qm;
        Atom                mm;
        bool                accepted;
        std::string         reason;
    };
    //! A nonbonded exclusion between a QM atom and another atom
    struct Exclusion
    {
        Atom qm;
        Atom other;
        int  bondDistance;  //!< number of bonds between the two atoms, -1 if more than 7
        bool presentBefore; //!< already excluded by the force field (nrexcl), or a link atom
    };
    std::vector<Atom>         qmAtoms;
    std::vector<BoundaryBond> boundaryBonds;
    std::vector<Term>         removedBonded;
    std::vector<Term>         keptRestraints;
    std::vector<Term>         convertedBonds;
    std::vector<Term>         removedPairs;
    std::vector<Exclusion>    qmQmExclusions;
    std::vector<Exclusion>    qmMmExclusions;
    //! QM--MM exclusions of the final topology; presentBefore marks a pair with a link atom
    std::vector<Exclusion> finalQmMmExclusions;
    //! QM--MM LJ-14 pairs kept in the final topology; presentBefore marks a link atom
    std::vector<Exclusion> keptQmMmPairs;
    //! Whether the LJ and LJ-14 of the boundary MM atoms with all QM atoms are excluded
    bool excludeBoundaryLJ = false;
};

//! "RESnr NAME" label of a local atom of a molecule type
static std::string qmmmAtomLabel(const gmx_moltype_t& molt, int localIndex)
{
    const t_atoms&   atoms = molt.atoms;
    const t_resinfo& ri    = atoms.resinfo[atoms.atom[localIndex].resind];
    return gmx::formatString("%s%d %s", *ri.name, ri.nr, *atoms.atomname[localIndex]);
}

//! Bond distances (up to \p maxDepth) from local atom \p start over the chemical bonds, -1 beyond
static std::vector<int> qmmmBondDistances(const std::vector<std::vector<int>>& graph, int start, int maxDepth)
{
    std::vector<int> depth(graph.size(), -1);
    std::vector<int> frontier{ start };
    depth[start] = 0;
    for (int d = 1; d <= maxDepth && !frontier.empty(); d++)
    {
        std::vector<int> next;
        for (int a : frontier)
        {
            for (int b : graph[a])
            {
                if (depth[b] == -1)
                {
                    depth[b] = d;
                    next.push_back(b);
                }
            }
        }
        frontier = std::move(next);
    }
    return depth;
}

//! Name of a QM/MM mode for the output
static const char* qmmmSchemeName(QmmmModeType qmmmMode)
{
    switch (qmmmMode)
    {
        case QmmmModeType::MiMiC: return "MiMiC";
        case QmmmModeType::Amber: return "amber";
        default: return "classic";
    }
}

/*! \brief
 * Exclude molecular interactions for QM atoms in QM/MM
 *
 * Update the exclusion lists to include all QM atoms of this molecule,
 * replace bonds between QM atoms with InteractionFunction::ConnectBonds and
 * set charges of QM atoms to 0.
 *
 * The bonded interactions that involve both QM and MM atoms are treated
 * according to one of two conventions, selected with \p qmmmMode
 * (GMX_QMMM_BONDED_SCHEME in grompp):
 *
 * QmmmModeType::Original ("classic", the default): a bonded interaction is
 *   removed as soon as all but one of its atoms are QM (a QM-QM-MM angle and a
 *   QM-QM-QM-MM dihedral are removed), because the QM calculation with the link
 *   atom describes it and keeping the force-field term would count it twice.
 *   The LJ between QM and MM atoms follows the exclusion rules of the force field
 *   (nrexcl and [ pairs ]). With \p excludeBoundaryLJ (GMX_QMMM_LJ_SCHEME=exclude)
 *   the LJ and LJ-14 of every QM atom with the MM atoms covalently bound to the QM
 *   region are excluded as well.
 *
 * QmmmModeType::Amber ("amber"): only interactions whose atoms are all QM are
 *   removed; every term with at least one MM atom is kept at the force-field
 *   level, as in the QM/MM implementation of AMBER. No exclusions of the boundary
 *   MM atoms are generated.
 *
 * QmmmModeType::MiMiC uses the rule of the amber scheme for the bonded terms,
 *   because MiMiC treats the link atoms as quantum atoms.
 *
 * Restraints (position, flat-bottomed position, distance, orientation, angle,
 *   dihedral restraints and restraint potentials) are kept with every scheme, also
 *   on QM atoms: they are set by the user and not part of the force field that the
 *   QM calculation replaces.
 *
 * \param[in,out] molt molecule type with QM atoms
 * \param[in] grpnr group informatio
 * \param[in,out] ir input record
 * \param[in] qmmmMode QM/MM mode: classic (Original), amber or MiMiC
 * \param[in] logger Handle to logging interface.
 * \param[in,out] removed Counters of the removed interactions, summed over the molecule types
 * \param[in] atomOffset Global index of the first atom of this molecule, for the report
 * \param[in,out] report Per-atom record of the changes
 * \param[in] excludeBoundaryLJ Exclude the LJ and LJ-14 of the boundary MM atoms with every
 *                              QM atom (classic scheme only)
 */
static void generate_qmexcl_moltype(gmx_moltype_t*           molt,
                                    const unsigned char*     grpnr,
                                    t_inputrec*              ir,
                                    QmmmModeType             qmmmMode,
                                    const gmx::MDLogger&     logger,
                                    QmmmRemovedInteractions* removed,
                                    int                      atomOffset,
                                    QmmmTopologyReport*      report,
                                    bool                     excludeBoundaryLJ)
{
    /* This routine expects molt->ilist to be of size InteractionFunction::Count and ordered. */

    /* generates the exclusions between the individual QM atoms, as
     * these interactions should be handled by the QM subroutines and
     * not by the gromacs routines
     */
    int   qm_max = 0, qm_nr = 0, link_nr = 0, link_max = 0;
    int * qm_arr = nullptr, *link_arr = nullptr;
    bool *bQMMM, *blink;

    /* First we search and select the QM atoms in an qm_arr array that
     * we use to create the exclusions.
     *
     * we take the possibility into account that a user has defined more
     * than one QM group:
     *
     * for that we also need to do this an ugly work-about just in case
     * the QM group contains the entire system...
     */

    /* we first search for all the QM atoms and put them in an array
     */
    for (int j = 0; j < ir->opts.ngQM; j++)
    {
        for (int i = 0; i < molt->atoms.nr; i++)
        {
            if (qm_nr >= qm_max)
            {
                qm_max += 100;
                srenew(qm_arr, qm_max);
            }
            if ((grpnr ? grpnr[i] : 0) == j)
            {
                qm_arr[qm_nr++]        = i;
                molt->atoms.atom[i].q  = 0.0;
                molt->atoms.atom[i].qB = 0.0;
            }
        }
    }
    /* bQMMM[..] is an array containin TRUE/FALSE for atoms that are
     * QM/not QM. We first set all elements to false. Afterwards we use
     * the qm_arr to change the elements corresponding to the QM atoms
     * to TRUE.
     */
    snew(bQMMM, molt->atoms.nr);
    for (int i = 0; i < molt->atoms.nr; i++)
    {
        bQMMM[i] = FALSE;
    }
    for (int i = 0; i < qm_nr; i++)
    {
        bQMMM[qm_arr[i]] = TRUE;
    }

    const auto reportAtom = [&](int local) {
        return QmmmTopologyReport::Atom{ atomOffset + local, qmmmAtomLabel(*molt, local),
                                         static_cast<bool>(bQMMM[local]) };
    };
    const auto reportTerm = [&](InteractionFunction ftype, const int* iatoms, int nratoms) {
        QmmmTopologyReport::Term term{ ftype, {} };
        for (int k = 0; k < nratoms; k++)
        {
            term.atoms.push_back(reportAtom(iatoms[k]));
        }
        return term;
    };
    for (int i = 0; i < qm_nr; i++)
    {
        report->qmAtoms.push_back(reportAtom(qm_arr[i]));
    }

    /* The atoms each virtual site is constructed from. A link atom is a virtual
     * site in the QM group, and the bonds that start from it are connections
     * (funct 5) that only serve to generate exclusions.
     */
    std::vector<std::vector<int>> vsiteConstructingAtoms(molt->atoms.nr);
    for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        if (!(interaction_function[ftype].flags & IF_VSITE))
        {
            continue;
        }
        const int              nratoms = interaction_function[ftype].nratoms;
        const InteractionList& il      = molt->ilist[ftype];
        for (int i = 0; i < il.size(); i += 1 + nratoms)
        {
            for (int k = 2; k <= nratoms; k++)
            {
                vsiteConstructingAtoms[il.iatoms[i + 1]].push_back(il.iatoms[i + k]);
            }
        }
    }

    /* Chemical-bond graph of the molecule, for the bond distances in the report */
    std::vector<std::vector<int>> bondGraph(molt->atoms.nr);
    for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        if (IS_CHEMBOND(ftype))
        {
            const InteractionList& il = molt->ilist[ftype];
            for (int i = 0; i < il.size(); i += 3)
            {
                bondGraph[il.iatoms[i + 1]].push_back(il.iatoms[i + 2]);
                bondGraph[il.iatoms[i + 2]].push_back(il.iatoms[i + 1]);
            }
        }
    }

    /* We remove the bonded interactions (i.e. bonds, angles, dihedrals,
     * 1-4's) that are described by the QM calculation. An interaction between
     * two atoms is removed if both atoms are QM atoms. An interaction of three
     * or more atoms is removed, with the classic scheme, if at most one of its
     * atoms is an MM atom, and with the amber scheme (and MiMiC) only if all of
     * its atoms are QM. Restraints are never removed. Since this routine is
     * called once before any forces are computed, the top->idef.il[N].iatom[]
     * array (see idef.h) can be rewritten at this point without any problem.
     */

    /* first check whether we already have CONNBONDS.
     * Note that if we don't, we don't add a param entry and set ftype=0,
     * which is ok, since CONNBONDS does not use parameters.
     */
    int ftype_connbond = 0;
    int ind_connbond   = 0;
    if (!molt->ilist[InteractionFunction::ConnectBonds].empty())
    {
        GMX_LOG(logger.info)
                .asParagraph()
                .appendTextFormatted("nr. of CONNBONDS present already: %d",
                                     molt->ilist[InteractionFunction::ConnectBonds].size() / 3);
        ftype_connbond = molt->ilist[InteractionFunction::ConnectBonds].iatoms[0];
        ind_connbond   = molt->ilist[InteractionFunction::ConnectBonds].size();
    }
    /* now we delete all bonded interactions, except the ones describing
     * a chemical bond. These are converted to CONNBONDS
     */
    for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        if (!(interaction_function[ftype].flags & IF_BOND) || ftype == InteractionFunction::ConnectBonds)
        {
            continue;
        }
        int nratoms = interaction_function[ftype].nratoms;
        if (IS_RESTRAINT_TYPE(ftype))
        {
            /* Restraints are set by the user and are not part of the force field
             * that the QM calculation replaces: keep them, also on QM atoms.
             */
            const InteractionList& il = molt->ilist[ftype];
            for (int j = 0; j < il.size(); j += nratoms + 1)
            {
                bool hasQm = false;
                for (int k = 0; k < nratoms; k++)
                {
                    hasQm = hasQm || bQMMM[il.iatoms[j + 1 + k]];
                }
                if (hasQm)
                {
                    report->keptRestraints.push_back(reportTerm(ftype, il.iatoms.data() + j + 1, nratoms));
                }
            }
            continue;
        }
        int j = 0;
        while (j < molt->ilist[ftype].size())
        {
            bool bexcl;

            int numQmAtoms = 0;
            for (int jj = j + 1; jj < j + 1 + nratoms; jj++)
            {
                if (bQMMM[molt->ilist[ftype].iatoms[jj]])
                {
                    numQmAtoms++;
                }
            }

            if (nratoms <= 2)
            {
                /* Remove an interaction of one or two atoms when all of them are
                 * in the QM region. Note that we don't have to worry about
                 * link atoms here, as they won't have 2-atom interactions.
                 */
                bexcl = (numQmAtoms == nratoms);
                /* A chemical bond between two QM atoms will be copied to
                 * the InteractionFunction::ConnectBonds list, for reasons mentioned above.
                 */
                if (bexcl && IS_CHEMBOND(ftype))
                {
                    int a1 = molt->ilist[ftype].iatoms[1 + j + 0];
                    int a2 = molt->ilist[ftype].iatoms[1 + j + 1];
                    report->convertedBonds.push_back(
                            reportTerm(ftype, molt->ilist[ftype].iatoms.data() + j + 1, nratoms));
                    InteractionList& ilist = molt->ilist[InteractionFunction::ConnectBonds];
                    ilist.iatoms.resize(ind_connbond + 3);
                    ilist.iatoms[ind_connbond++] = ftype_connbond;
                    ilist.iatoms[ind_connbond++] = a1;
                    ilist.iatoms[ind_connbond++] = a2;
                }
            }
            else
            {
                /* With the classic scheme, MM interactions have to be excluded if
                 * they are included in the QM already. Because we use a link atom
                 * (H atom) when the QM/MM boundary runs through a chemical bond,
                 * this means that as long as one atom is MM, we still exclude, as
                 * the interaction is included in the QM via:
                 * QMatom1-QMatom2-QMatom-3-Linkatom.
                 * MiMiC treats link atoms as quantum atoms, and the amber scheme
                 * keeps every term with an MM atom at the force-field level, so
                 * these only remove interactions all of whose atoms are QM.
                 */
                if (qmmmMode == QmmmModeType::Original)
                {
                    bexcl = (numQmAtoms >= nratoms - 1);
                }
                else
                {
                    bexcl = (numQmAtoms == nratoms);
                }

                if (bexcl && ftype == InteractionFunction::SETTLE)
                {
                    gmx_fatal(FARGS,
                              "Can not apply QM to molecules with SETTLE, replace the moleculetype "
                              "using QM and SETTLE by one without SETTLE");
                }
            }
            if (bexcl)
            {
                /* keep track of what is being removed, for the report at the end */
                if (interaction_function[ftype].flags & IF_PAIR)
                {
                    report->removedPairs.push_back(
                            reportTerm(ftype, molt->ilist[ftype].iatoms.data() + j + 1, nratoms));
                }
                else if (!(nratoms == 2 && IS_CHEMBOND(ftype)))
                {
                    report->removedBonded.push_back(
                            reportTerm(ftype, molt->ilist[ftype].iatoms.data() + j + 1, nratoms));
                }
                if (numQmAtoms == nratoms)
                {
                    removed->allQm[ftype]++;
                }
                else if (numQmAtoms > 0)
                {
                    removed->boundary[ftype]++;
                }
                else
                {
                    removed->mmOnly[ftype]++;
                }
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList& ilist = molt->ilist[ftype];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms + 1; /* the +1 is for the functype */
            }
        }
    }
    /* Now, we search for atoms bonded to a QM atom because we may also want
     * to exclude their nonbonded interactions with the QM atoms. The
     * reason for this is that this interaction is accounted for in the
     * linkatoms interaction with the QMatoms and would be counted
     * twice. This is only done with the classic scheme and on request
     * (GMX_QMMM_LJ_SCHEME=exclude); the boundary atoms are collected and
     * reported in any case.
     */
    snew(blink, molt->atoms.nr);
    for (int i = 0; i < molt->atoms.nr; i++)
    {
        blink[i] = FALSE;
    }
    if (qmmmMode == QmmmModeType::Original)
    {
        for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
        {
            if (IS_CHEMBOND(ftype))
            {
                int j = 0;
                while (j < molt->ilist[ftype].size())
                {
                    int a1 = molt->ilist[ftype].iatoms[j + 1];
                    int a2 = molt->ilist[ftype].iatoms[j + 2];
                    if ((bQMMM[a1] && !bQMMM[a2]) || (!bQMMM[a1] && bQMMM[a2]))
                    {
                        const int qmAtom = bQMMM[a1] ? a1 : a2;
                        const int mmAtom = bQMMM[a1] ? a2 : a1;
                        /* A bond from a link atom (a virtual site of the QM group) is
                         * a connection that only generates exclusions. It marks a
                         * boundary MM atom only if that atom is one the link atom is
                         * constructed from, i.e. the MM atom of the cut bond. Bonds
                         * from a link atom to the further neighbours (MM2) must not
                         * extend the LJ and LJ-14 exclusions to those atoms.
                         */
                        const std::vector<int>& constructing = vsiteConstructingAtoms[qmAtom];
                        const bool              fromLinkAtom = !constructing.empty();
                        const bool              accepted =
                                !fromLinkAtom
                                || std::find(constructing.begin(), constructing.end(), mmAtom)
                                           != constructing.end();
                        report->boundaryBonds.push_back(
                                { ftype, reportAtom(qmAtom), reportAtom(mmAtom), accepted,
                                  !fromLinkAtom ? "bond of a QM atom"
                                                : (accepted ? "link atom constructed from this MM atom"
                                                            : "link atom NOT constructed from this MM atom: ignored") });
                        if (accepted)
                        {
                            if (link_nr >= link_max)
                            {
                                link_max += 10;
                                srenew(link_arr, link_max);
                            }
                            link_arr[link_nr++] = mmAtom;
                        }
                    }
                    j += 3;
                }
            }
        }
    }
    /* The boundary MM atoms lose their LJ and LJ-14 interactions with every QM
     * atom only on request. By default the LJ between the QM and the MM atoms
     * follows the exclusion rules of the force field: the pairs within nrexcl
     * bonds are excluded by grompp anyway, the 1-4 pairs keep their LJ-14, and
     * every other pair keeps its LJ.
     */
    const int numLinkExclusions = (qmmmMode == QmmmModeType::Original && excludeBoundaryLJ) ? link_nr : 0;
    for (int i = 0; i < numLinkExclusions; i++)
    {
        blink[link_arr[i]] = TRUE;
    }
    /* creating the exclusion block for the QM atoms. Each QM atom has
     * as excluded elements all the other QMatoms (and itself).
     */
    t_blocka qmexcl;
    qmexcl.nr  = molt->atoms.nr;
    qmexcl.nra = qm_nr * (qm_nr + numLinkExclusions) + numLinkExclusions * qm_nr;
    snew(qmexcl.index, qmexcl.nr + 1);
    snew(qmexcl.a, qmexcl.nra);
    int l = 0;
    for (int i = 0; i < qmexcl.nr; i++)
    {
        qmexcl.index[i] = l;
        if (bQMMM[i])
        {
            for (int k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k + l] = qm_arr[k];
            }
            for (int k = 0; k < numLinkExclusions; k++)
            {
                qmexcl.a[qm_nr + k + l] = link_arr[k];
            }
            l += (qm_nr + numLinkExclusions);
        }
        if (blink[i])
        {
            for (int k = 0; k < qm_nr; k++)
            {
                qmexcl.a[k + l] = qm_arr[k];
            }
            l += qm_nr;
        }
    }
    qmexcl.index[qmexcl.nr] = l;

    /* record the exclusions for the report, and whether the force field
     * (nrexcl) had excluded the pair already
     */
    {
        const auto excludedBefore = [&](int a, int b) {
            const auto list = molt->excls[a];
            return std::find(list.begin(), list.end(), b) != list.end();
        };
        std::vector<bool> reported(molt->atoms.nr, false);
        for (int a = 0; a < qm_nr; a++)
        {
            const int              qa    = qm_arr[a];
            const std::vector<int> depth = qmmmBondDistances(bondGraph, qa, 7);
            for (int b = a + 1; b < qm_nr; b++)
            {
                const int qb = qm_arr[b];
                report->qmQmExclusions.push_back(
                        { reportAtom(qa), reportAtom(qb), depth[qb], excludedBefore(qa, qb) });
            }
            std::fill(reported.begin(), reported.end(), false);
            for (int k = 0; k < numLinkExclusions; k++)
            {
                const int mm = link_arr[k];
                if (reported[mm])
                {
                    continue; // the same MM atom can be reached by more than one bond
                }
                reported[mm] = true;
                report->qmMmExclusions.push_back(
                        { reportAtom(qa), reportAtom(mm), depth[mm], excludedBefore(qa, mm) });
            }
        }
    }

    /* and merging with the exclusions already present in sys.
     */

    std::vector<gmx::ExclusionBlock> qmexcl2(molt->atoms.nr);
    gmx::blockaToExclusionBlocks(&qmexcl, qmexcl2);
    gmx::mergeExclusions(&(molt->excls), qmexcl2);
    sfree(qmexcl.index);
    sfree(qmexcl.a);

    /* Finally, we also need to get rid of the pair interactions of the
     * classical atom bonded to the boundary QM atoms with the QMatoms,
     * as this interaction is already accounted for by the QM, so also
     * here we run the risk of double counting! We proceed in a similar
     * way as we did above for the other bonded interactions.
     * (Unless the boundary LJ is excluded, blink[] is false everywhere,
     * and only the pairs between two QM atoms are removed here.) */
    for (InteractionFunction i : { InteractionFunction::LennardJones14, InteractionFunction::Coulomb14 })
    {
        int nratoms = interaction_function[i].nratoms;
        int j       = 0;
        while (j < molt->ilist[i].size())
        {
            int  a1 = molt->ilist[i].iatoms[j + 1];
            int  a2 = molt->ilist[i].iatoms[j + 2];
            bool bexcl =
                    ((bQMMM[a1] && bQMMM[a2]) || (blink[a1] && bQMMM[a2]) || (bQMMM[a1] && blink[a2]));
            if (bexcl)
            {
                report->removedPairs.push_back(reportTerm(i, molt->ilist[i].iatoms.data() + j + 1, 2));
                if (bQMMM[a1] && bQMMM[a2])
                {
                    removed->allQm[i]++;
                }
                else
                {
                    removed->boundary[i]++;
                }
                /* since the interaction involves QM atoms, these should be
                 * removed from the MM ilist
                 */
                InteractionList& ilist = molt->ilist[i];
                for (int k = j; k < ilist.size() - (nratoms + 1); k++)
                {
                    ilist.iatoms[k] = ilist.iatoms[k + (nratoms + 1)];
                }
                ilist.iatoms.resize(ilist.size() - (nratoms + 1));
            }
            else
            {
                j += nratoms + 1; /* the +1 is for the functype */
            }
        }
    }

    /* the QM--MM exclusions and LJ-14 pairs that end up in the tpr, for the report */
    {
        const auto isLinkAtom = [&](int a) { return !vsiteConstructingAtoms[a].empty(); };
        for (int a = 0; a < molt->atoms.nr; a++)
        {
            if (!bQMMM[a])
            {
                continue;
            }
            const std::vector<int> depth = qmmmBondDistances(bondGraph, a, 7);
            std::vector<int>       partners;
            for (int b : molt->excls[a])
            {
                if (!bQMMM[b])
                {
                    partners.push_back(b);
                }
            }
            std::sort(partners.begin(), partners.end());
            for (int b : partners)
            {
                report->finalQmMmExclusions.push_back(
                        { reportAtom(a), reportAtom(b), depth[b], isLinkAtom(a) });
            }
        }
        const InteractionList& il = molt->ilist[InteractionFunction::LennardJones14];
        for (int k = 0; k < il.size(); k += 3)
        {
            const int a1 = il.iatoms[k + 1], a2 = il.iatoms[k + 2];
            if (bQMMM[a1] != bQMMM[a2])
            {
                const int qa = bQMMM[a1] ? a1 : a2, mb = bQMMM[a1] ? a2 : a1;
                report->keptQmMmPairs.push_back({ reportAtom(qa), reportAtom(mb),
                                                  qmmmBondDistances(bondGraph, qa, 7)[mb],
                                                  isLinkAtom(qa) });
            }
        }
    }
    report->excludeBoundaryLJ = (numLinkExclusions > 0) || (qmmmMode == QmmmModeType::Original && excludeBoundaryLJ);

    std::free(qm_arr);
    std::free(bQMMM);
    std::free(link_arr);
    std::free(blink);
} /* generate_qmexcl_moltype */

/*! \brief Report the force-field terms that were removed around the QM region.
 *
 * The interesting column is the middle one: those are the terms that connect the
 * two regions, and they are the ones whose treatment differs between the classic
 * and the amber scheme. The terms with only QM atoms are removed by either scheme,
 * because the QM calculation describes them.
 */
static void reportQmmmRemovedInteractions(const QmmmRemovedInteractions& removed,
                                          const QmmmTopologyReport&      report,
                                          QmmmModeType                   qmmmMode,
                                          const gmx::MDLogger&           logger)
{
    int totalAllQm = 0, totalBoundary = 0, totalMmOnly = 0;
    for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        totalAllQm += removed.allQm[ftype];
        totalBoundary += removed.boundary[ftype];
        totalMmOnly += removed.mmOnly[ftype];
    }

    GMX_LOG(logger.info)
            .appendTextFormatted(
                    "\nQM/MM: force-field terms removed with the '%s' scheme, by interaction type:",
                    qmmmSchemeName(qmmmMode));
    if (totalAllQm + totalBoundary + totalMmOnly == 0)
    {
        GMX_LOG(logger.info).appendTextFormatted("  (none)");
    }
    else
    {
        GMX_LOG(logger.info)
                .appendTextFormatted("  %-22s %10s %10s %10s", "interaction", "all-QM", "QM--MM", "MM-only");
        for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
        {
            if (removed.allQm[ftype] + removed.boundary[ftype] + removed.mmOnly[ftype] == 0)
            {
                continue;
            }
            GMX_LOG(logger.info)
                    .appendTextFormatted("  %-22s %10d %10d %10d", interaction_function[ftype].longname,
                                         removed.allQm[ftype], removed.boundary[ftype],
                                         removed.mmOnly[ftype]);
        }
        GMX_LOG(logger.info)
                .appendTextFormatted("  %-22s %10d %10d %10d", "total", totalAllQm, totalBoundary, totalMmOnly);
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "  all-QM  = every atom of the term is QM: described by the QM calculation\n"
                        "  QM--MM  = the term spans the QM/MM boundary\n"
                        "  MM-only = no QM atom in the term at all");
    }
    if (totalBoundary == 0)
    {
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "No term spanning the boundary was removed: every force-field term with at "
                        "least one MM atom is kept.");
    }
    else
    {
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "The QM--MM terms above are assumed to be described by the QM calculation "
                        "with the link atom;\nswitch the scheme with GMX_QMMM_BONDED_SCHEME=amber to "
                        "keep them at the force-field level.");
    }
    GMX_LOG(logger.info)
            .appendTextFormatted(
                    "Restraints kept on QM atoms: %zu. Chemical bonds between two QM atoms are not "
                    "lost but converted to\nconnections (%zu), and the removed pair interactions "
                    "(LJ-14) are counted above as well.\n",
                    report.keptRestraints.size(), report.convertedBonds.size());
}

/*! \brief Whether the per-atom QM/MM report files are written.
 *
 * On by default; GMX_QMMM_REPORTS set to 0, no, off or false switches off the
 * reports of both grompp and mdrun, together with the lines that point to them.
 */
static bool qmmmReportsEnabled()
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

//! Writes the detailed per-atom report of the QM/MM changes to the topology
static void writeQmmmTopologyReport(const QmmmTopologyReport& report, QmmmModeType qmmmMode, const char* fileName)
{
    FILE* fp = std::fopen(fileName, "w");
    if (fp == nullptr)
    {
        return;
    }
    const auto atomText = [](const QmmmTopologyReport::Atom& a) {
        return gmx::formatString("%7d %-14s %s", a.index + 1, a.label.c_str(), a.isQm ? "QM" : "MM");
    };
    const auto distanceText = [](int d) {
        return d < 0 ? std::string("> 1-8") : gmx::formatString("1-%d", d + 1);
    };
    const auto termText = [&atomText](const QmmmTopologyReport::Term& term) {
        std::string line;
        for (const auto& a : term.atoms)
        {
            line += " |" + atomText(a);
        }
        return line;
    };

    std::fprintf(fp, "; QM/MM changes to the force-field topology, written by gmx grompp\n");
    std::fprintf(fp, "; scheme for the bonded terms at the QM/MM boundary: %s (GMX_QMMM_BONDED_SCHEME)\n",
                 qmmmSchemeName(qmmmMode));
    std::fprintf(fp, "; LJ between QM and MM atoms: %s (GMX_QMMM_LJ_SCHEME)\n",
                 report.excludeBoundaryLJ
                         ? "exclude -- boundary MM atoms lose LJ and LJ-14 with every QM atom"
                         : "forcefield -- exclusion rules of the force field (nrexcl, [ pairs ])");
    std::fprintf(fp, "; atom numbers are global and 1-based, i.e. the numbering of the input .gro file\n");
    std::fprintf(fp, "; labels are RESIDUEnumber ATOMNAME from the topology; QM/MM marks the region\n");
    std::fprintf(fp, "; the QM/MM electrostatics at the boundary are set by mdrun, not here\n\n");

    std::fprintf(fp, "[ qm_atoms ]\n; %zu atoms; their charges are set to zero in the tpr\n",
                 report.qmAtoms.size());
    for (const auto& a : report.qmAtoms)
    {
        std::fprintf(fp, "%s\n", atomText(a).c_str());
    }

    std::fprintf(fp, "\n[ boundary_bonds ]\n");
    std::fprintf(fp, "; chemical bonds and connections with exactly one QM atom (classic scheme). The MM\n");
    std::fprintf(fp, "; atom of an accepted bond is a boundary MM atom. Only with GMX_QMMM_LJ_SCHEME=exclude\n");
    std::fprintf(fp, "; its LJ and LJ-14 interactions with every QM atom are excluded; %s.\n",
                 report.excludeBoundaryLJ ? "this is the case here" : "not the case here");
    std::fprintf(fp, "; %-26s %-26s %-12s %s\n", "QM atom", "MM atom", "bond type", "boundary MM atom?");
    for (const auto& b : report.boundaryBonds)
    {
        std::fprintf(fp, "%s %s  %-12s %s (%s)\n", atomText(b.qm).c_str(), atomText(b.mm).c_str(),
                     interaction_function[b.ftype].name, b.accepted ? "yes" : "no", b.reason.c_str());
    }
    if (report.boundaryBonds.empty())
    {
        std::fprintf(fp, "; none\n");
    }

    std::fprintf(fp, "\n[ removed_bonded_terms ]\n");
    std::fprintf(fp, "; force-field terms removed from the topology, per interaction type.\n");
    std::fprintf(fp, "; class: boundary = QM and MM atoms, all-QM = described by the QM calculation,\n");
    std::fprintf(fp, ";        MM-only = no QM atom\n");
    for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        for (const char* cls : { "boundary", "all-QM", "MM-only" })
        {
            std::vector<const QmmmTopologyReport::Term*> terms;
            for (const auto& term : report.removedBonded)
            {
                if (term.ftype != ftype)
                {
                    continue;
                }
                int nQm = 0;
                for (const auto& a : term.atoms)
                {
                    nQm += a.isQm ? 1 : 0;
                }
                const char* c = nQm == static_cast<int>(term.atoms.size())
                                        ? "all-QM"
                                        : (nQm > 0 ? "boundary" : "MM-only");
                if (std::strcmp(c, cls) == 0)
                {
                    terms.push_back(&term);
                }
            }
            if (terms.empty())
            {
                continue;
            }
            std::fprintf(fp, "; %s, %s: %zu\n", interaction_function[ftype].longname, cls, terms.size());
            for (const auto* term : terms)
            {
                std::fprintf(fp, "  %-8s%s\n", cls, termText(*term).c_str());
            }
        }
    }
    if (report.removedBonded.empty())
    {
        std::fprintf(fp, "; none\n");
    }

    std::fprintf(fp, "\n[ kept_restraints ]\n");
    std::fprintf(fp, "; restraints with QM atoms: never removed, they are not part of the force field\n");
    for (const auto& term : report.keptRestraints)
    {
        std::fprintf(fp, "  %-12s%s\n", interaction_function[term.ftype].name, termText(term).c_str());
    }
    if (report.keptRestraints.empty())
    {
        std::fprintf(fp, "; none\n");
    }

    std::fprintf(fp, "\n[ qm_qm_bonds_converted_to_connections ]\n");
    std::fprintf(fp, "; chemical bonds between two QM atoms: no force any more, but still used for exclusions\n");
    for (const auto& term : report.convertedBonds)
    {
        std::fprintf(fp, "  %-12s%s\n", interaction_function[term.ftype].name, termText(term).c_str());
    }

    for (const bool boundary : { true, false })
    {
        std::vector<const QmmmTopologyReport::Term*> pairs;
        for (const auto& term : report.removedPairs)
        {
            const bool isBoundary = !(term.atoms[0].isQm && term.atoms[1].isQm);
            if (isBoundary == boundary)
            {
                pairs.push_back(&term);
            }
        }
        std::fprintf(fp, "\n[ removed_lj14_pairs_%s ]\n", boundary ? "qm_boundary_mm" : "qm_qm");
        std::fprintf(fp, "; pair interactions (LJ-14 with its Coulomb-14) removed from the topology: %zu\n",
                     pairs.size());
        if (boundary)
        {
            std::fprintf(fp, "; a QM atom with a boundary MM atom -- counted in the QM calculation via the link atom\n");
        }
        for (const auto* term : pairs)
        {
            std::fprintf(fp, "  %-8s%s\n", interaction_function[term->ftype].name, termText(*term).c_str());
        }
    }

    for (const bool boundary : { true, false })
    {
        const auto& list = boundary ? report.qmMmExclusions : report.qmQmExclusions;
        int         nNew = 0;
        for (const auto& e : list)
        {
            nNew += e.presentBefore ? 0 : 1;
        }
        std::fprintf(fp, "\n[ lj_exclusions_%s ]\n", boundary ? "qm_boundary_mm" : "qm_qm");
        std::fprintf(fp, "; nonbonded exclusions generated for QM/MM: %zu pairs, %d of them new, %zu already\n",
                     list.size(), nNew, list.size() - nNew);
        std::fprintf(fp, "; excluded by the force field (nrexcl). An excluded pair has neither LJ nor Coulomb;\n");
        std::fprintf(fp, "; the Coulomb of a QM atom is zero anyway, so what is removed here is the LJ.\n");
        std::fprintf(fp, "; %-26s %-26s %-8s %s\n", "QM atom", boundary ? "boundary MM atom" : "QM atom",
                     "bonds", "status");
        for (const auto& e : list)
        {
            std::fprintf(fp, "%s %s  %-8s %s\n", atomText(e.qm).c_str(), atomText(e.other).c_str(),
                         distanceText(e.bondDistance).c_str(),
                         e.presentBefore ? "already excluded by the force field" : "NEW: LJ removed by QM/MM");
        }
        if (boundary && !report.excludeBoundaryLJ)
        {
            std::fprintf(fp, "; none: the LJ of the boundary MM atoms follows the force field\n");
        }
    }

    for (const bool excl : { true, false })
    {
        const auto& list  = excl ? report.finalQmMmExclusions : report.keptQmMmPairs;
        int         nLink = 0;
        for (const auto& e : list)
        {
            nLink += e.presentBefore ? 1 : 0;
        }
        std::fprintf(fp, "\n[ %s ]\n", excl ? "final_lj_exclusions_qm_mm" : "final_lj14_pairs_qm_mm");
        std::fprintf(fp,
                     excl ? "; every QM--MM pair without LJ (nor Coulomb) in the tpr: %zu, %d of them with a link atom\n"
                          : "; every QM--MM pair with LJ-14 in the tpr: %zu, %d of them with a link atom\n",
                     list.size(), nLink);
        std::fprintf(fp, "; %-26s %-26s %-8s %s\n", "QM atom", "MM atom", "bonds", "");
        for (const auto& e : list)
        {
            std::fprintf(fp, "%s %s  %-8s %s\n", atomText(e.qm).c_str(), atomText(e.other).c_str(),
                         distanceText(e.bondDistance).c_str(), e.presentBefore ? "link atom" : "");
        }
    }
    std::fclose(fp);
}

void generate_qmexcl(gmx_mtop_t* sys, t_inputrec* ir, WarningHandler* wi, QmmmModeType qmmmMode, const gmx::MDLogger& logger)
{
    /* This routine expects molt->molt[m].ilist to be of size InteractionFunction::Count and ordered.
     */

    unsigned char*  grpnr;
    int             mol, nat_mol, nr_mol_with_qm_atoms = 0;
    gmx_molblock_t* molb;
    bool            bQMMM;
    int             index_offset = 0;
    // Counters of the removed force-field terms, summed over the molecule types
    //   that contain QM atoms, reported at the end of this routine.
    QmmmRemovedInteractions removed;
    // Per-atom record of the same changes, written to a separate file.
    QmmmTopologyReport report;

    // LJ between the QM and the MM atoms: by the exclusion rules of the force field
    //   (default), or with the boundary MM atoms excluded from every QM atom.
    bool excludeBoundaryLJ = false;
    if (qmmmMode != QmmmModeType::MiMiC)
    {
        const char* ljEnv = std::getenv("GMX_QMMM_LJ_SCHEME");
        if (ljEnv != nullptr && gmx_strcasecmp(ljEnv, "exclude") == 0)
        {
            excludeBoundaryLJ = true;
        }
        else if (ljEnv != nullptr && gmx_strcasecmp(ljEnv, "forcefield") != 0)
        {
            gmx_fatal(FARGS,
                      "Unknown value '%s' of the environment variable GMX_QMMM_LJ_SCHEME. "
                      "Use 'forcefield' (the default) or 'exclude'.",
                      ljEnv);
        }
        if (excludeBoundaryLJ && qmmmMode == QmmmModeType::Amber)
        {
            wi->addWarning(
                    "GMX_QMMM_LJ_SCHEME=exclude has no effect with GMX_QMMM_BONDED_SCHEME=amber: "
                    "the amber scheme keeps every interaction with an MM atom at the force-field "
                    "level, the LJ of the boundary MM atoms included.");
        }
        if (qmmmMode == QmmmModeType::Original)
        {
            GMX_LOG(logger.info)
                    .asParagraph()
                    .appendTextFormatted(
                            excludeBoundaryLJ
                                    ? "QM/MM: GMX_QMMM_LJ_SCHEME=exclude -- the LJ and LJ-14 "
                                      "interactions of every QM atom with the MM atoms bound to the "
                                      "QM region are excluded."
                                    : "QM/MM: the LJ between the QM and the MM atoms follows the "
                                      "exclusion rules of the force field (nrexcl and [ pairs ]); "
                                      "only the LJ within the QM region is excluded. To exclude also "
                                      "the LJ of the MM atoms bound to the QM region with every QM "
                                      "atom, set GMX_QMMM_LJ_SCHEME=exclude.");
        }
    }

    grpnr = sys->groups.groupNumbers[SimulationAtomGroupType::QuantumMechanics].data();

    for (size_t mb = 0; mb < sys->molblock.size(); mb++)
    {
        molb    = &sys->molblock[mb];
        nat_mol = sys->moltype[molb->type].atoms.nr;
        for (mol = 0; mol < molb->nmol; mol++)
        {
            bQMMM = FALSE;
            for (int i = 0; i < nat_mol; i++)
            {
                if ((grpnr ? grpnr[i] : 0) < (ir->opts.ngQM))
                {
                    bQMMM = TRUE;
                }
            }

            if (bQMMM)
            {
                nr_mol_with_qm_atoms++;
                if (molb->nmol > 1)
                {
                    /* We need to split this molblock */
                    if (mol > 0)
                    {
                        /* Split the molblock at this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        sys->molblock[mb].nmol = mol;
                        sys->molblock[mb + 1].nmol -= mol;
                        mb++;
                        molb = &sys->molblock[mb];
                    }
                    if (molb->nmol > 1)
                    {
                        /* Split the molblock after this molecule */
                        auto pos = sys->molblock.begin() + mb + 1;
                        sys->molblock.insert(pos, sys->molblock[mb]);
                        molb                   = &sys->molblock[mb];
                        sys->molblock[mb].nmol = 1;
                        sys->molblock[mb + 1].nmol -= 1;
                    }

                    /* Create a copy of a moltype for a molecule
                     * containing QM atoms and append it in the end of the list
                     */
                    std::vector<gmx_moltype_t> temp(sys->moltype.size());
                    for (size_t i = 0; i < sys->moltype.size(); ++i)
                    {
                        copy_moltype(&sys->moltype[i], &temp[i]);
                    }
                    sys->moltype.resize(sys->moltype.size() + 1);
                    for (size_t i = 0; i < temp.size(); ++i)
                    {
                        copy_moltype(&temp[i], &sys->moltype[i]);
                    }
                    copy_moltype(&sys->moltype[molb->type], &sys->moltype.back());
                    /* Copy the exclusions to a new array, since this is the only
                     * thing that needs to be modified for QMMM.
                     */
                    sys->moltype.back().excls = sys->moltype[molb->type].excls;
                    /* Set the molecule type for the QMMM molblock */
                    molb->type = sys->moltype.size() - 1;
                }
                generate_qmexcl_moltype(&sys->moltype[molb->type], grpnr, ir, qmmmMode, logger,
                                        &removed, index_offset, &report, excludeBoundaryLJ);
            }
            if (grpnr)
            {
                grpnr += nat_mol;
            }
            index_offset += nat_mol;
        }
    }
    if (qmmmMode != QmmmModeType::MiMiC && nr_mol_with_qm_atoms > 0)
    {
        reportQmmmRemovedInteractions(removed, report, qmmmMode, logger);
        if (qmmmReportsEnabled())
        {
            const char* reportFile = std::getenv("GMX_QMMM_TOPOLOGY_REPORT");
            if (reportFile == nullptr)
            {
                reportFile = "qmmm_topology_report.txt";
            }
            writeQmmmTopologyReport(report, qmmmMode, reportFile);
            GMX_LOG(logger.info)
                    .appendTextFormatted(
                            "QM/MM: every removed term, pair and exclusion is listed atom by atom in %s\n"
                            "       (file name set with GMX_QMMM_TOPOLOGY_REPORT, switched off with "
                            "GMX_QMMM_REPORTS=off).\n",
                            reportFile);
        }
    }
    if (qmmmMode != QmmmModeType::MiMiC && nr_mol_with_qm_atoms > 1)
    {
        /* generate a warning is there are QM atoms in different topologies.
         * In this case, it is not possible at this stage to mutualy exclude
         * the non-bonded interactions via the exclusions (AFAIK). Instead,
         * the user is advised to use the energy group exclusions in the mdp file
         */
        wi->addNote( "\nThe QM subsystem is divided over multiple topologies. "
                     "The mutual non-bonded interactions cannot be excluded. "
                     "There are two ways to achieve this:\n\n"
                     "1) merge the topologies, such that the atoms of the QM "
                     "subsystem are all present in one single topology file. "
                     "In this case this warning will dissappear\n\n"
                     "2) exclude the non-bonded interactions explicitly via the "
                     "energygrp-excl option in the mdp file. if this is the case "
                     "this warning may be ignored"
                     "\n\n");
    }
}
