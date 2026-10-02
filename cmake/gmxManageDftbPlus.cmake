#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# Locate DFTB+ for the QM/MM interface (GMX_QMMM_PROGRAM=dftbplus).
#
# Any DFTB+ installation with the C API (WITH_API=ON) and its exported CMake
# package is accepted -- releases 21.x to 25.x as well as forks. The package is
# found with find_package(DftbPlus CONFIG) below GMX_QMMM_DFTBPLUS_LIB, which may
# be the install prefix of DFTB+ or its lib/ subdirectory.
#
# The C API changed between the versions, so it is not assumed but read from the
# installed header dftbplus.h:
#   GMX_QMMM_DFTBPLUS_ATOM_LIST      dftbp_process_input() takes a DftbPlusAtomList,
#                                    i.e. the QM atoms and species can be passed from
#                                    the GROMACS topology (some forks). Otherwise they
#                                    are taken from the Geometry block of dftb_in.hsd.
#   GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS  dftbp_get_atomic_shifts() exists (some forks),
#                                    needed for the output GMX_DFTB_ATOMIC_SHIFTS.
#
# On return the imported target DftbPlus::DftbPlus exists and
# GMX_DFTBPLUS_INCLUDE_DIR holds the directory of dftbplus.h.

macro(gmx_manage_dftbplus)
    set(GMX_QMMM_DFTBPLUS_LIB "" CACHE PATH
        "Install prefix of DFTB+ (the directory that contains lib/cmake/dftbplus); its lib/ subdirectory is accepted, too")

    # The exported package of DFTB+ depends on OpenMP::OpenMP_Fortran, which
    # FindOpenMP only provides when Fortran is an enabled language.
    enable_language(Fortran)

    set(_gmx_dftbp_hints "")
    if(GMX_QMMM_DFTBPLUS_LIB)
        get_filename_component(_gmx_dftbp_given "${GMX_QMMM_DFTBPLUS_LIB}" ABSOLUTE)
        get_filename_component(_gmx_dftbp_parent "${_gmx_dftbp_given}" DIRECTORY)
        list(APPEND _gmx_dftbp_hints "${_gmx_dftbp_given}" "${_gmx_dftbp_parent}")
    endif()
    find_package(DftbPlus CONFIG REQUIRED HINTS ${_gmx_dftbp_hints})

    # Not every version puts dftbplus.h into the include directories of its target
    get_filename_component(_gmx_dftbp_prefix "${DftbPlus_DIR}/../../.." ABSOLUTE)
    find_path(GMX_DFTBPLUS_INCLUDE_DIR dftbplus.h
              HINTS "${_gmx_dftbp_prefix}/include" "${_gmx_dftbp_prefix}/include/dftbplus"
              NO_DEFAULT_PATH)
    if(NOT GMX_DFTBPLUS_INCLUDE_DIR)
        message(FATAL_ERROR "DFTB+ was found in ${_gmx_dftbp_prefix}, but its C API header dftbplus.h "
                            "was not. DFTB+ has to be built with WITH_API=ON.")
    endif()
    mark_as_advanced(GMX_DFTBPLUS_INCLUDE_DIR)

    file(READ "${GMX_DFTBPLUS_INCLUDE_DIR}/dftbplus.h" _gmx_dftbp_header)
    if(_gmx_dftbp_header MATCHES "dftbp_process_input[ \t\r\n]*\\([^)]*DftbPlusAtomList")
        set(GMX_QMMM_DFTBPLUS_ATOM_LIST 1)
    else()
        set(GMX_QMMM_DFTBPLUS_ATOM_LIST 0)
    endif()
    if(_gmx_dftbp_header MATCHES "dftbp_get_atomic_shifts")
        set(GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS 1)
    else()
        set(GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS 0)
    endif()
    unset(_gmx_dftbp_header)

    if(GMX_QMMM_DFTBPLUS_ATOM_LIST)
        set(_gmx_dftbp_atoms "from the GROMACS topology")
    else()
        set(_gmx_dftbp_atoms "from the Geometry block of dftb_in.hsd")
    endif()
    if(GMX_QMMM_DFTBPLUS_ATOMIC_SHIFTS)
        set(_gmx_dftbp_shifts "available")
    else()
        set(_gmx_dftbp_shifts "not available")
    endif()
    message(STATUS "Found DFTB+ (C API ${DftbPlus_VERSION}) for QM/MM in ${_gmx_dftbp_prefix}: "
                   "QM atoms and species ${_gmx_dftbp_atoms}, atomic shifts ${_gmx_dftbp_shifts}")
endmacro()
