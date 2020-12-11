/*
 * This file is part of FNFT.
 *
 * FNFT is free software; you can redistribute it and/or
 * modify it under the terms of the version 2 of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * FNFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 * Sander Wahls (TU Delft) 2018, 2020.
 */

/**
 * \file fnft__kdv_finvscatter.h
 * @brief Recovers the signal from a scattering matrix.
 * @ingroup kdv
 */

#ifndef FNFT__KDV_FINVSCATTER_H
#define FNFT__KDV_FINVSCATTER_H

#include "fnft__kdv_discretization.h"

/**
 * @brief Recovers the samples that corresponding to a transfer matrix fast.
 *
 * @ingroup kdv
 * @param deg Degree of the polynomials in the transfer matrix.
 * @param [in] transfer_matrix A transfer matrix in the same format as used by
 *   \link fnft__kdv_fscatter \endlink.
 * @param [out] q Array with D=deg/base_deg entries in which the samples are
 *   stored, where base_deg is the output of
 *   \link fnft__kdv_discretization_degree \endlink.
 * @param[in] eps_t See \link fnft__kdv_fscatter \endlink.
 * @param[in] discretization See \link fnft__kdv_fscatter \endlink. Currently,
 *   only the 2SPLIT2_MODAL and 2SPLIT2A discretizations are supported.
 */
FNFT_INT fnft__kdv_finvscatter(
    const FNFT_UINT deg,
    FNFT_COMPLEX * const transfer_matrix,
    FNFT_COMPLEX * const q,
    const FNFT_REAL eps_t,
    const fnft_kdv_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_finvscatter(...) fnft__kdv_finvscatter(__VA_ARGS__)
#endif

#endif
