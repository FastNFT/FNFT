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
 * Sander Wahls (TU Delft) 2018.
 */

/**
 * \file fnft__nse_finvscatter.h
 * @brief Recovers the signal from a scattering matrix.
 * @ingroup nse
 */

#ifndef FNFT__NSE_FINVSCATTER_H
#define FNFT__NSE_FINVSCATTER_H

#include "fnft_nse_discretization_t.h"

/**
 * @brief Recovers the samples that corresponding to a transfer matrix fast.
 *
 * Recovers samples q[0], q[1], ..., q[D-1] that, when applied to \link
 * fnft__nse_fscatter \endlink, the transfer matrix provided by the user is
 * reconstructed. Note that the transfer matrix has to be (at least close to)
 * realizable, i.e., such samples have to exist. Also, the transfer matrix
 * has to be sufficiently nice to keep numerical errors tolerable, which
 * in my experience means that the q[n] it corresponds to are "smooth enough"
 * and such that |eps_t*q[n]|<<1 for all n.
 *
 * More information about the algorithm used here can be found in
 * <a href="http://dx.doi.org/10.1109/ISIT.2015.7282741">Wahls and Poor (Proc.
 * IEEE ISIT 2015)</a> and <a href="https://doi.org/10.1190/1.1441417">McClary
 * (Geophysics 48(10), 1983)</a>.
 *
 * @ingroup nse
 * @param deg Degree of the polynomials in the transfer matrix.
 * @param [in] transfer_matrix A transfer matrix in the same format as used by
 *   \link fnft__nse_fscatter \endlink.
 * @param [out] q Array with D=deg/base_deg entries in which the samples are
 *   stored, where base_deg is the output of
 *   \link fnft__nse_discretization_degree \endlink.
 * @param eps_t See \link fnft__nse_fscatter \endlink.
 * @param kappa See \link fnft__nse_fscatter \endlink.
 * @param discretization See \link fnft__nse_fscatter \endlink. Currently,
 *   only the 2SPLIT2_MODAL discretization is supported.
 */
FNFT_INT fnft__nse_finvscatter(
    const FNFT_UINT deg,
    FNFT_COMPLEX * const transfer_matrix,
    FNFT_COMPLEX * const q,
    const FNFT_REAL eps_t,
    const FNFT_INT kappa,
    const fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_finvscatter(...) fnft__nse_finvscatter(__VA_ARGS__)
#endif

#endif
