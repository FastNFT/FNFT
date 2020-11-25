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
* Sander Wahls (TU Delft) 2017-2018.
* Peter J Prins (TU Delft) 2020.
*/

/**
 * \file fnft__kdv_fscatter.h
 * @brief Computes the polynomial approximation of the combined scattering
 * matrix.
 * @ingroup kdv
 */

#ifndef FNFT__KDV_FSCATTER_H
#define FNFT__KDV_FSCATTER_H


#include "fnft__kdv_discretization.h"
#include "fnft__akns_fscatter.h"

/**
 * @brief Returns the length of transfer_matrix to be allocated based on the number
 * of samples and discretization.
 *
 * @ingroup kiv
 * This routine returns the length 4*D*(kdv_discretization_degree(discretization) + 1)
 * to be allocated based on the number
 * of samples and discretization of type discretization.
 * @param[in] D Number of samples.
 * @param[in] discretization Type of discretization from \link fnft_kdv_discretization_t \endlink.
 * @returns Returns the length to be allocated. Returns 0 for unknown discretizations.
 */
FNFT_UINT fnft__kdv_fscatter_numel(FNFT_UINT D,
    fnft_kdv_discretization_t discretization);

/**
 * @brief Fast computation of polynomial approximation of the combined scattering
 * matrix.
 *
 * This routine computes the polynomial approximation of the combined scattering
 * matrix by multipying together individual scattering matrices.\n
 * Individual scattering matrices depend on the chosen discretization.
 * \n
 * The main reference is Wahls and Poor
 * (<a href="http://dx.doi.org/10.1109/ICASSP.2013.6638772">Proc. ICASSP 2013 </a>).
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] eps_t Step-size, eps_t \f$= (T[1]-T[0])/(D-1) \f$.
 * @param[in] kappa unused.
 * @param[out] result array of length `kdv_fscatter_numel(D,discretization)`,
 * will contain the combined scattering matrix. Result needs to be pre-allocated
 * with `malloc(kdv_fscatter_numel(D,discretization)*sizeof(COMPLEX))`.
 * @param[out] deg_ptr Pointer to variable containing degree of the discretization.
 * Determined based on discretization by \link fnft__kdv_discretization_degree \endlink.
 * @param[in] W_ptr Pointer to normalization flag \link fnft_kdvv_opts_t::normalization_flag \endlink.
 * @param[in] discretization The type of discretization to be used. Should be of type
 * \link fnft_kdv_discretization_t \endlink. Not all fnft_kdv_discretization_t discretizations are supported.
 * Check \link fnft_kdv_discretization_t \endlink for list of supported types.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup kdv
 */
FNFT_INT fnft__kdv_fscatter(const FNFT_UINT D,
                            FNFT_COMPLEX const * const q,
                            FNFT_COMPLEX const * const r,
                            const FNFT_REAL eps_t,
                            const FNFT_INT kappa,
                            FNFT_COMPLEX * const result,
                            FNFT_UINT * const deg_ptr,
                            FNFT_INT * const W_ptr,
                            fnft_kdv_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_fscatter_numel(...) fnft__kdv_fscatter_numel(__VA_ARGS__)
#define kdv_fscatter(...) fnft__kdv_fscatter(__VA_ARGS__)
#endif

#endif
