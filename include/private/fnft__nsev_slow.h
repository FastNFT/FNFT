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
* Shrinivas Chimmalgi (TU Delft) 2020.
*/

/**
 * @file fnft__nsev_slow.h
 * @brief Auxiliary functions used in fnft_nsev and fnft_nsev_slow.
 * @ingroup nse
 */

#ifndef FNFT__NSEV_SLOW_H
#define FNFT__NSEV_SLOW_H

#include "fnft__nse_discretization.h"
#include "fnft__akns_scatter.h"
#include "fnft__errwarn.h"

/**
 * @brief Refines initial guesses of bound states using Newton's method
 * 
 * The function performs Newton's iterations to refine initial guesses of bound states
 * for the nonlinear Schroedinger equation.
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q_0, q_1, \dots, q_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in,out] r Array of length D, contains samples \f$ r_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ r_0, r_1, \dots, r_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization. 
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be T[0]<T[1].
 * @param[in] K Number of bound-states.
 * @param[in,out] bound_states Array of length K, on entry should contain 
 * the initial guesses for the bound-states \f$\lambda\f$. On return it 
 * contains the refined bound states.
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft_nse_discretization_t \endlink. Not all nse_discretization_t discretizations are supported.
 * Check \link fnft_nse_discretization_t \endlink for list of supported types.
 * @param[in] niter Positive integer. It is used as the upper limit
 * on the number of Newton interations that will be performed per bound state.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup nse
 */
FNFT_INT fnft__nse_refine_roots_newton(const FNFT_UINT D,
        FNFT_COMPLEX const * const q,
        FNFT_COMPLEX * r,
        FNFT_REAL const * const T,
        const FNFT_UINT K,
        FNFT_COMPLEX * bound_states,
        fnft_nse_discretization_t discretization,
        const FNFT_UINT niter);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_refine_roots_newton(...) fnft__nse_refine_roots_newton(__VA_ARGS__)
#endif

#endif
