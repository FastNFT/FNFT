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
* Peter J Prins (TU Delft) 2018
*/

/**
 * \file fnft__kdv_fscatter.h
 * @brief Computes the polynomial approximation of the combined scattering 
 * matrix.
 * @ingroup kdv
 */

#ifndef FNFT__KDV_FSCATTER_H
#define FNFT__KDV_FSCATTER_H

#include "fnft_kdv_discretization_t.h"

/**
 * Calculate M = [cos(eps_t*sqrt(q)),         q * eps * sinc(eps_t*sqrt(q) ;
 *                -eps_t*sinc(eps_t*sqrt(q)), cos(eps_t*sqrt(q))           ];
 */
FNFT_INT fnft__kdv_fscatter_zero_freq_scatter_matrix(COMPLEX *M,
                                               const REAL eps_t, const REAL q);

/**
 * @brief Returns the length of vector to be allocated based on the number
 * of samples and discretization.
 * \n
 * This routine returns the length of vector to be allocated based on the number
 * of samples and discretization of type 
 * \link fnft_kdv_discretization_t \endlink. 
 * @ingroup kdv
 */
FNFT_UINT fnft__kdv_fscatter_length(FNFT_UINT D,
        fnft_kdv_discretization_t discretization);

/**
 * @brief Fast computation of polynomial approximation of the combined scattering 
 * matrix.
 * @ingroup kdv
 */
//FNFT_INT fnft__kdv_fscatter(const FNFT_UINT D, FNFT_COMPLEX *Q,
//    const FNFT_REAL eps_t, FNFT_COMPLEX *result, FNFT_UINT *deg_ptr,
//    fnft_kdv_discretization_t discretization);

FNFT_INT fnft__kdv_fscatter(const FNFT_UINT D, FNFT_COMPLEX const * const q,
                 const FNFT_REAL eps_t, FNFT_COMPLEX * const result, FNFT_UINT * const deg_ptr,
                            fnft_kdv_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_fscatter_length(...) fnft__kdv_fscatter_length(__VA_ARGS__)
#define kdv_fscatter(...) fnft__kdv_fscatter(__VA_ARGS__)
#endif

#endif
