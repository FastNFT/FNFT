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
* Shrinivas Chimmalgi (TU Delft) 2017-2020.
*/

/**
 * @file fnft__akns_scatter.h
 * @brief Slow forward scattering.
 * @ingroup akns
 */

#ifndef FNFT__AKNS_SCATTER_H
#define FNFT__AKNS_SCATTER_H

#include "fnft__akns_discretization.h"
#include <stdio.h>
#include <string.h> // for memcpy
#include "fnft__errwarn.h"
#include "fnft__misc.h"// for square_matrix_mult


/**
 * @brief Computes the scattering matrix and its derivative.
 * 
 * The function computes the scattering matrix and the derivative of the scattering matrix with 
 * respect to \f$\lambda\f$. The function performs slow direct scattering and is primarily based on the reference 
 * Boffetta and Osborne 
 * (<a href="http://dx.doi.org/10.1016/0021-9991(92)90370-E">J. Comput. Physics 1992 </a>).
 * 
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q_0, q_1, \dots, q_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in,out] r Array of length D, contains samples \f$ r_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ r_0, r_1, \dots, r_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in] eps_t Step-size, eps_t \f$= (T[1]-T[0])/(D-1) \f$.
 * @param[in] K Number of values of \f$\lambda\f$.
 * @param[in] lambda Array of length K, contains the values of \f$\lambda\f$.
 * @param[out] result Array of length 8*K or 4*K,  If derivative_flag=0 returns 
 * [S11 S12 S21 S22] in result where S = [S11, S12; S21, S22] is the 
 * scattering matrix computed using the chosen discretization.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in 
 * result where S11' is the derivative of S11 w.r.t to \f$\lambda\f$.
 * Should be preallocated with size 4*K or 8*K accordingly.
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft__akns_discretization_t \endlink. Not all akns_discretization_t 
 * discretizations are supported. Check \link fnft__akns_discretization_t \endlink 
 * for list of supported types.
 * @param[in] derivative_flag Should be set to either 0 or 1. If set to 1
 * the derivatives [S11' S12' S21' S22'] are calculated. result should be
 * preallocated with size 8*K if flag is set to 1.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup akns
 */
FNFT_INT fnft__akns_scatter_matrix(const FNFT_UINT D, FNFT_COMPLEX const * const q, 
    FNFT_COMPLEX const * const r, const FNFT_REAL eps_t,
    const FNFT_UINT K, FNFT_COMPLEX const * const lambda,
    FNFT_COMPLEX * const result, fnft__akns_discretization_t discretization,
    const FNFT_UINT derivative_flag);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_scatter_matrix(...) fnft__akns_scatter_matrix(__VA_ARGS__)
#endif

#endif
