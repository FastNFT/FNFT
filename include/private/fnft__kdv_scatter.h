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
* Shrinivas Chimmalgi (TU Delft) 2017.
*/

/**
 * @file fnft__kdv_scatter.h
 * @brief Slow forward scattering.
 * @ingroup kdv
 */

#ifndef FNFT__KDV_SCATTER_H
#define FNFT__KDV_SCATTER_H

#include "fnft__kdv_discretization.h"
#include "fnft__akns_scatter.h"


/**
 * @brief Computes the scattering matrix and its derivative.
 * 
 * The function computes the scattering matrix and the derivative of the scattering matrix with 
 * respect to \f$\lambda\f$. The function performs slow direct scattering and is primarily based on the reference 
 * Boffetta and Osborne 
 * (<a href="http://dx.doi.org/10.1016/0021-9991(92)90370-E">J. Comput. Physics 1992 </a>).
 * 
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] eps_t Step-size, eps_t \f$= (T[1]-T[0])/(D-1) \f$.
 * @param[in] K Number of values of \f$\lambda\f$.
 * @param[in] lambda Array of length K, contains the values of \f$\lambda\f$.
 * @param[out] result Array of length 8*K, contains the values [S11 S12 S21 S22 S11' S12' S21' S22'] 
 * where S = [S11, S12; S21, S22] is the scattering matrix. 
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft_kdv_discretization_t \endlink. Not all kdv_discretization_t discretizations are supported.
 * Check \link fnft_kdv_discretization_t \endlink for list of supported types.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup kdv
 */
FNFT_INT fnft__kdv_scatter_matrix(const FNFT_UINT D, FNFT_COMPLEX const * const q,
    const FNFT_REAL eps_t, const FNFT_UINT K, 
    FNFT_COMPLEX const * const lambda,
    FNFT_COMPLEX * const result, fnft_kdv_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_scatter_matrix(...) fnft__kdv_scatter_matrix(__VA_ARGS__)
#endif

#endif
