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
 * @file fnft__nse_scatter.h
 * @brief Slow forward scattering.
 * @ingroup nse
 */

#ifndef FNFT__NSE_SCATTER_H
#define FNFT__NSE_SCATTER_H

#include "fnft__nse_discretization.h"
#include "fnft__akns_scatter.h"

/**
 * @brief Computes \f$a(\lambda)\f$, \f$ a'(\lambda) = \frac{\partial a(\lambda)}{\partial \lambda}\f$
 * and \f$b(\lambda)\f$ for complex values \f$\lambda\f$ assuming that they are very close to the true 
 * bound-states.
 * 
 * The function performs slow direct scattering and is primarily based on the reference 
 * Boffetta and Osborne 
 * (<a href="http://dx.doi.org/10.1016/0021-9991(92)90370-E">J. Comput. Physics 1992 </a>).
 * A forward-backward scheme as mentioned by Aref in 
 * (<a href="https://arxiv.org/pdf/1605.06328.pdf"> Unpublished</a>)
 * is used to compute the norming constants \f$b(\lambda)\f$.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be T[0]<T[1].
 * @param[in] K Number of bound-states.
 * @param[in] bound_states Array of length K, contains the bound-states \f$\lambda\f$.
 * @param[out] a_vals Array of length K, contains the values of \f$a(\lambda)\f$.
 * @param[out] aprime_vals Array of length K, contains the values of
 * \f$ a'(\lambda) = \frac{\partial a(\lambda)}{\partial \lambda}\f$.
 * @param[out] b Array of length K, contains the values of \f$b(\lambda)\f$.
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft_nse_discretization_t \endlink. Not all nse_discretization_t discretizations are supported.
 * Check \link fnft_nse_discretization_t \endlink for list of supported types.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup nse
 */
FNFT_INT fnft__nse_scatter_bound_states(const FNFT_UINT D, FNFT_COMPLEX const *const q,
    FNFT_REAL const *const T, FNFT_UINT K,
    FNFT_COMPLEX *bound_states, FNFT_COMPLEX *a_vals,
    FNFT_COMPLEX *aprime_vals, FNFT_COMPLEX *b,
    fnft_nse_discretization_t discretization);

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
 * @param[in] kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one
 * @param[in] K Number of values of \f$\lambda\f$.
 * @param[in] lambda Array of length K, contains the values of \f$\lambda\f$.
 * @param[out] result Array of length 8*K, contains the values [S11 S12 S21 S22 S11' S12' S21' S22'] 
 * where S = [S11, S12; S21, S22] is the scattering matrix. 
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft_nse_discretization_t \endlink. Not all nse_discretization_t discretizations are supported.
 * Check \link fnft_nse_discretization_t \endlink for list of supported types.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup nse
 */
FNFT_INT fnft__nse_scatter_matrix(const FNFT_UINT D, FNFT_COMPLEX const * const q,
    const FNFT_REAL eps_t, const FNFT_INT kappa, const FNFT_UINT K, 
    FNFT_COMPLEX const * const lambda,
    FNFT_COMPLEX * const result, fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_scatter_bound_states(...) fnft__nse_scatter_bound_states(__VA_ARGS__)
#define nse_scatter_matrix(...) fnft__nse_scatter_matrix(__VA_ARGS__)
#endif

#endif
