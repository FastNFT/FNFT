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
* Lianne de Vries (TU Delft student) 2021.
* Sander Wahls (TU Delft) 2021.
*/

/**
 * @file fnft__manakov_scatter.h
 * @brief Slow forward scattering.
 * @ingroup manakov
 */

#ifndef FNFT__MANAKOV_SCATTER_H
#define FNFT__MANAKOV_SCATTER_H

#include "fnft__manakov_discretization.h"
#include <stdio.h>
#include <string.h> // for memcpy
#include "fnft__errwarn.h"
#include "fnft__misc.h"// for square_matrix_mult

/**
 * @brief Computes the scattering matrix
 * 
 * The function computes the scattering matrix with respect to \f$\lambda\f$. The bo method 
 * is based on the paper by Boffetta and Osborne (<a href="http://dx.doi.org/10.1016/0021-9991(92)90370-E">J. Comput. Physics 1992 </a>).
 * The CF4_2 method is based on the article by Chimmalgo, Prins and Wahls
 * (<a href="https://repository.tudelft.nl/islandora/object/uuid%3A5cbb52a1-6ec7-4368-bec3-a9d5137112b9"> IEEE 2018 </a>).
 * 
 * @param[in] D Number of samples
 * @param[in] q1 Array of length D, contains samples of the first element of the potential function \f$ q1_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q1_0, q1_1, \dots, q1_{D-1} \f$).
 * @param[in] q2 Array of length D, contains samples of the first element of the potential function \f$ q2_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q2_0, q2_1, \dots, q2_{D-1} \f$).
 * @param[in] eps_t timestep size
 * @param[in] K Number of values of \f$\lambda\f$.
 * @param[in] lambda Array with values of \f$\lambda\f$.
 * @param[in] kappa dispersion constant. kappa=-1 for normal dispersion, kappa=+1 for anomalous dispersion
 * @param[out] result Array of length 9*K containing S=[S11(lambda[0]) S12(lambda[0]) S13(lambda[0]) ... S33(lambda[K-1])]
 * upon exit where S = [S11, S12, S13; S21, S22, S23; S31, S32, S33] is the scattering matrix computes using the chosen discretization
 * @param[in] discretization Discretization to be used. Should be of type \link fnft_manakov_discretization_t \endlink.
 * Check \link fnft_manakov_discretization_t \endlink to see which discretizations correspond to slow methods
 */
FNFT_INT fnft__manakov_scatter_matrix(FNFT_UINT const D,
                        FNFT_COMPLEX const * const q1,
                        FNFT_COMPLEX const * const q2,
                        FNFT_REAL const eps_t,
                        FNFT_UINT const K,
                        FNFT_COMPLEX const * const lambda,
                        FNFT_INT const kappa,
                        FNFT_COMPLEX * const result,
                        fnft_manakov_discretization_t const discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_scatter_matrix(...) fnft__manakov_scatter_matrix(__VA_ARGS__)
#endif

#endif
