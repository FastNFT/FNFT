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
*/

/**
 * @file fnft__nse_discretization.h
 * @brief Properties of the discretizations for the nonlinear Schroedinger equation.
 * @ingroup nse
 */

#ifndef FNFT__NSE_DISCRETIZATION_H
#define FNFT__NSE_DISCRETIZATION_H

#include "fnft_nse_discretization_t.h"

/**
 * @brief Returns the degree of the discretization scheme.
 * \n
 * This routine returns the degree of the discretization scheme of type
 * \link fnft_nse_discretization_t \endlink. \n
 * Returns 0 for discretizations not supported by \link fnft__nse_fscatter \endlink.
 * @ingroup nse
 */
FNFT_UINT fnft__nse_discretization_degree(fnft_nse_discretization_t
        discretization);

/**
 * @brief Returns the mapping coefficient based on discretization.
 * \n
 * This routine returns the mapping coefficient map_coeff based on the discretization of type 
 * \link fnft_nse_discretization_t \endlink. Then \f$ z=e^{map\_coeff.j.xi.eps_t} \f$.\n
 * Returns NAN for discretizations not supported by \link fnft__nse_fscatter \endlink.
 * @ingroup nse
 */
FNFT_REAL fnft__nse_discretization_mapping_coeff(fnft_nse_discretization_t discretization);

/**
 * @brief Returns the boundary coefficient based on discretization.
 * \n
 * This routine returns the boundary coefficient boundary_coeff based on the discretization of type 
 * \link fnft_nse_discretization_t \endlink. Then T[end] = T[1] + eps_t * boundary_coeff.
 * @ingroup nse
 */
FNFT_REAL fnft__nse_discretization_boundary_coeff(fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_discretization_degree(...) fnft__nse_discretization_degree(__VA_ARGS__)
#define nse_discretization_mapping_coeff(...) fnft__nse_discretization_mapping_coeff(__VA_ARGS__)
#define nse_discretization_boundary_coeff(...) fnft__nse_discretization_boundary_coeff(__VA_ARGS__)
#endif

#endif
