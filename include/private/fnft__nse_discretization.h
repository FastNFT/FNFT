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
 * @brief Properties of the discretizations for the nonlinear Schroedinger
 * equation.
 *
 * @file fnft__nse_discretization.h
 * @ingroup nse
 */

#ifndef FNFT__NSE_DISCRETIZATION_H
#define FNFT__NSE_DISCRETIZATION_H

#include "fnft_nse_discretization_t.h"

/**
 * @brief This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 *
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @returns polynomial degree, or 0 for discretizations not supported by \link fnft__nse_fscatter \endlink.
 *
 * @ingroup nse
 */
FNFT_UINT fnft__nse_discretization_degree(fnft_nse_discretization_t
        discretization);

/**
 * @brief Returns the mapping coefficient based on discretization.
 *
 * This routine returns the mapping coefficient map_coeff based on the discretization of type 
 * \link fnft_nse_discretization_t \endlink. Then \f$ z=e^{map\_coeff.j.xi.eps_t} \f$.\n
 * Returns NAN for discretizations not supported by \link fnft__nse_fscatter \endlink.
 *
 * @ingroup nse
 */
FNFT_REAL fnft__nse_discretization_mapping_coeff(fnft_nse_discretization_t discretization);

/**
 * @brief This routine returns the boundary coefficient based on the
 * discretization.
 *
 * The boundary coefficient is the fraction of the step size that a discretized
 * potential extends beyond the last sample. This routine returns this value
 * based on the discretization of type \link fnft_nse_discretization_t \endlink.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @returns the boundary coefficient, or NAN for discretizations not supported
 * by \link fnft__nse_fscatter \endlink.
 *
 * @ingroup nse
 */
FNFT_REAL fnft__nse_discretization_boundary_coeff(fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_discretization_degree(...) fnft__nse_discretization_degree(__VA_ARGS__)
#define nse_discretization_mapping_coeff(...) fnft__nse_discretization_mapping_coeff(__VA_ARGS__)
#define nse_discretization_boundary_coeff(...) fnft__nse_discretization_boundary_coeff(__VA_ARGS__)
#endif

#endif
