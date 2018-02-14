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
* Sander Wahls (TU Delft) 2017.
* Shrinivas Chimmalgi (TU Delft) 2017.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_discretization.h"

/**
 * Returns the max degree of the polynomials in a single scattering
 * matrix or zero if the discretization is unknown.
 */
UINT fnft__nse_discretization_degree(nse_discretization_t
        discretization)
{
    switch (discretization) {
        case nse_discretization_2SPLIT2_MODAL:
        case nse_discretization_2SPLIT2A:
            return 1;
        case nse_discretization_2SPLIT4A:
            return 4;
        case nse_discretization_2SPLIT4B:
            return 2;
            
        default:  // Unknown discretization
            return 0;
    }
}

/**
 * Returns the mapping coefficient map_coeff in the mapping
 * z = exp(map_coeff*I*\lam*eps_t) based on the discretization type.
 * Defined only for discretizations supported by nse_fscatter.c
 */
REAL fnft__nse_discretization_mapping_coeff(nse_discretization_t discretization)
{

    switch (discretization) {
        case nse_discretization_2SPLIT2_MODAL:
        case nse_discretization_2SPLIT2A:
            return 2.0;
        case nse_discretization_2SPLIT4A:
            return -0.5;
        case nse_discretization_2SPLIT4B:
            return -1.0;
            
        default: // Unknown discretization
            return NAN;
            
    }
}
/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__nse_discretization_boundary_coeff(nse_discretization_t discretization)
{

    switch (discretization) {
        case nse_discretization_2SPLIT2_MODAL:
        case nse_discretization_2SPLIT2A:
            return 0.0; // TODO: this value should be 0.5 for every staircase approximation that uses midpoint values.
        case nse_discretization_2SPLIT4A:
        case nse_discretization_2SPLIT4B:
            return 0.5;
            
        default: // Unknown discretization
            return NAN;
    }
}
