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
* Peter J. Prins (TU Delft) 2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdv_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 */
UINT fnft__kdv_discretization_degree(kdv_discretization_t
        discretization)
{
    switch (discretization) {
        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2A: // With change of base trick
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
            return 1;
        case kdv_discretization_2SPLIT4B:
            return 2;
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
            return 3;
        case kdv_discretization_2SPLIT4A:
            return 4;
        case kdv_discretization_2SPLIT6B:
            return 6;
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT8B:
            return 12;
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
            return 15;
        case kdv_discretization_2SPLIT8A:
            return 24;
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
            return 105;
            
        default: // Unknown discretization
            return 0;
    }
}

/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__kdv_discretization_boundary_coeff(kdv_discretization_t discretization)
{

    switch (discretization) {
        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2A:
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
        case kdv_discretization_2SPLIT4A:
        case kdv_discretization_2SPLIT4B:
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT6B:
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
        case kdv_discretization_2SPLIT8A:
        case kdv_discretization_2SPLIT8B:
            return 0.5;
            
        default: // Unknown discretization
            return NAN;
    }
}
