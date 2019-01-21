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
* Shrinivas Chimmalgi (TU Delft) 2017-2019.
* Peter J. Prins (TU Delft) 2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__akns_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 */
UINT fnft__akns_discretization_degree(akns_discretization_t
        discretization)
{
    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A: // With change of base trick
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT2_MODAL:
            return 1;
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_4SPLIT4B:
            return 2;
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
            return 3;
        case akns_discretization_2SPLIT4A:
        case akns_discretization_4SPLIT4A:
            return 4;
        case akns_discretization_2SPLIT6B:
            return 6;
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT8B:
            return 12;
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
            return 15;
        case akns_discretization_2SPLIT8A:
            return 24;
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
            return 105;
            
        default: // Unknown discretization
            return 0;
    }
}

/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__akns_discretization_boundary_coeff(akns_discretization_t discretization)
{

    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
            return 0.5;
            
        default: // Unknown discretization
            return NAN;
    }
}

/**
 * This routine returns the scaling for effective number of samples based on the discretization.
 */
UINT fnft__akns_discretization_D_scale(akns_discretization_t discretization)
{

    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
        case akns_discretization_BO:
            return 1;
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_CF4_2:
            return 2;
            
        default: // Unknown discretization
            return 0;
    }
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 */
INT fnft__akns_lambda_to_z(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, akns_discretization_t discretization)
{
    REAL degree1step;
    UINT i, D_scale;
    D_scale = akns_discretization_D_scale(discretization);
    if (D_scale == 0)
         return E_INVALID_ARGUMENT(discretization);
    degree1step = akns_discretization_degree(discretization);
    if (degree1step == 0)
         return E_INVALID_ARGUMENT(discretization);
    degree1step = degree1step * D_scale;
    for (i = 0; i < n; i++)
        vals[i] = CEXP(2*I*vals[i]*eps_t/degree1step);
    return SUCCESS;
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 */
INT fnft__akns_z_to_lambda(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, akns_discretization_t discretization)
{
    REAL degree1step;
    UINT i, D_scale;
    D_scale = akns_discretization_D_scale(discretization);
    if (D_scale == 0)
         return E_INVALID_ARGUMENT(discretization);
    degree1step = akns_discretization_degree(discretization);
    if (degree1step == 0)
         return E_INVALID_ARGUMENT(discretization);
    degree1step = degree1step * D_scale;
    for (i = 0; i < n; i++)
        vals[i] = CLOG(vals[i])/(2*I*eps_t/degree1step);
    return SUCCESS;    
}
