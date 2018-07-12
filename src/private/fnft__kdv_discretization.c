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
#include "fnft__akns_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 */
UINT fnft__kdv_discretization_degree(kdv_discretization_t
        kdv_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = kdv_discretization_to_akns_discretization(kdv_discretization);
    return akns_discretization_degree(akns_discretization);
}

/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__kdv_discretization_boundary_coeff(kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = kdv_discretization_to_akns_discretization(kdv_discretization);
    return akns_discretization_boundary_coeff(akns_discretization);
}


/**
 * This akns_discretization related to the given kdv_discretization
 */
akns_discretization_t fnft__kdv_discretization_to_akns_discretization(kdv_discretization_t kdv_discretization)
{

    switch (kdv_discretization) {
        case kdv_discretization_2SPLIT1A:
            return akns_discretization_2SPLIT1A;
        case kdv_discretization_2SPLIT1B:
            return akns_discretization_2SPLIT1B;
        case kdv_discretization_2SPLIT2A:
            return akns_discretization_2SPLIT2A;
        case kdv_discretization_2SPLIT2B:
            return akns_discretization_2SPLIT2B;
        case kdv_discretization_2SPLIT2S:
            return akns_discretization_2SPLIT2S;
        case kdv_discretization_2SPLIT3S:
            return akns_discretization_2SPLIT3S;
        case kdv_discretization_2SPLIT4B:
            return akns_discretization_2SPLIT4B;
        case kdv_discretization_2SPLIT3A:
            return akns_discretization_2SPLIT3A;
        case kdv_discretization_2SPLIT3B:
            return akns_discretization_2SPLIT3B;
        case kdv_discretization_2SPLIT4A:
            return akns_discretization_2SPLIT4A;
        case kdv_discretization_2SPLIT6B:
            return akns_discretization_2SPLIT6B;
        case kdv_discretization_2SPLIT6A:
            return akns_discretization_2SPLIT6A;
        case kdv_discretization_2SPLIT8B:
            return akns_discretization_2SPLIT8B;
        case kdv_discretization_2SPLIT5A:
            return akns_discretization_2SPLIT5A;
        case kdv_discretization_2SPLIT5B:
            return akns_discretization_2SPLIT5B;
        case kdv_discretization_2SPLIT8A:
            return akns_discretization_2SPLIT8A;
        case kdv_discretization_2SPLIT7A:
            return akns_discretization_2SPLIT7A;
        case kdv_discretization_2SPLIT7B:
            return akns_discretization_2SPLIT7B;
        case kdv_discretization_BO:
            return akns_discretization_BO;

            
            
        default: // Unknown discretization
            return NAN;
    }
    
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 */
COMPLEX fnft__kdv_lambda_to_z(const COMPLEX lambda, const REAL eps_t, kdv_discretization_t
        kdv_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = kdv_discretization_to_akns_discretization(kdv_discretization);
    return akns_lambda_to_z(lambda,eps_t,akns_discretization);
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 */
COMPLEX fnft__kdv_z_to_lambda(const COMPLEX z, const REAL eps_t, kdv_discretization_t
        kdv_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = kdv_discretization_to_akns_discretization(kdv_discretization);
    return akns_z_to_lambda(z,eps_t,akns_discretization);     
}
