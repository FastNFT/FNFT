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
#include "fnft__akns_discretization.h"

/**
 * Returns the max degree of the polynomials in a single scattering
 * matrix or zero if the discretization is unknown.
 */
UINT fnft__nse_discretization_degree(nse_discretization_t
        nse_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = nse_discretization_to_akns_discretization(nse_discretization);
    return akns_discretization_degree(akns_discretization);
}


/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__nse_discretization_boundary_coeff(nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = nse_discretization_to_akns_discretization(nse_discretization);
    return akns_discretization_boundary_coeff(akns_discretization);
}

/**
 * This akns_discretization related to the given nse_discretization
 */
akns_discretization_t fnft__nse_discretization_to_akns_discretization(nse_discretization_t nse_discretization)
{

    switch (nse_discretization) {
        case nse_discretization_2SPLIT2_MODAL:
            return akns_discretization_2SPLIT2_MODAL;
        case nse_discretization_2SPLIT1A:
            return akns_discretization_2SPLIT1A;
        case nse_discretization_2SPLIT1B:
            return akns_discretization_2SPLIT1B;
        case nse_discretization_2SPLIT2A:
            return akns_discretization_2SPLIT2A;
        case nse_discretization_2SPLIT2B:
            return akns_discretization_2SPLIT2B;
        case nse_discretization_2SPLIT2S:
            return akns_discretization_2SPLIT2S;
        case nse_discretization_2SPLIT3S:
            return akns_discretization_2SPLIT3S;
        case nse_discretization_2SPLIT4B:
            return akns_discretization_2SPLIT4B;
        case nse_discretization_2SPLIT3A:
            return akns_discretization_2SPLIT3A;
        case nse_discretization_2SPLIT3B:
            return akns_discretization_2SPLIT3B;
        case nse_discretization_2SPLIT4A:
            return akns_discretization_2SPLIT4A;
        case nse_discretization_2SPLIT6B:
            return akns_discretization_2SPLIT6B;
        case nse_discretization_2SPLIT6A:
            return akns_discretization_2SPLIT6A;
        case nse_discretization_2SPLIT8B:
            return akns_discretization_2SPLIT8B;
        case nse_discretization_2SPLIT5A:
            return akns_discretization_2SPLIT5A;
        case nse_discretization_2SPLIT5B:
            return akns_discretization_2SPLIT5B;
        case nse_discretization_2SPLIT8A:
            return akns_discretization_2SPLIT8A;
        case nse_discretization_2SPLIT7A:
            return akns_discretization_2SPLIT7A;
        case nse_discretization_2SPLIT7B:
            return akns_discretization_2SPLIT7B;
        case nse_discretization_BO:
            return akns_discretization_BO;

            
            
        default: // Unknown discretization
            return NAN;
    }
    
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 */
COMPLEX fnft__nse_lambda_to_z(const COMPLEX lambda, const REAL eps_t, nse_discretization_t
        nse_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = nse_discretization_to_akns_discretization(nse_discretization);
    return akns_lambda_to_z(lambda,eps_t,akns_discretization);
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 */
COMPLEX fnft__nse_z_to_lambda(const COMPLEX z, const REAL eps_t, nse_discretization_t
        nse_discretization)
{
    akns_discretization_t akns_discretization;
    akns_discretization = nse_discretization_to_akns_discretization(nse_discretization);
    return akns_z_to_lambda(z,eps_t,akns_discretization);     
}
