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
    INT ret_code;
    UINT degree1step = 0;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    degree1step = akns_discretization_degree(akns_discretization); 
    leave_fun:    
        return degree1step; 
}


/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__nse_discretization_boundary_coeff(nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization;
    REAL bnd_coeff = NAN;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);    
    bnd_coeff = akns_discretization_boundary_coeff(akns_discretization);
    leave_fun:    
        return bnd_coeff;
}

/**
 * This akns_discretization related to the given nse_discretization
 */
INT fnft__nse_discretization_to_akns_discretization(nse_discretization_t nse_discretization, 
        akns_discretization_t * const akns_discretization)
{

    switch (nse_discretization) {
        case nse_discretization_2SPLIT2_MODAL:
            *akns_discretization = akns_discretization_2SPLIT2_MODAL;
            break;
        case nse_discretization_2SPLIT1A:
            *akns_discretization = akns_discretization_2SPLIT1A;
            break;
        case nse_discretization_2SPLIT1B:
            *akns_discretization = akns_discretization_2SPLIT1B;
            break;
        case nse_discretization_2SPLIT2A:
            *akns_discretization = akns_discretization_2SPLIT2A;
            break;
        case nse_discretization_2SPLIT2B:
            *akns_discretization = akns_discretization_2SPLIT2B;
            break;
        case nse_discretization_2SPLIT2S:
            *akns_discretization = akns_discretization_2SPLIT2S;
            break;
        case nse_discretization_2SPLIT3S:
            *akns_discretization = akns_discretization_2SPLIT3S;
            break;
        case nse_discretization_2SPLIT4B:
            *akns_discretization = akns_discretization_2SPLIT4B;
            break;
        case nse_discretization_2SPLIT3A:
            *akns_discretization = akns_discretization_2SPLIT3A;
            break;
        case nse_discretization_2SPLIT3B:
            *akns_discretization = akns_discretization_2SPLIT3B;
            break;
        case nse_discretization_2SPLIT4A:
            *akns_discretization = akns_discretization_2SPLIT4A;
            break;
        case nse_discretization_2SPLIT6B:
            *akns_discretization = akns_discretization_2SPLIT6B;
            break;
        case nse_discretization_2SPLIT6A:
            *akns_discretization = akns_discretization_2SPLIT6A;
            break;
        case nse_discretization_2SPLIT8B:
            *akns_discretization = akns_discretization_2SPLIT8B;
            break;
        case nse_discretization_2SPLIT5A:
            *akns_discretization = akns_discretization_2SPLIT5A;
            break;
        case nse_discretization_2SPLIT5B:
            *akns_discretization = akns_discretization_2SPLIT5B;
            break;
        case nse_discretization_2SPLIT8A:
            *akns_discretization = akns_discretization_2SPLIT8A;
            break;
        case nse_discretization_2SPLIT7A:
            *akns_discretization = akns_discretization_2SPLIT7A;
            break;
        case nse_discretization_2SPLIT7B:
            *akns_discretization = akns_discretization_2SPLIT7B;
            break;
        case nse_discretization_BO:
            *akns_discretization = akns_discretization_BO;
            break;            
            
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(nse_discretization);
    }
    return SUCCESS;    
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 */
INT fnft__nse_lambda_to_z(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization;
    INT ret_code = SUCCESS;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__akns_lambda_to_z(n, eps_t, vals, akns_discretization);
    leave_fun:    
        return ret_code;
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 */
INT fnft__nse_z_to_lambda(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization;
    INT ret_code = SUCCESS;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__akns_z_to_lambda(n, eps_t, vals, akns_discretization);
    leave_fun:    
        return ret_code;   
}
