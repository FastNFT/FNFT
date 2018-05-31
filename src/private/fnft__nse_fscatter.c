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
 * Shrinivas Chimmalgi (TU Delft) 2018.
 *
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__nse_fscatter.h"
#include "fnft__nse_discretization.h"
#include "fnft__discretization_t.h"
#include "fnft__misc.h"

/**
 * Returns the length of array to be allocated based on the number
 * of samples and discretization.
 */

UINT nse_fscatter_numel(UINT D, nse_discretization_t discretization)
{

    const UINT deg = nse_discretization_degree(discretization);
    if (deg == 0)
        return 0; // unknown discretization
    else
        return poly_fmult2x2_numel(deg, D);
}

INT nse_fscatter(const UINT D, COMPLEX const * const q,
        const REAL eps_t, const INT kappa,
        COMPLEX * const result, UINT * const deg_ptr,
        INT * const W_ptr, nse_discretization_t discretization)
{
    INT ret_code;
    UINT i;
    fnft__discretization_t akns_discretization;
    COMPLEX *r = NULL;
    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(eps_t);
    if (abs(kappa) != 1)
        return E_INVALID_ARGUMENT(kappa);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (deg_ptr == NULL)
        return E_INVALID_ARGUMENT(deg_ptr);
    
    switch (discretization) {
        case nse_discretization_2SPLIT2_MODAL:
            akns_discretization = discretization_2SPLIT2_MODAL;
            break;
        case nse_discretization_2SPLIT1A:
            akns_discretization = discretization_2SPLIT1A;
            break;
        case nse_discretization_2SPLIT1B:
            akns_discretization = discretization_2SPLIT1B;
            break;
        case nse_discretization_2SPLIT2A:
            akns_discretization = discretization_2SPLIT2A;
            break;
        case nse_discretization_2SPLIT2B:
            akns_discretization = discretization_2SPLIT2B;
            break;
        case nse_discretization_2SPLIT2S:
            akns_discretization = discretization_2SPLIT2S;
            break;
        case nse_discretization_2SPLIT3S:
            akns_discretization = discretization_2SPLIT3S;
            break;
        case nse_discretization_2SPLIT4B:
            akns_discretization = discretization_2SPLIT4B;
            break;
        case nse_discretization_2SPLIT3A:
            akns_discretization = discretization_2SPLIT3A;
            break;
        case nse_discretization_2SPLIT3B:
            akns_discretization = discretization_2SPLIT3B;
            break;
        case nse_discretization_2SPLIT4A:
            akns_discretization = discretization_2SPLIT4A;
            break;
        case nse_discretization_2SPLIT6B:
            akns_discretization = discretization_2SPLIT6B;
            break;
        case nse_discretization_2SPLIT6A:
            akns_discretization = discretization_2SPLIT6A;
            break;
        case nse_discretization_2SPLIT8B:
            akns_discretization = discretization_2SPLIT8B;
            break;
        case nse_discretization_2SPLIT5A:
            akns_discretization = discretization_2SPLIT5A;
            break;
        case nse_discretization_2SPLIT5B:
            akns_discretization = discretization_2SPLIT5B;
            break;
        case nse_discretization_2SPLIT8A:
            akns_discretization = discretization_2SPLIT8A;
            break;
        case nse_discretization_2SPLIT7A:
            akns_discretization = discretization_2SPLIT7A;
            break;
        case nse_discretization_2SPLIT7B:
            akns_discretization = discretization_2SPLIT7B;
            break;
            
            
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(discretization);
    }
    
    r = malloc(D*sizeof(COMPLEX));
    if (r == NULL) {
        ret_code = E_NOMEM;
    }
    
    if (kappa == 1){
        for (i = 0; i < D; i++)
            r[i] = -CONJ(q[i]);}
    else{
        for (i = 0; i < D; i++)
            r[i] = CONJ(q[i]);}
    
    ret_code = akns_fscatter(D, q, r, eps_t, result, deg_ptr, W_ptr, akns_discretization);
    return ret_code;
    
    
}
