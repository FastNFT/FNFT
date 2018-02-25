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
    // 2x2 matrix of degree+1 elements for each sample
    return 4*(nse_discretization_degree(discretization) + 1)*D;
}

/**
 * Returns the scattering matrix for a single step at frequency zero.
 */
/*
INT kdv_fscatter_zero_freq_scatter_matrix(COMPLEX *M,
                                                const REAL eps_t, const REAL q)
{
    REAL Delta = eps_t * CSQRT(q);
    M[0] = CCOS(Delta);                // M(1,1) = M(2,2)
    M[2] = -eps_t * misc_CSINC(Delta); // M(2,1)
    M[1] = -q * M[2];                  // M(1,2)
    return 0;
}
*/

INT nse_fscatter(const UINT D, COMPLEX const * const q,
    const REAL eps_t, const INT kappa,
    COMPLEX * const result, UINT * const deg_ptr,
    INT * const W_ptr, nse_discretization_t discretization)
{
    INT ret_code;
    fnft__discretization_opts_t opts = {
	.discretization = discretization_2SPLIT2A,
	.evolution_equation = evolution_equation_fnse
	};

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
    
    if (kappa == 1)
	opts.evolution_equation = evolution_equation_fnse;
    else
	opts.evolution_equation = evolution_equation_dnse;

    switch (discretization) {
        case nse_discretization_2SPLIT1A:
	   opts.discretization = discretization_2SPLIT1A;
	   break;
        case nse_discretization_2SPLIT1B:
	   opts.discretization = discretization_2SPLIT1B;
	   break;
        case nse_discretization_2SPLIT2A: 
	   opts.discretization = discretization_2SPLIT2A;
	   break;
        case nse_discretization_2SPLIT2B:
	   opts.discretization = discretization_2SPLIT2B;
	   break;
        case nse_discretization_2SPLIT2S:
	   opts.discretization = discretization_2SPLIT2S;
	   break;
        case nse_discretization_2SPLIT3S:
	   opts.discretization = discretization_2SPLIT3S;
	   break;
        case nse_discretization_2SPLIT4B:
	   opts.discretization = discretization_2SPLIT4B;
	   break;
        case nse_discretization_2SPLIT3A:
	   opts.discretization = discretization_2SPLIT3A;
	   break;
        case nse_discretization_2SPLIT3B:
	   opts.discretization = discretization_2SPLIT3B;
	   break;
        case nse_discretization_2SPLIT4A:
	   opts.discretization = discretization_2SPLIT4A;
	   break;
        case nse_discretization_2SPLIT6B:
	   opts.discretization = discretization_2SPLIT6B;
	   break;
        case nse_discretization_2SPLIT6A:
	   opts.discretization = discretization_2SPLIT6A;
	   break;
        case nse_discretization_2SPLIT8B:
	   opts.discretization = discretization_2SPLIT8B;
	   break;
        case nse_discretization_2SPLIT5A:
	   opts.discretization = discretization_2SPLIT1A;
	   break;
        case nse_discretization_2SPLIT5B:
	   opts.discretization = discretization_2SPLIT1A;
	   break;
        case nse_discretization_2SPLIT8A:
	   opts.discretization = discretization_2SPLIT8A;
	   break;
        case nse_discretization_2SPLIT7A:
	   opts.discretization = discretization_2SPLIT7A;
	   break;
        case nse_discretization_2SPLIT7B:
	   opts.discretization = discretization_2SPLIT7B;
	   break;
 
            
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(discretization);
    }

    ret_code = akns_fscatter(D, q, eps_t, result, deg_ptr, W_ptr, &opts);
    return ret_code;
    

}
