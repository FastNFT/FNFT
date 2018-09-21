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
* Peter J Prins (TU Delft) 2017-2018.
* Shrinivas Chimmalgi (TU Delft) 2018.

*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__kdv_fscatter.h"
#include "fnft__kdv_discretization.h"
#include "fnft__akns_discretization.h"
#include "fnft__misc.h"

/**
 * Returns the length of array to be allocated based on the number
 * of samples and discretization.
 */
UINT kdv_fscatter_numel(UINT D, kdv_discretization_t discretization)
{
    const UINT deg = kdv_discretization_degree(discretization);
    if (deg == 0)
        return 0; // unknown discretization
    else
        return poly_fmult2x2_numel(deg, D);
}

INT kdv_fscatter(const UINT D, COMPLEX const * const q,
                 const REAL eps_t, COMPLEX * const result, UINT * const deg_ptr,
                 INT * const W_ptr, kdv_discretization_t discretization)
{
    INT ret_code;
    UINT i;
    fnft__akns_discretization_t akns_discretization;
    COMPLEX *r = NULL;
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(u);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(eps_t);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (deg_ptr == NULL)
        return E_INVALID_ARGUMENT(deg_ptr);

    ret_code = kdv_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);   
    
    
    r = malloc(D*sizeof(COMPLEX));
    if (r == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    for (i = 0; i < D; i++)
        r[i] = -1;
    
    ret_code = akns_fscatter(D, q, r, eps_t, result, deg_ptr, W_ptr, akns_discretization);

leave_fun:
    free(r);
    return ret_code;
}
