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
* Shrinivas Chimmalgi (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__kdv_scatter.h"
#include <stdio.h>

/**
 * Returns [S11 S12 S21 S22 S11' S12' S21' S22'] in result 
 * where S = [S11, S12; S21, S22] is the scattering matrix.
 * computed using the chosen scheme.
 * Result should be preallocated with size 8 * K
 * Default scheme should be set as BO.
 */
INT kdv_scatter_matrix(const UINT D, COMPLEX const * const q,
    const REAL eps_t, const UINT K, 
    COMPLEX const * const lambda,
    COMPLEX * const result, kdv_discretization_t discretization)
{
     
    INT ret_code = SUCCESS;
    UINT i;
    akns_discretization_t akns_discretization;
    COMPLEX *r = NULL;
    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (!(eps_t > 0))
        return E_INVALID_ARGUMENT(eps_t);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (lambda == NULL)
        return E_INVALID_ARGUMENT(lambda);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);

    ret_code = kdv_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);   
    
    
    r = malloc(D*sizeof(COMPLEX));
    if (r == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    for (i = 0; i < D; i++)
        r[i] = -1;

    
    ret_code = akns_scatter_matrix(D, q, r, eps_t, K, lambda, result, akns_discretization);

leave_fun:
    free(r);
    return ret_code;
}
