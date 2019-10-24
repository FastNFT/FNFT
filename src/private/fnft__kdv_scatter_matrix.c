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
    
    UINT D_scale = akns_discretization_D_scale(akns_discretization);
    
    r = malloc(D*sizeof(COMPLEX));
    if (r == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    
    
    switch (akns_discretization) {
        
        case akns_discretization_BO: // Bofetta-Osborne scheme
            for (i = 0; i < D; i++)
                r[i] = -1;
            
            break;
            
        case akns_discretization_CF4_2:
            // commutator-free fourth-order two exponentials
            if (D%2 != 0){
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            for (i = 0; i < D; i++)
                r[i] = -0.5;
            break;
        case akns_discretization_CF4_3:
            // commutator-free fourth-order three exponentials
            if (D%3 != 0){
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            for (i = 0; i < D; i = i+3)
                r[i] = -0.275;
            for (i = 1; i < D; i = i+3)
                r[i] = -0.45;
            for (i = 2; i < D; i = i+3)
                r[i] = -0.275;
            break;
            
        case akns_discretization_CF5_3:
            // commutator-free fifth-order three exponentials
            if (D%3 != 0){
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            for (i = 0; i < D; i = i+3)
                r[i] = -0.3-0.1*I;
            for (i = 1; i < D; i = i+3)
                r[i] = -0.4;
            for (i = 2; i < D; i = i+3)
                r[i] = -0.3+0.1*I;
            break;
            
        case akns_discretization_CF6_4:
            // commutator-free sixth-order four exponentials
            if (D%4 != 0){
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            for (i = 0; i < D; i = i+4)
                r[i] = -0.210073786808785 - 0.046600721949282*I;
            for (i = 1; i < D; i = i+4)
                r[i] = -0.289926213191215 + 0.046600721949282*I;
            for (i = 2; i < D; i = i+4)
                r[i] = -0.289926213191215 + 0.046600721949282*I;
            for (i = 3; i < D; i = i+4)
                r[i] = -0.210073786808785 - 0.046600721949282*I;
            break;
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(discretization);
            goto leave_fun;
    }
    
    
    ret_code = akns_scatter_matrix(D, q, r, eps_t, K, lambda, result, akns_discretization);
    
    leave_fun:
        free(r);
        return ret_code;
}
