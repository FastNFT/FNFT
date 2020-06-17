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
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__akns_scatter.h"
#include <stdio.h>

/**
 * If derivative_flag=0 returns [S11 S12 S21 S22] in result where
 * S = [S11, S12; S21, S22] is the scattering matrix computed using the 
 * chosen scheme.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in 
 * result where S11' is the derivative of S11 w.r.t to lambda.
 * Result should be preallocated with size 4*K or 8*K accordingly.
 */
INT akns_scatter_matrix(const UINT D, COMPLEX const * const q,
        COMPLEX const * const r, const REAL eps_t,
        const UINT K, COMPLEX const * const lambda,
        COMPLEX * const result, akns_discretization_t discretization,
        const UINT derivative_flag)
{
    
    INT ret_code = SUCCESS;
    UINT i, n;
    INT  c1, c2, c3;
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (!(eps_t > 0))
        return E_INVALID_ARGUMENT(eps_t);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (lambda == NULL)
        return E_INVALID_ARGUMENT(lambda);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (derivative_flag != 0 && derivative_flag != 1)
        return E_INVALID_ARGUMENT(derivative_flag);
    
    COMPLEX *l = NULL, qn, rn, ks, k, sum=0, ch, sh, u1, l_curr,chi, ud1, 
            TM[4][4], ud2;
    REAL scl_factor = 0;    
    
    l = malloc(D*sizeof(COMPLEX));
    if (l == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < K; i++) { // iterate over lambda
        l_curr = lambda[i];
        switch (discretization) {
            
            case akns_discretization_BO: //  bofetta-osborne scheme
                scl_factor = 1;
                for (n = 0; n < D; n++)
                    l[n] = l_curr;
                break;
                
            case akns_discretization_CF4_2: // commutator-free fourth-order
                if (D%2 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 0.5;
                for (n = 0; n < D; n++)
                    l[n] = l_curr/2.0;
                break;
                
            case akns_discretization_CF4_3: // commutator-free fourth-order
                if (D%3 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 1.0/3.0;
                for (n = 0; n < D; n = n+3)
                    l[n] = l_curr*0.275;
                for (n = 1; n < D; n = n+3)
                    l[n] = l_curr*0.45;
                for (n = 2; n < D; n = n+3)
                    l[n] = l_curr*0.275;
                break;
                
            case akns_discretization_CF5_3: // commutator-free fifth-order
                if (D%3 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 1.0/3.0;
                for (n = 0; n < D; n = n+3)
                    l[n] = (0.3+0.1*I)*l_curr;
                for (n = 1; n < D; n = n+3)
                    l[n] = l_curr*0.4;
                for (n = 2; n < D; n = n+3)
                    l[n] = (0.3-0.1*I)*l_curr;
                break;
                
            case akns_discretization_CF6_4: // commutator-free sixth-order
                if (D%4 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 0.25;
                for (n = 0; n < D; n = n+4)
                    l[n] = (0.210073786808785 + 0.046600721949282*I)*l_curr;
                for (n = 1; n < D; n = n+4)
                    l[n] = (0.289926213191215 - 0.046600721949282*I)*l_curr;
                for (n = 2; n < D; n = n+4)
                    l[n] = (0.289926213191215 - 0.046600721949282*I)*l_curr;
                for (n = 3; n < D; n = n+4)
                    l[n] = (0.210073786808785 + 0.046600721949282*I)*l_curr;
                break;
                
            default: // Unknown discretization
                
                ret_code = E_INVALID_ARGUMENT(discretization);
                goto leave_fun;
                
        }
        
        if (derivative_flag == 1){
            COMPLEX T[4][4] =
            { {1,0,0,0}, {0,1,0,0},{0,0,1,0},{0,0,0,1} };
            COMPLEX U[4][4] = {{ 0 }};
            
            for (n = 0; n < D; n++){
                qn = q[n];
                rn = r[n];
                ks = ((q[n]*r[n])-(l[n]*l[n]));
                k = CSQRT(ks);
                ch = CCOSH(k*eps_t);
                chi = ch/ks;
                if (ks != 0)
                    sh = CSINH(k*eps_t)/k;
                else
                    sh = eps_t;
                u1 = l[n]*sh*I;
                ud1 = eps_t*l[n]*l[n]*chi*I;
                ud2 = l[n]*(eps_t*ch-sh)/ks;
                
                U[0][0] = ch-u1;
                U[0][1] = qn*sh;
                U[1][0] = rn*sh;
                U[1][1] = ch + u1;
                U[2][0] = ud1-(l[n]*eps_t+I+(l[n]*l[n]*I)/ks)*sh;
                U[2][1] = -qn*ud2;
                U[3][0] = -rn*ud2;
                U[3][1] = -ud1-(l[n]*eps_t-I-(l[n]*l[n]*I)/ks)*sh;
                U[2][2] = ch-u1;
                U[2][3] = qn*sh;
                U[3][2] = rn*sh;
                U[3][3] = ch+u1;
                
                
                for (c1 = 0; c1 < 4; c1++) {
                    for (c2 = 0; c2 < 4; c2++) {
                        for (c3 = 0; c3 < 4; c3++) {
                            sum = sum + U[c1][c3]*T[c3][c2];
                        }
                        TM[c1][c2] = sum;
                        sum = 0;
                    }
                }
                for (c1 = 0; c1 < 4; c1++) {
                    for (c2 = 0; c2 < 4; c2++)
                        T[c1][c2] = TM[c1][c2];
                }
            }
            
            result[i*8] = T[0][0];
            result[i*8 + 1] = T[0][1];
            result[i*8 + 2] = T[1][0];
            result[i*8 + 3] = T[1][1];
            result[i*8 + 4] = T[2][0]*scl_factor;
            result[i*8 + 5] = T[2][1]*scl_factor;
            result[i*8 + 6] = T[3][0]*scl_factor;
            result[i*8 + 7] = T[3][1]*scl_factor;
            
        }else{
            COMPLEX T[2][2] =
            { {1,0}, {0,1}};
            COMPLEX U[2][2] = {{ 0 }};
            for (n = 0; n < D; n++){
                qn = q[n];
                rn = r[n];
                ks = ((q[n]*r[n])-(l[n]*l[n]));
                k = CSQRT(ks);
                ch = CCOSH(k*eps_t);
                if (ks != 0)
                    sh = CSINH(k*eps_t)/k;
                else
                    sh = eps_t;
                
                u1 = l[n]*sh*I;
                U[0][0] = ch-u1;
                U[0][1] = qn*sh;
                U[1][0] = rn*sh;
                U[1][1] = ch + u1;
                
                for (c1 = 0; c1 < 2; c1++) {
                    for (c2 = 0; c2 < 2; c2++) {
                        for (c3 = 0; c3 < 2; c3++) {
                            sum = sum + U[c1][c3]*T[c3][c2];
                        }
                        TM[c1][c2] = sum;
                        sum = 0;
                    }
                }
                for (c1 = 0; c1 < 2; c1++) {
                    for (c2 = 0; c2 < 2; c2++)
                        T[c1][c2] = TM[c1][c2];
                }
            }
            
            result[i*4] = T[0][0];
            result[i*4 + 1] = T[0][1];
            result[i*4 + 2] = T[1][0];
            result[i*4 + 3] = T[1][1];
        }
    }

    
    leave_fun:
        free(l);
        return ret_code;
}
