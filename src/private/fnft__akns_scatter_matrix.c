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
 * Peter J Prins (TU Delft) 2020.
 */
#define FNFT_ENABLE_SHORT_NAMES


#include "fnft__akns_scatter.h"

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
        COMPLEX * const result, akns_discretization_t const discretization,
        const UINT derivative_flag)
{
    
    INT ret_code = SUCCESS;
    UINT i, n, j;
    UINT disc_flag = 0;
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
    
    COMPLEX *l = NULL, qn, rn, ks, k, ch, sh, u1, l_curr,chi, ud1,
            ud2, TM[2][2] = {{0}}, TMD[2][2] = {{0}}, UD[2][2] = {{0}}, UN[2][2] = {{0}};
    REAL scl_factor = 0;
    COMPLEX a1, a2, a3, s, c, w;
    COMPLEX *tmp1 = NULL, *tmp2 = NULL;
    COMPLEX *weights = NULL;

    l = malloc(D*sizeof(COMPLEX));
    if (l == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    COMPLEX l_weights[4] = {0};
    UINT M = 0, N = 0;
    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    switch (discretization) {
        
        case akns_discretization_BO: //  bofetta-osborne scheme
            M = 1; N = 1;
            break;
        case akns_discretization_CF4_2: // commutator-free fourth-order
            M = 2; N = 2;
            break;
        case akns_discretization_CF4_3: // commutator-free fourth-order
            M = 3; N = 3;
            break;
        case akns_discretization_CF5_3: // commutator-free fifth-order
            M = 3; N = 3;
            break;
        case akns_discretization_CF6_4: // commutator-free sixth-order
            M = 4; N = 3;
            break;
            
        default:
            disc_flag = 1;
            break;
            
    }    
    if (disc_flag == 0){
        ret_code = akns_discretization_method_weights(&weights,discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        
        for  (i = 0; i < M; i++){
            for  (j = 0; j < N; j++)
                l_weights[i] = l_weights[i] + weights[i*N+j];
        }
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
                    l[n] = l_curr*l_weights[0];
                break;
                
            case akns_discretization_CF4_3: // commutator-free fourth-order
            case akns_discretization_CF5_3: // commutator-free fifth-order
                if (D%3 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 1.0/3.0;
                for (n = 0; n < D; n = n+3){
                    l[n] = l_curr*l_weights[0];
                    l[n+1] = l_curr*l_weights[1];
                    l[n+2] = l_curr*l_weights[2];
                }
                break;
                
            case akns_discretization_CF6_4: // commutator-free sixth-order
                if (D%4 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 0.25;
                for (n = 0; n < D; n = n+4){
                    l[n] = l_curr*l_weights[0];
                    l[n+1] = l_curr*l_weights[1];
                    l[n+2] = l_curr*l_weights[2];
                    l[n+3] = l_curr*l_weights[3];
                }
                break;
                
            default:
                disc_flag = 1;
                break;
                
        }
        if (disc_flag == 1)
            break;
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
                
                misc_mat_mult_4x4(&U[0][0], &T[0][0]);
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
                misc_mat_mult_2x2(&U[0][0], &T[0][0]);
            }
            
            result[i*4] = T[0][0];
            result[i*4 + 1] = T[0][1];
            result[i*4 + 2] = T[1][0];
            result[i*4 + 3] = T[1][1];
        }
    }
    
    
    // The following two methods have similar structure but different from
    // the ones above. Hence they are grouped together to ensure the
    // implementations are equally efficient.
    
    if (disc_flag == 1){
        REAL eps_t_3 = eps_t*eps_t*eps_t;
        REAL eps_t_2 = eps_t*eps_t;
        COMPLEX w_d, s_d, c_d;
        // Computing values that are functions of q and r but not l.
        // They are stored so that they can be reused for all l.
        switch (discretization) {
            //  Fourth-order exponential method which requires
            // one matrix exponential. The matrix exponential is
            // implmented by using the expansion of the 2x2 matrix
            // in terms of Pauli matrices.
            case akns_discretization_ES4:
                scl_factor = 1;
                tmp1 = malloc(D*sizeof(COMPLEX));
                if (tmp1 == NULL) {
                    ret_code = E_NOMEM;
                    goto leave_fun;
                }
                for (n = 0; n < D; n=n+3){
                    tmp1[n] = eps_t_3*(q[n+2]+r[n+2])/48.0 + (eps_t*(q[n]+r[n]))*0.5;
                    tmp1[n+1] = (eps_t*(q[n]-r[n])*I)*0.5 + (eps_t_3*(q[n+2]-r[n+2])*I)/48.0;
                    tmp1[n+2] = -eps_t_3*(q[n]*r[n+1]- q[n+1]*r[n])/12.0;
                }
                if (derivative_flag == 1){
                    tmp2 = malloc(D*sizeof(COMPLEX));
                    if (tmp2 == NULL) {
                        ret_code = E_NOMEM;
                        goto leave_fun;
                    }
                    for (n = 0; n < D; n=n+3){
                        tmp2[n] = I*eps_t_3*(q[n+1]-r[n+1])/12.0;
                        tmp2[n+1] = -eps_t_3*(q[n+1]+r[n+1])/12.0;
                        tmp2[n+2] = -I*eps_t;
                    }
                    
                }
                break;
                
                //  Fourth-order exponential method which requires
                // three matrix exponentials. The matrix exponential is
                // implmented by using the expansion of the 2x2 matrix
                // in terms of Pauli matrices.
            case akns_discretization_TES4:
                scl_factor = 1;
                tmp1 = malloc(D*sizeof(COMPLEX));
                tmp2 = malloc(D*sizeof(COMPLEX));
                if (tmp1 == NULL || tmp2 == NULL) {
                    ret_code = E_NOMEM;
                    goto leave_fun;
                }
                for (n = 0; n < D; n=n+3){
                    tmp1[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp1[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                    tmp2[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp2[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
                }
                break;
            default: // Unknown discretization
                
                ret_code = E_INVALID_ARGUMENT(discretization);
                goto leave_fun;
                
        }
        
        for (i = 0; i < K; i++) { // iterate over lambda
            l_curr = lambda[i];
            
            if (derivative_flag == 1){
                COMPLEX T[4][4] =
                { {1,0,0,0}, {0,1,0,0},{0,0,1,0},{0,0,0,1} };
                COMPLEX U[4][4] = {{ 0 }};
                
                switch (discretization) {
                    //  Fourth-order exponential method which requires
                    // one matrix exponential. The matrix exponential is
                    // implmented by using the expansion of the 2x2 matrix
                    // in terms of Pauli matrices.
                    case akns_discretization_ES4:
                        for (n = 0; n < D; n=n+3){
                            a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                            a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                            a3 = - eps_t*I*l_curr +tmp1[n+2];
                            w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            w_d = -(1/w)*(a1*tmp2[n]+a2*tmp2[n+1]+a3*tmp2[n+2]);
                            c_d = -CSIN(w)*w_d;
                            s_d = w_d*(c-s)/w;
                            U[0][0] = (c+s*a3);
                            U[0][1] = s*(a1-I*a2);
                            U[1][0] = s*(a1+I*a2);
                            U[1][1] = (c-s*a3);
                            U[2][0] = c_d+s_d*a3+s*tmp2[n+2];
                            U[2][1] = s_d*a1+s*tmp2[n]-I*s_d*a2-I*s*tmp2[n+1];
                            U[3][0] = s_d*a1+s*tmp2[n]+I*s_d*a2+I*s*tmp2[n+1];
                            U[3][1] = c_d-s_d*a3-s*tmp2[n+2];
                            U[2][2] = U[0][0];
                            U[2][3] = U[0][1];
                            U[3][2] = U[1][0];
                            U[3][3] = U[1][1];
                            misc_mat_mult_4x4(&U[0][0], &T[0][0]);
                        }
                        break;
                        // Fourth-order exponential method which requires
                        // three matrix exponentials.
                    case akns_discretization_TES4:
                        for (n = 0; n < D; n=n+3){
                            
                            a1 = tmp1[n];
                            a2 = tmp1[n+1];
//                             a3 = 0;
                            
                            w = CSQRT(-(a1*a1)-(a2*a2));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            TM[0][0] = c;
                            TM[0][1] = s*(a1-I*a2);
                            TM[1][0] = s*(a1+I*a2);
                            TM[1][1] = c;
                            TMD[0][0] = TM[0][0];
                            TMD[0][1] = TM[0][1];
                            TMD[1][0] = TM[1][0];
                            TMD[1][1] = TM[1][1];
                            
                            a1 = (eps_t*(q[n] + r[n]))*0.5;
                            a2 = (eps_t*(q[n]*I - r[n]*I))*0.5;
                            a3 = -eps_t*l_curr*I;
                            w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            UN[0][0] = (c+s*a3);
                            UN[0][1] = s*(a1-I*a2);
                            UN[1][0] = s*(a1+I*a2);
                            UN[1][1] = (c-s*a3);
                            s_d = CSIN(w*eps_t)/w;
                            c_d = -eps_t*l_curr*s_d;
                            w_d = l_curr*(eps_t*w*CCOS(w*eps_t)-CSIN(w*eps_t))/(w*w*w);
                            UD[0][0] = c_d-I*s_d;
                            UD[0][1] = w_d*q[n];
                            UD[1][0] = w_d*r[n];
                            UD[1][1] = c_d+I*s_d;
                            
                            misc_mat_mult_2x2(&UN[0][0], &TM[0][0]);
                            misc_mat_mult_2x2(&UD[0][0], &TMD[0][0]);
                            
                            a1 = tmp2[n];
                            a2 = tmp2[n+1];
//                             a3 = 0;
                            w = CSQRT(-(a1*a1)-(a2*a2));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            UN[0][0] = c;
                            UN[0][1] = s*(a1-I*a2);
                            UN[1][0] = s*(a1+I*a2);
                            UN[1][1] = c;
                            UD[0][0] = UN[0][0];
                            UD[0][1] = UN[0][1];
                            UD[1][0] = UN[1][0];
                            UD[1][1] = UN[1][1];
                            
                            misc_mat_mult_2x2(&UN[0][0], &TM[0][0]);
                            misc_mat_mult_2x2(&UD[0][0], &TMD[0][0]);
                            
                            U[0][0] = TM[0][0];
                            U[0][1] = TM[0][1];
                            U[1][0] = TM[1][0];
                            U[1][1] = TM[1][1];
                            U[2][0] = TMD[0][0];
                            U[2][1] = TMD[0][1];
                            U[3][0] = TMD[1][0];
                            U[3][1] = TMD[1][1];
                            U[2][2] = U[0][0];
                            U[2][3] = U[0][1];
                            U[3][2] = U[1][0];
                            U[3][3] = U[1][1];
                            
                            misc_mat_mult_4x4(&U[0][0], &T[0][0]);
                        }
                        break;
                        
                    default: // Unknown discretization
                        
                        ret_code = E_INVALID_ARGUMENT(discretization);
                        goto leave_fun;
                        
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
                switch (discretization) {
                    //  Fourth-order exponential method which requires
                    // one matrix exponential. The matrix exponential is
                    // implmented by using the expansion of the 2x2 matrix
                    // in terms of Pauli matrices.
                    case akns_discretization_ES4:
                        for (n = 0; n < D; n=n+3){
                            a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                            a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                            a3 = - eps_t*I*l_curr +tmp1[n+2];
                            w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            U[0][0] = (c+s*a3);
                            U[0][1] = s*(a1-I*a2);
                            U[1][0] = s*(a1+I*a2);
                            U[1][1] = (c-s*a3);
                            misc_mat_mult_2x2(&U[0][0], &T[0][0]);
                        }
                        break;
                        // Fourth-order exponential method which requires
                        // three matrix exponentials.
                    case akns_discretization_TES4:
                        j = 0;
                        for (n = 0; n < D; n++){
                            if (j == 0){
                                a1 = tmp1[n];
                                a2 = tmp1[n+1];
                                a3 = 0;
                                j++;
                            }else if (j == 1){
                                a1 = (eps_t*(q[n-1] + r[n-1]))*0.5;
                                a2 = (eps_t*(q[n-1]*I - r[n-1]*I))*0.5;
                                a3 = -eps_t*l_curr*I;
                                j++;
                            }else{
                                a1 = tmp2[n-2];
                                a2 = tmp2[n-1];
                                a3 = 0;
                                j = 0;
                            }
                            w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
                            if (w != 0)
                                s = CSIN(w)/w;
                            else
                                s = 1;
                            c = CCOS(w);
                            U[0][0] = (c+s*a3);
                            U[0][1] = s*(a1-I*a2);
                            U[1][0] = s*(a1+I*a2);
                            U[1][1] = (c-s*a3);
                            misc_mat_mult_2x2(&U[0][0], &T[0][0]);
                        }
                        break;
                        
                    default: // Unknown discretization
                        
                        ret_code = E_INVALID_ARGUMENT(discretization);
                        goto leave_fun;
                        
                }
                
                result[i*4] = T[0][0];
                result[i*4 + 1] = T[0][1];
                result[i*4 + 2] = T[1][0];
                result[i*4 + 3] = T[1][1];
            }
        }
    }
    
    leave_fun:
        free(l);
        free(tmp1);
        free(tmp2);
        free(weights);
        return ret_code;
}
