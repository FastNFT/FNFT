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
#include "fnft__kdv_scatter.h"
#include "fnft__nse_scatter.h"

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

static inline void akns_scatter_bound_states_U_ES4(COMPLEX const a1,
                                                  COMPLEX const a2,
                                                  COMPLEX const a3,
                                                  UINT const derivative_flag,
                                                  COMPLEX * const U,
                                                  COMPLEX * const w,
                                                  COMPLEX * const s,
                                                  COMPLEX * const c)
{
    *w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
    *s = misc_CSINC(*w);
    *c = CCOS(*w);
    U[0]                   = *c + *s*a3;
    U[1]                   = *s*(a1 - I*a2);
    U[derivative_flag?4:2] = *s*(a1 + I*a2);
    U[derivative_flag?5:3] = *c - *s*a3;
}

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 */
INT akns_scatter_bound_states(UINT const D,
                              COMPLEX const * const q,
                              COMPLEX const * const r,
                              REAL const *const T,
                              UINT const K,
                              COMPLEX const * const bound_states,
                              COMPLEX * const a_vals,
                              COMPLEX * const aprime_vals,
                              COMPLEX * const b,
                              akns_discretization_t const discretization,
                              akns_pde_t const PDE,
                              UINT const vanilla_flag,
                              UINT const skip_b_flag)
{
    INT ret_code = SUCCESS;

    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (a_vals == NULL)
        return E_INVALID_ARGUMENT(a);
    if (aprime_vals == NULL)
        return E_INVALID_ARGUMENT(a_prime);
    if (b == NULL)
        return E_INVALID_ARGUMENT(b);
    UINT const upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    REAL const boundary_coeff = akns_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(>discretization);
    if (D%upsampling_factor != 0)
        return E_ASSERTION_FAILED;
    UINT const D_given = D/upsampling_factor;
    if (PDE!=akns_pde_KdV && PDE!=akns_pde_NSE)
        return E_INVALID_ARGUMENT(PDE);

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp4 = NULL;

    // Allocating memory for the array l, which will be used to store xi-samples
    // with the appropriate weight for the discretization.
    // Also allocating memory for storing PHI and PSI at all D_given points as
    // there are required to find the right value of b.
    // First, we will store the values of PHI and its xi-derivative as follows:
    // PSIPHI = [*,*,PHI1[0],PHI2[0],PHI1_D[0],PHI2_D[0],PHI1[1],PHI2[1],PHI1_D[1],PHI2_D[1], ... ,,PHI1[D_given-1],PHI2[D_given-1],PHI1_D[D_given-1],PHI2_D[D_given-1]]
    // Next, we will overwrite the derivatives that we don't need anymore
    // (all except for those at D_given) to store PSI:
    // PSIPHI = [PSI1[0],PSI2[0],PHI1[0],PHI2[0],PSI1[1],PSI2[1],PHI1[1],PHI2[1], ... PSI1[D_given-1],PSI2[D_given-1],PHI1[D_given-1],PHI2[D_given-1],PHI1_D[D_given-1],PHI2_D[D_given-1]]
    // This keeps all vectors in adjacent memory locations, such that we can
    // use matrix-vector multiplication.
    COMPLEX * const l = malloc(D*sizeof(COMPLEX));
    COMPLEX * const PSIPHI = malloc((4*(D_given+1)+2) * sizeof(COMPLEX));
    if (l == NULL || PSIPHI == NULL) {
        ret_code = E_NOMEM;
        CHECK_RETCODE(ret_code, leave_fun);
    }
    COMPLEX * const PSI = &PSIPHI[0];
    COMPLEX * const PHI = &PSIPHI[2];

    // Define stepsize constants that are often needed
    REAL const eps_t = (T[1] - T[0])/(D_given - 1);
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.

    REAL scl_factor = 1.0;
    COMPLEX l_weights[4] = {0};
    UINT N = 0;
    switch (discretization) {
        //  Fourth-order exponential method which requires
        // one matrix exponential. The matrix exponential is
        // implmented by using the expansion of the 2x2 matrix
        // in terms of Pauli matrices.
        case akns_discretization_ES4:
            tmp1 = malloc(2*D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            tmp2 = &tmp1[D];
            for (UINT n=0; n<D; n+=3){
                tmp1[n] = eps_t_3*(q[n+2]+r[n+2])/48.0 + (eps_t*(q[n]+r[n]))*0.5;
                tmp1[n+1] = (eps_t*(q[n]-r[n])*I)*0.5 + (eps_t_3*(q[n+2]-r[n+2])*I)/48.0;
                tmp1[n+2] = -eps_t_3*(q[n]*r[n+1]- q[n+1]*r[n])/12.0;

                tmp2[n] = I*eps_t_3*(q[n+1]-r[n+1])/12.0;
                tmp2[n+1] = -eps_t_3*(q[n+1]+r[n+1])/12.0;
                tmp2[n+2] = -I*eps_t;
            }
            break;

            //  Fourth-order exponential method which requires
            // three matrix exponentials. The matrix exponential is
            // implmented by using the expansion of the 2x2 matrix
            // in terms of Pauli matrices.
        case akns_discretization_TES4:
            tmp1 = skip_b_flag ? malloc(2*D*sizeof(COMPLEX)) : malloc(4*D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            tmp2 = &tmp1[D];
            for (UINT n=0; n<D; n+=3){
                tmp1[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp1[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                tmp2[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp2[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
            }
            if (!skip_b_flag){
                tmp3 = &tmp1[2*D];
                tmp4 = &tmp1[3*D];
                if (tmp3 == NULL || tmp4 == NULL) {
                    ret_code = E_NOMEM;
                    CHECK_RETCODE(ret_code, leave_fun);
                }
                for (UINT n = 0; n < D; n+=3){
                    tmp3[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp3[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                    tmp4[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0  + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp4[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0  + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
                }
            }
            break;

        case akns_discretization_CF4_3:         // commutator-free fourth-order
        case akns_discretization_CF5_3:         // commutator-free fifth-order
        case akns_discretization_CF6_4:         // commutator-free sixth-order
            N++;                               // The previous three discretizations require N=3
            // fall through
        case akns_discretization_CF4_2:         // commutator-free fourth-order
            N++;                               // The previous discretization requires N=2
            // fall through
        case akns_discretization_BO:            // bofetta-osborne scheme
            N++;                               // The previous discretization requires N=1
            scl_factor /= upsampling_factor;
            COMPLEX *weights = NULL;
            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, leave_fun); // if ret_code != SUCCESS, akns_discretization_method_weights frees weights if needed

            for (UINT i = 0; i < upsampling_factor; i++){
                for (UINT j = 0; j < N; j++)
                    l_weights[i] += weights[i*N+j];
            }

            free(weights);
            break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

    for (UINT neig=0; neig<K; neig++) { // iterate over bound states
        COMPLEX l_curr = bound_states[neig];

        for (UINT n=0; n<D; n+=upsampling_factor){
            for (UINT i=0; i<upsampling_factor; i++)
                l[n+i] = l_curr*l_weights[i];
        }

        // Scattering PHI and PHI_D from T[0]-eps_t/2 to T[1]+eps_t/2
        // PHI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        // Set initial condition for PHI in S basis:
        COMPLEX f_S[4];
        f_S[0] = 1.0*CEXP(-I*l_curr*(T[0]-eps_t*boundary_coeff));
        f_S[1] = 0.0;
        f_S[2] = f_S[0]*(-I*(T[0]-eps_t*boundary_coeff));
        f_S[3] = 0.0;

        // Fetch the change of basis matrix from S to the basis of the discretization
        COMPLEX Tmx[4][4];
        UINT derivative_flag = 1;
        ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
        CHECK_RETCODE(ret_code, leave_fun);

        // Calculate the initial condition for PHI in the basis of the discretization
        misc_matrix_mult(4,4,1,&Tmx[0][0],&f_S[0],&PHI[4*0 + 0]);

        // Declaring chonge of state matrices here, to avoid letting them be
        // overwritten with zeros in every loop iteration.
        COMPLEX U[4][4] = {{0}}, UN[4][4] = {{0}};
        switch (discretization) {

            case akns_discretization_BO:
            case akns_discretization_CF4_2:
            case akns_discretization_CF4_3:
            case akns_discretization_CF5_3:
            case akns_discretization_CF6_4:
                {} // Stop the compiler from complaining about starting a case with a declaration rather than a statement.
                COMPLEX phi_temp[4];
                memcpy(phi_temp, PHI, 4 * sizeof(COMPLEX));
                for (UINT n_given=0; n_given<D_given; n_given++) {
                    for (UINT count=0; count<upsampling_factor; count++) {
                        UINT n = n_given * upsampling_factor + count;
                        COMPLEX ks = ((q[n]*r[n])-(l[n]*l[n]));
                        COMPLEX k = CSQRT(ks);
                        COMPLEX ch = CCOSH(k*eps_t);
                        COMPLEX chi = ch/ks;
                        COMPLEX sh = eps_t * misc_CSINC(I*k*eps_t);
                        COMPLEX u1 = l[n]*sh*I;
                        COMPLEX ud1 = eps_t*l[n]*l[n]*chi*I;
                        COMPLEX ud2 = l[n]*(eps_t*ch-sh)/ks;
                        U[0][0] = ch-u1;
                        U[0][1] = q[n]*sh;
                        U[1][0] = r[n]*sh;
                        U[1][1] = ch + u1;
                        U[2][0] = ud1-(l[n]*eps_t+I+(l[n]*l[n]*I)/ks)*sh;
                        U[2][1] = -q[n]*ud2;
                        U[3][0] = -r[n]*ud2;
                        U[3][1] = -ud1-(l[n]*eps_t-I-(l[n]*l[n]*I)/ks)*sh;

                        COMPLEX c = U[2][0]*phi_temp[0] + U[2][1]*phi_temp[1] + U[0][0]*phi_temp[2] + U[0][1]*phi_temp[3];
                        phi_temp[3] = U[3][0]*phi_temp[0] + U[3][1]*phi_temp[1] + U[1][0]*phi_temp[2] + U[1][1]*phi_temp[3];
                        phi_temp[2] = c;

                        c = U[1][0]*phi_temp[0] + U[1][1]*phi_temp[1];
                        phi_temp[0] = U[0][0]*phi_temp[0] + U[0][1]*phi_temp[1];
                        phi_temp[1] = c;
                    }
                    memcpy(&PHI[4*(n_given+1)], phi_temp, 4 * sizeof(COMPLEX));
                }
                break;
                //  Fourth-order exponential method which requires
                // one matrix exponential. The matrix exponential is
                // implmented by using the expansion of the 2x2 matrix
                // in terms of Pauli matrices.
            case akns_discretization_ES4:
                for (UINT n = 0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX w, s, c;
                    COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                    COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                    COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                    akns_scatter_bound_states_U_ES4(a1,a2,a3,1,*U,&w,&s,&c);
                    memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                    COMPLEX w_d = -(a1*tmp2[n]+a2*tmp2[n+1]+a3*tmp2[n+2]);
                    COMPLEX c_d = -misc_CSINC(w)*w_d;
                    w_d /= w;
                    COMPLEX s_d = w_d*(c-s)/w;
                    U[2][0] = c_d+s_d*a3+s*tmp2[n+2];
                    U[2][1] = s_d*a1+s*tmp2[n]-I*s_d*a2-I*s*tmp2[n+1];
                    U[3][0] = s_d*a1+s*tmp2[n]+I*s_d*a2+I*s*tmp2[n+1];
                    U[3][1] = c_d-s_d*a3-s*tmp2[n+2];
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                }
                break;
                // Fourth-order exponential method which requires
                // three matrix exponentials. The transfer metrix
                // needs to be built differently compared to the CF schemes.
            case akns_discretization_TES4:
                for (UINT n=0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX phi_temp[4], w, s, c;

                    akns_scatter_bound_states_U_ES4(tmp1[n],tmp1[n+1],0.0,1,*U,&w,&s,&c);
                    memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);

                    COMPLEX a1 = (eps_t*(q[n] + r[n]))*0.5;
                    COMPLEX a2 = (eps_t*(q[n]*I - r[n]*I))*0.5;
                    COMPLEX a3 = -eps_t*l_curr*I;
                    akns_scatter_bound_states_U_ES4(a1,a2,a3,1,*UN,&w,&s,&c);
                    memcpy(&UN[2][2],&UN[0][0],6 * sizeof(COMPLEX)); // lower right block
                    COMPLEX s_d = eps_t * misc_CSINC(w);
                    COMPLEX c_d = -eps_t*l_curr*s_d;
                    COMPLEX w_d = l_curr*(eps_t*w*CCOS(w*eps_t)-CSIN(w*eps_t))/(w*w*w);
                    UN[2][0] = c_d-I*s_d;
                    UN[2][1] = w_d*q[n];
                    UN[3][0] = w_d*r[n];
                    UN[3][1] = c_d+I*s_d;
                    misc_matrix_mult(4,4,1,*UN,&PHI[4*(n_given+1)],phi_temp);

                    akns_scatter_bound_states_U_ES4(tmp2[n],tmp2[n+1],0.0,1,*U,&w,&s,&c);
                    memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                    misc_matrix_mult(4,4,1,*U,phi_temp,&PHI[4*(n_given+1)]);
                }
                break;

            default: // Unknown discretization
                ret_code = E_INVALID_ARGUMENT(discretization);
                CHECK_RETCODE(ret_code, leave_fun);
        }

        // If b-coefficient is requested skip_b_flag will not be set.
        // Scattering PSI from T[1]+eps_t/2 to T[0]-eps_t/2.
        // PSI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        if (skip_b_flag == 0){

            // Set final condition for PSI in S basis:
            f_S[0] = 0.0;
            f_S[1] = 1.0*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));

            // Fetch the change of basis matrix from S to the basis of the discretization
            derivative_flag = 0;
            ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Calculate the final condition for PSI in the basis of the discretization
            misc_matrix_mult(2,2,1,&Tmx[0][0],&f_S[0],&PSI[4*D_given+0]);

            // Inverse transfer matrix at each step is built by taking
            // negative step -eps_t.
            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    {} // Stop the compiler from complaining about starting a case with a declaration rather than a statement.
                    COMPLEX psi_temp[2];
                    memcpy(psi_temp, &PSI[4*D_given], 2 * sizeof(COMPLEX));
                    for (UINT n_given=D_given; n_given-->0; ) {
                        for (UINT count=upsampling_factor; count-->0; ) {
                            COMPLEX U[2][2];
                            UINT n = n_given * upsampling_factor + count;
                            COMPLEX ks = ((q[n]*r[n])-(l[n]*l[n]));
                            COMPLEX k = CSQRT(ks);
                            COMPLEX ch = CCOSH(k*eps_t);
                            COMPLEX sh = -eps_t * misc_CSINC(I*k*eps_t);
                            COMPLEX u1 = l[n]*sh*I;

                            U[0][0] = ch - u1;
                            U[0][1] = q[n]*sh;
                            U[1][0] = r[n]*sh;
                            U[1][1] = ch + u1;

                            COMPLEX c = U[1][0]*psi_temp[0] + U[1][1]*psi_temp[1];
                            psi_temp[0] = U[0][0]*psi_temp[0] + U[0][1]*psi_temp[1];
                            psi_temp[1] = c;
                        }
                        memcpy(&PSI[4*n_given], psi_temp, 2 * sizeof(COMPLEX));
                    }
                    break;

                    //  Fourth-order exponential method which requires
                    // one matrix exponential. The matrix exponential is
                    // implmented by using the expansion of the 2x2 matrix
                    // in terms of Pauli matrices.
                case akns_discretization_ES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2], w, s, c;
                        COMPLEX a1 = -tmp1[n]- eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = -tmp1[n+1] + eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 =  eps_t*I*l_curr -tmp1[n+2];
                        akns_scatter_bound_states_U_ES4(a1,a2,a3,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                    }
                    break;
                    // Fourth-order exponential method which requires
                    // three matrix exponentials. The transfer metrix cannot
                    // needs to be built differently compared to the CF schemes.
                case akns_discretization_TES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2], psi_temp[2], w, s, c;

                        akns_scatter_bound_states_U_ES4(tmp3[n],tmp3[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);

                        akns_scatter_bound_states_U_ES4(-0.5*eps_t*(q[n]+r[n]), -0.5*I*eps_t*(q[n]-r[n]),eps_t*l_curr*I,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*n_given],psi_temp);

                        akns_scatter_bound_states_U_ES4(tmp4[n],tmp4[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,psi_temp,&PSI[4*n_given]);
                    }
                    break;

                default: // Unknown discretization

                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
        }

        // Fetch the change of basis matrix from the basis of the discretization to S
        derivative_flag = 1;
        ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
        CHECK_RETCODE(ret_code, leave_fun);

        // Calculate the final state for PHI in S basis
        misc_matrix_mult(4,4,1,&Tmx[0][0],&PHI[4*D_given + 0],&f_S[0]);

        // Calculate the final state for PHI in E basis
        a_vals[neig] = f_S[0] * CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));
        if (PDE==akns_pde_KdV)
            a_vals[neig] = CREAL(a_vals[neig]);
        aprime_vals[neig] = scl_factor*(f_S[2]*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff))+(I*(T[1]+eps_t*boundary_coeff))* a_vals[neig]);
        if (PDE==akns_pde_KdV)
            aprime_vals[neig] = I * CIMAG(aprime_vals[neig]);

        if (skip_b_flag == 0){
            // Calculation of b assuming a=0
            // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
            // computation point

            // Fetch the change of basis matrix from the basis of the discretization to S
            derivative_flag = 0;
            ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            REAL error_metric = INFINITY, tmp = INFINITY;
            for (UINT n = 0; n <= D_given; n++){
                COMPLEX * const psi_S = &f_S[0], * const phi_S = &f_S[2], b_temp[2];
                // Calculate phi and psi for this sample in S basis
                misc_matrix_mult(2,2,1,&Tmx[0][0],&PHI[4*n + 0], phi_S);
                misc_matrix_mult(2,2,1,&Tmx[0][0],&PSI[4*n + 0], psi_S);
                if (PDE==akns_pde_KdV) {
                    for (UINT i=0; i<4; i++)
                        f_S[i] = CREAL(f_S[i]);
                }
                for (UINT i=0; i<2; i++)
                    b_temp[i] = phi_S[i]/psi_S[i];
                if (PDE!=akns_pde_KdV || CREAL(b_temp[0]*b_temp[1])>0) {
                    tmp = FABS( 0.5* LOG( (REAL)CABS( b_temp[1]/b_temp[0] ) ) );
                    if (tmp < error_metric){
                        b[neig] = b_temp[0];
                        error_metric = tmp;
                    }
                }
            }
        }

    }
leave_fun:
    free(tmp1);
    free(l);
    free(PSIPHI);
    return ret_code;
}
