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
 * Shrinivas Chimmalgi (TU Delft) 2017, 2019-2020.
 * Marius Brehler (TU Dortmund) 2018.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"


/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 */
INT nse_scatter_bound_states(const UINT D, COMPLEX const *const q,
        COMPLEX *r, REAL const *const T, UINT K,
        COMPLEX *bound_states, COMPLEX *a_vals,
        COMPLEX *aprime_vals, COMPLEX *b,
        nse_discretization_t discretization, UINT skip_b_flag)
{
    
    INT ret_code = SUCCESS;
    UINT neig;
    UINT n, upsampling_factor, D_given, n_given, count;
    COMPLEX * l = NULL, qn, rn, ks, k,ch,chi,sh,u1,ud1,ud2, l_curr;
    REAL eps_t_n = 0, eps_t = 0, scl_factor = 0;
    COMPLEX * PHI1 = NULL, * PHI2 = NULL;
    COMPLEX * PSI1 = NULL, * PSI2 = NULL;
    COMPLEX PHI1_D = 0, PHI2_D = 0;
    COMPLEX phi1 = 0, phi2 = 0, psi1 = 0, psi2 = 0;
    
    COMPLEX a1, a2, a3, s, c, w;
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp4 = NULL;
    COMPLEX w_d, s_d, c_d, TM[2][2] = {{0}}, TMD[2][2] = {{0}}, UD[2][2] = {{0}}, UN[2][2] = {{0}};
    UINT disc_flag = 0;
    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL)
        return E_INVALID_ARGUMENT(eps_t);
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
    
    
    if (r == NULL) {
        r = malloc(D*sizeof(COMPLEX));
        if (r == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (n = 0; n < D; n++)
            r[n] = -CONJ(q[n]);
        
    }
    l = malloc(D*sizeof(COMPLEX));
    if (l == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    upsampling_factor = nse_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0){
        ret_code =  E_INVALID_ARGUMENT(discretization);
        goto leave_fun;
    }
    D_given = D/upsampling_factor;
    
    // Allocating memory for storing PHI and PSI at all D_given points as
    // there are required to find the right value of b.
    PHI1 = malloc((D_given+1) * sizeof(COMPLEX));
    PHI2 = malloc((D_given+1) * sizeof(COMPLEX));
    PSI1 = malloc((D_given+1) * sizeof(COMPLEX));
    PSI2 = malloc((D_given+1) * sizeof(COMPLEX));
    if (PHI1 == NULL || PHI2 == NULL || PSI1 == NULL || PSI2 == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    
    eps_t = (T[1] - T[0])/(D_given - 1);
    eps_t_n = -eps_t;
    REAL eps_t_3 = eps_t*eps_t*eps_t;
    REAL eps_t_2 = eps_t*eps_t;
    
    REAL boundary_coeff;
    boundary_coeff = nse_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN){
        ret_code = E_INVALID_ARGUMENT(>discretization);
        goto leave_fun;
    }
    
    COMPLEX *weights = NULL;
    COMPLEX l_weights[4] = {0};
    UINT M = 0, N = 0, j = 0, i = 0;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l. 
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.
    
    switch (discretization) {
        //  Fourth-order exponential method which requires
        // one matrix exponential. The matrix exponential is
        // implmented by using the expansion of the 2x2 matrix
        // in terms of Pauli matrices.
        case akns_discretization_ES4:
            scl_factor = 1.0;
            tmp1 = malloc(D*sizeof(COMPLEX));
            tmp2 = malloc(D*sizeof(COMPLEX));
            if (tmp1 == NULL ||tmp2 == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            for (n = 0; n < D; n=n+3){
                tmp1[n] = eps_t_3*(q[n+2]+r[n+2])/48.0 + (eps_t*(q[n]+r[n]))*0.5;
                tmp1[n+1] = (eps_t*(q[n]-r[n])*I)*0.5 + (eps_t_3*(q[n+2]-r[n+2])*I)/48.0;
                tmp1[n+2] = -eps_t_3*(q[n]*r[n+1]- q[n+1]*r[n])/12.0;
                
                tmp2[n] = I*eps_t_3*(q[n+1]-r[n+1])/12.0;
                tmp2[n+1] = -eps_t_3*(q[n+1]+r[n+1])/12.0;
                tmp2[n+2] = -I*eps_t;
            }
            
            disc_flag = 1;
            break;
            
            //  Fourth-order exponential method which requires
            // three matrix exponentials. The matrix exponential is
            // implmented by using the expansion of the 2x2 matrix
            // in terms of Pauli matrices.
        case akns_discretization_TES4:
            scl_factor = 1.0;
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
            if (skip_b_flag == 0){
                tmp3 = malloc(D*sizeof(COMPLEX));
                tmp4 = malloc(D*sizeof(COMPLEX));
                if (tmp3 == NULL || tmp4 == NULL) {
                    ret_code = E_NOMEM;
                    goto leave_fun;
                }
                for (n = 0; n < D; n=n+3){
                    tmp3[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp3[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                    tmp4[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0  + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp4[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0  + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
                }
            }
            disc_flag = 1;
            break;
            
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
            
        default: // Unknown discretization
            
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
    
    
    for (neig = 0; neig < K; neig++) { // iterate over bound states
        l_curr = bound_states[neig];
        
        switch (discretization) {
            
            case nse_discretization_BO: //  bofetta-osborne scheme
                scl_factor = 1;
                for (n = 0; n < D; n++)
                    l[n] = l_curr;
                break;
                
            case nse_discretization_CF4_2: // commutator-free fourth-order
                if (D%2 != 0){
                    ret_code = E_ASSERTION_FAILED;
                    goto leave_fun;
                }
                scl_factor = 0.5;
                for (n = 0; n < D; n++)
                    l[n] = l_curr*l_weights[0];
                break;
                
            case nse_discretization_CF4_3: // commutator-free fourth-order
            case nse_discretization_CF5_3: // commutator-free fifth-order

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
  
            case nse_discretization_CF6_4: // commutator-free sixth-order
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
                if (disc_flag == 0){
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    goto leave_fun;
                }
                
        }
        
        // Scattering PHI and PHI_D from T[0] to T[1]
        // PHI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        COMPLEX U[4][4] = {{0}};        
        PHI1[0] = 1.0*CEXP(-I*l_curr*(T[0]-eps_t*boundary_coeff));
        PHI2[0] = 0.0;
        PHI1_D = PHI1[0]*(-I*(T[0]-eps_t*boundary_coeff));
        PHI2_D = 0.0;
        
        n_given = 0;
        switch (discretization) {
            
            case nse_discretization_BO:
            case nse_discretization_CF4_2:
            case nse_discretization_CF4_3:
            case nse_discretization_CF5_3:
            case nse_discretization_CF6_4:
                count = upsampling_factor - 1;
                phi1 = PHI1[0];
                phi2 = PHI2[0];
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
                    c = U[2][0]*phi1 + U[2][1]*phi2 + U[0][0]*PHI1_D + U[0][1]*PHI2_D;
                    PHI2_D = U[3][0]*phi1 + U[3][1]*phi2 + U[1][0]*PHI1_D + U[1][1]*PHI2_D;
                    PHI1_D = c;
                    
                    c = U[1][0]*phi1 + U[1][1]*phi2;
                    phi1 = U[0][0]*phi1 + U[0][1]*phi2;
                    phi2 = c;
                    
                    if (count == 0){
                        count = upsampling_factor - 1;                        
                        PHI1[n_given+1] = phi1;
                        PHI2[n_given+1] = phi2;
                        n_given++;
                    }
                    else
                        count--;
                    
                }
                break;
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

                    c = U[2][0]*PHI1[n_given] + U[2][1]*PHI2[n_given] + U[0][0]*PHI1_D + U[0][1]*PHI2_D;
                    PHI2_D = U[3][0]*PHI1[n_given] + U[3][1]*PHI2[n_given] + U[1][0]*PHI1_D + U[1][1]*PHI2_D;
                    PHI1_D = c;
                    
                    c = U[1][0]*PHI1[n_given] + U[1][1]*PHI2[n_given];
                    PHI1[n_given+1] = U[0][0]*PHI1[n_given] + U[0][1]*PHI2[n_given];
                    PHI2[n_given+1] = c;
                    
                    n_given++;
                }
                break;
                // Fourth-order exponential method which requires
                // three matrix exponentials. The transfer metrix
                // needs to be built differently compared to the CF schemes.
            case akns_discretization_TES4:
                for (n = 0; n < D; n=n+3){
                    
                    a1 = tmp1[n];
                    a2 = tmp1[n+1];
//                     a3 = 0;
                    
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
//                     a3 = 0;
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

                    c = U[2][0]*PHI1[n_given] + U[2][1]*PHI2[n_given] + U[0][0]*PHI1_D + U[0][1]*PHI2_D;
                    PHI2_D = U[3][0]*PHI1[n_given] + U[3][1]*PHI2[n_given] + U[1][0]*PHI1_D + U[1][1]*PHI2_D;
                    PHI1_D = c;
                    
                    c = U[1][0]*PHI1[n_given] + U[1][1]*PHI2[n_given];
                    PHI1[n_given+1] = U[0][0]*PHI1[n_given] + U[0][1]*PHI2[n_given];
                    PHI2[n_given+1] = c;
                    
                    n_given++;

                }
                break;
                
            default: // Unknown discretization
                
                ret_code = E_INVALID_ARGUMENT(discretization);
                goto leave_fun;
                
        }
        
        
        // If b-coefficient is requested skip_b_flag will not be set.
        // Scattering PSI from T[1] to T[0].
        // PSI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        if (skip_b_flag == 0){
            
            PSI1[D_given] = 0.0;
            PSI2[D_given] = 1.0*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));
            psi1 = PSI1[D_given];
            psi2 = PSI2[D_given];
            // Inverse transfer matrix at each step is built by taking
            // negative step eps_t_n=-eps_t_n. Some quantities pre-calculated
            // for eps_t can be used for eps_t_n with only a sign change.
            switch (discretization) {
                case nse_discretization_BO:
                case nse_discretization_CF4_2:
                case nse_discretization_CF4_3:
                case nse_discretization_CF5_3:
                case nse_discretization_CF6_4:
                    n = D;
                    n_given = D_given; 
                    count = upsampling_factor-1;
                    do{
                        n--;
                        qn = q[n];
                        rn = r[n];
                        ks = ((q[n]*r[n])-(l[n]*l[n]));
                        k = CSQRT(ks);
                        ch = CCOSH(k*eps_t_n);
                        chi = ch/ks;
                        if (ks != 0)
                            sh = CSINH(k*eps_t_n)/k;
                        else
                            sh = eps_t_n;
                        u1 = l[n]*sh*I;
                        U[0][0] = ch-u1;
                        U[0][1] = qn*sh;
                        U[1][0] = rn*sh;
                        U[1][1] = ch + u1;
                        
                        c = U[1][0]*psi1 + U[1][1]*psi2;
                        psi1 = U[0][0]*psi1 + U[0][1]*psi2;
                        psi2 = c;
                        
                        if (count == 0){
                            count = upsampling_factor - 1;
                            PSI1[n_given-1] = psi1;
                            PSI2[n_given-1] = psi2;
                            n_given--;
                        }
                        else
                            count--;
                    } while (n > 0);
                    break;
                    //  Fourth-order exponential method which requires
                    // one matrix exponential. The matrix exponential is
                    // implmented by using the expansion of the 2x2 matrix
                    // in terms of Pauli matrices.
                case akns_discretization_ES4:
                    n = D;
                    n_given = D_given;
                    do{
                        n = n-3;
                        a1 = -tmp1[n]- eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        a2 = -tmp1[n+1] + eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        a3 =  eps_t*I*l_curr -tmp1[n+2];
                        w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
                        if (w != 0)
                            s = CSIN(w)/w;
                        else
                            s = 1;
                        c = CCOS(w);
                        w_d = -(1/w)*(-a1*tmp2[n]-a2*tmp2[n+1]-a3*tmp2[n+2]);
                        c_d = -CSIN(w)*w_d;
                        s_d = w_d*(c-s)/w;
                        U[0][0] = (c+s*a3);
                        U[0][1] = s*(a1-I*a2);
                        U[1][0] = s*(a1+I*a2);
                        U[1][1] = (c-s*a3);
                        
                    c = U[1][0]*PSI1[n_given] + U[1][1]*PSI2[n_given];
                    PSI1[n_given-1] = U[0][0]*PSI1[n_given] + U[0][1]*PSI2[n_given];
                    PSI2[n_given-1] = c;                    
                    n_given--;
                    
                    } while (n > 0);
                    break;
                    // Fourth-order exponential method which requires
                    // three matrix exponentials. The transfer metrix cannot
                    // needs to be built differently compared to the CF schemes.
                case akns_discretization_TES4:
                     n = D;
                    n_given = D_given;
                    do{
                        n = n-3;
                        
                        a1 = tmp3[n];
                        a2 = tmp3[n+1];
//                         a3 = 0;                        
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

                        a1 = (eps_t_n*(q[n] + r[n]))*0.5;
                        a2 = (eps_t_n*(q[n]*I - r[n]*I))*0.5;
                        a3 = -eps_t_n*l_curr*I;
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
                        
                        misc_mat_mult_2x2(&UN[0][0], &TM[0][0]);
                        
                        a1 = tmp4[n];
                        a2 = tmp4[n+1];
//                         a3 = 0;
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
                        
                        misc_mat_mult_2x2(&UN[0][0], &TM[0][0]);
                        
                        U[0][0] = TM[0][0];
                        U[0][1] = TM[0][1];
                        U[1][0] = TM[1][0];
                        U[1][1] = TM[1][1];
   
                        c = U[1][0]*PSI1[n_given] + U[1][1]*PSI2[n_given];
                        PSI1[n_given-1] = U[0][0]*PSI1[n_given] + U[0][1]*PSI2[n_given];
                        PSI2[n_given-1] = c;
                        n_given--;
                    } while (n > 0);
                    break;
                    
                default: // Unknown discretization
                    
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    goto leave_fun;
                    
            }
        }
        
        a_vals[neig] = PHI1[D_given]*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));
        aprime_vals[neig] = scl_factor*(PHI1_D*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff))+(I*(T[1]+eps_t*boundary_coeff))* a_vals[neig]);
        
        if (skip_b_flag == 0){
            // Calculation of b assuming a=0
            // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
            // computation point
            REAL error_metric = INFINITY, tmp = INFINITY;
            for (n = 0; n <= D_given; n++){
                tmp = FABS((0.5*LOG(CABS((PHI2[n]/PSI2[n])/(PHI1[n]/PSI1[n])))));
                if (tmp < error_metric){
                    b[neig] = PHI1[n]/PSI1[n];
                    error_metric = tmp;
                }
            }
        }
        
    }
    leave_fun:
        free(PHI1);
        free(PSI1);
        free(PHI2);
        free(PSI2);
        free(tmp1);
        free(tmp2);
        free(l);
        return ret_code;
}

