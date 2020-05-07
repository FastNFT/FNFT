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
 * Shrinivas Chimmalgi (TU Delft) 2017, 2019.
 * Marius Brehler (TU Dortmund) 2018.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__nse_scatter.h"
#include <stdio.h>
#include "fnft__misc.h"

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 * Default scheme should be set as BO.
 */
INT nse_scatter_bound_states(const UINT D, COMPLEX const *const q,
        COMPLEX *r, REAL const *const T, UINT K,
        COMPLEX *bound_states, COMPLEX *a_vals,
        COMPLEX *aprime_vals, COMPLEX *b,
        nse_discretization_t discretization, UINT skip_b_flag)
{
    
    INT ret_code = SUCCESS;
    UINT neig;
    UINT c1, c2, c3;
    UINT n, D_scale, D_given, n_given, count;
    COMPLEX * l, qn, rn, ks, k, sum=0, TM[4][4],ch,chi,sh,u1,ud1,ud2, l_curr;
    REAL eps_t_n = 0, eps_t = 0, scl_factor = 0;
    COMPLEX * PHI1 = NULL, * PHI2 = NULL;
    COMPLEX * PSI1 = NULL, * PSI2 = NULL;
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
    
    D_scale = nse_discretization_D_scale(discretization);
    if (D_scale == 0){
        ret_code =  E_INVALID_ARGUMENT(discretization);
        goto leave_fun;
    }
    D_given = D/D_scale;
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
                    l[n] = l_curr/2.0;
                break;
                
            case nse_discretization_CF4_3: // commutator-free fourth-order
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
                
            case nse_discretization_CF5_3: // commutator-free fifth-order
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
                
            case nse_discretization_CF6_4: // commutator-free sixth-order
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
        
        
        COMPLEX SR[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        COMPLEX SL[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
        COMPLEX U[4][4] = {{0}};
        
        PSI1[D_given] = SR[0][1]*CEXP(I*l_curr*(T[1]+eps_t/2));
        PSI2[D_given] = SR[1][1]*CEXP(I*l_curr*(T[1]+eps_t/2));
        
        if (skip_b_flag == 0){
            n = D;
            n_given = D_given;
            count = D_scale-1;
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
                ud1 = eps_t*l[n]*l[n]*chi*I;
                ud2 = l[n]*(eps_t_n*ch-sh)/ks;
                U[0][0] = ch-u1;
                U[0][1] = qn*sh;
                U[1][0] = rn*sh;
                U[1][1] = ch + u1;
                U[2][0] = ud1-(l[n]*eps_t_n+I+(l[n]*l[n]*I)/ks)*sh;
                U[2][1] = -qn*ud2;
                U[3][0] = -rn*ud2;
                U[3][1] = -ud1-(l[n]*eps_t_n-I-(l[n]*l[n]*I)/ks)*sh;
                U[2][2] = ch-u1;
                U[2][3] = qn*sh;
                U[3][2] = rn*sh;
                U[3][3] = ch+u1;
                
                for (c1 = 0; c1 < 4; c1++) {
                    for (c2 = 0; c2 < 4; c2++) {
                        for (c3 = 0; c3 < 4; c3++) {
                            sum = sum + U[c1][c3]*SR[c3][c2];
                        }
                        TM[c1][c2] = sum;
                        sum = 0;
                    }
                }
                for (c1 = 0; c1 < 4; c1++) {
                    for (c2 = 0; c2 < 4; c2++)
                        SR[c1][c2] = TM[c1][c2];
                }
                
                if (count == 0){
                    count = D_scale - 1;
                    n_given--;
                    PSI1[n_given] = SR[0][1]*CEXP(I*l_curr*(T[1]+eps_t/2));
                    PSI2[n_given] = SR[1][1]*CEXP(I*l_curr*(T[1]+eps_t/2));
                }
                else
                    count--;
                
                
            } while (n > 0);
        }
        PHI1[0] = SL[0][0]*CEXP(-I*l_curr*(T[0]-eps_t/2));
        PHI2[0] = SL[1][0]*CEXP(-I*l_curr*(T[0]-eps_t/2));
        n_given = 0;
        count = D_scale - 1;
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
                sh = eps_t_n;
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
                        sum = sum + U[c1][c3]*SL[c3][c2];
                    }
                    TM[c1][c2] = sum;
                    sum = 0;
                }
            }
            for (c1 = 0; c1 < 4; c1++) {
                for (c2 = 0; c2 < 4; c2++)
                    SL[c1][c2] = TM[c1][c2];
            }
            
            if (count == 0){
                count = D_scale - 1;
                n_given++;
                PHI1[n_given] = SL[0][0]*CEXP(-I*l_curr*(T[0]-eps_t/2));
                PHI2[n_given] = SL[1][0]*CEXP(-I*l_curr*(T[0]-eps_t/2));
            }
            else
                count--;
        }
        
        a_vals[neig] = SL[0][0]*CEXP(I*l_curr*(-T[0]+T[1]+eps_t));
        aprime_vals[neig] = scl_factor*(SL[2][0]+I*(T[1]+eps_t-T[0])*SL[0][0])*CEXP(I*l_curr*(-T[0]+T[1]+eps_t));
                //printf("SL00=%1.15e+i%1.15e\n",CREAL(SL[0][0]),CIMAG(SL[0][0]));

       // printf("SL20=%1.15e+i%1.15e\n",CREAL(SL[2][0]),CIMAG(SL[2][0]));
        //printf("aprime=%1.15e+i%1.15e\n",CREAL(aprime_vals[neig]),CIMAG(aprime_vals[neig]));
        
        if (skip_b_flag == 0){
            // Calculation of b assuming a=0
            // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
            // computation point
            REAL error_metric = INFINITY, tmp = INFINITY;
            for (n = 0; n <= D_given; n++){
                tmp = CABS((0.5*LOG(CABS((PHI2[n]/PSI2[n])/(PHI1[n]/PSI1[n])))));
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
        return ret_code;
}
