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

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 * Default scheme should be set as BO.
 */
INT nse_scatter_bound_states(const UINT D, COMPLEX const *const q,
        REAL const *const T, UINT K,
        COMPLEX *bound_states, COMPLEX *a_vals,
        COMPLEX *aprime_vals, COMPLEX *b,
        nse_discretization_t discretization)
{
    
    INT ret_code = SUCCESS;
    UINT neig;
    UINT c1, c2, c3;
    UINT n;
    COMPLEX l, qn, qnc, ks, k, sum=0, TM[4][4],ch,chi,sh,u1,ud1,ud2;
    REAL eps_t_n = 0, eps_t = 0;
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
    
    
    PHI1 = malloc((D+1) * sizeof(COMPLEX));
    PHI2 = malloc((D+1) * sizeof(COMPLEX));
    PSI1 = malloc((D+1) * sizeof(COMPLEX));
    PSI2 = malloc((D+1) * sizeof(COMPLEX));
    if (PHI1 == NULL || PHI2 == NULL || PSI1 == NULL || PSI2 == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }
    switch (discretization) {
        
        case nse_discretization_BO: // forward-backward bofetta-osborne scheme
            
            
            eps_t = (T[1] - T[0])/(D - 1);
            eps_t_n = -eps_t;
            
            for (neig = 0; neig < K; neig++) { // iterate over bound states
                l = bound_states[neig];
                
                COMPLEX SR[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX SL[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX U[4][4] = {{0}};
                
                PSI1[D] = SR[0][1]*CEXP(I*l*(T[1]+eps_t/2));
                PSI2[D] = SR[1][1]*CEXP(I*l*(T[1]+eps_t/2));
                n = D;
                do{
                    n--;
                    qn = q[n];
                    qnc = CONJ(qn);
                    ks = (-(CABS(qn)*CABS(qn))-(l*l)); //TODO:ks==0 case
                    k = CSQRT(ks);
                    ch = CCOSH(k*eps_t_n);
                    chi = ch/ks;
                    sh = CSINH(k*eps_t_n)/k;
                    u1 = l*sh*I;
                    ud1 = eps_t_n*l*l*chi*I;
                    ud2 = l*(eps_t_n*ch-sh)/ks;
                    
                    U[0][0] = ch-u1;
                    U[0][1] = qn*sh;
                    U[1][0] = -qnc*sh;
                    U[1][1] = ch + u1;
                    U[2][0] = ud1-(l*eps_t_n+I+(l*l*I)/ks)*sh;
                    U[2][1] = -qn*ud2;
                    U[3][0] = qnc*ud2;
                    U[3][1] = -ud1-(l*eps_t_n-I-(l*l*I)/ks)*sh;
                    U[2][2] = ch-u1;
                    U[2][3] = qn*sh;
                    U[3][2] = -qnc*sh;
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
                    
                    PSI1[n] = SR[0][1]*CEXP(I*l*(T[1]+eps_t/2));
                    PSI2[n] = SR[1][1]*CEXP(I*l*(T[1]+eps_t/2));
                    
                } while (n > 0);
                
                PHI1[0] = SL[0][0]*CEXP(-I*l*(T[0]-eps_t/2));
                PHI2[0] = SL[1][0]*CEXP(-I*l*(T[0]-eps_t/2));
                for (n = 0; n < D; n++){
                    qn = q[n];
                    qnc = CONJ(qn);
                    ks = (-(CABS(qn)*CABS(qn))-(l*l));
                    k = CSQRT(ks);
                    ch = CCOSH(k*eps_t);
                    chi = ch/ks;
                    sh = CSINH(k*eps_t)/k;
                    u1 = l*sh*I;
                    ud1 = eps_t*l*l*chi*I;
                    ud2 = l*(eps_t*ch-sh)/ks;
                    U[0][0] = ch-u1;
                    U[0][1] = qn*sh;
                    U[1][0] = -qnc*sh;
                    U[1][1] = ch+u1;
                    U[2][0] = ud1-(l*eps_t+I+(l*l*I)/ks)*sh;
                    U[2][1] = -qn*ud2;
                    U[3][0] = qnc*ud2;
                    U[3][1] = -ud1-(l*eps_t-I-(l*l*I)/ks)*sh;
                    U[2][2] = ch-u1;
                    U[2][3] = qn*sh;
                    U[3][2] = -qnc*sh;
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
                    PHI1[n+1] = SL[0][0]*CEXP(-I*l*(T[0]-eps_t/2));
                    PHI2[n+1] = SL[1][0]*CEXP(-I*l*(T[0]-eps_t/2));
                }
                
                a_vals[neig] = SL[0][0]*CEXP(-I*l*(T[0]-eps_t/2))*CEXP(I*l*(T[1]+eps_t/2));
                aprime_vals[neig] = (SL[2][0]+I*(T[1]+eps_t-T[0])*SL[0][0])*CEXP(I*l*(-T[0]+T[1]+eps_t));
                
                // Calculation of b assuming a=0
                // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
                // computation point
                REAL error_metric = INFINITY, tmp = INFINITY;
                for (n = 0; n <= D; n++){
                    tmp = CABS((0.5*LOG(CABS((PHI2[n]/PSI2[n])/(PHI1[n]/PSI1[n])))));
                    if (tmp < error_metric){
                        b[neig] = PHI1[n]/PSI1[n];
                        error_metric = tmp;
                    }
                }
            }
            break;
            
        case nse_discretization_CF4_2: // forward-backward commutator-free fourth-order
            
            // The routine requires interleaved non-equispaced samples
            eps_t = (T[1] - T[0])/(D/2 - 1);
            eps_t_n = -eps_t;
            COMPLEX l_curr = 0;
            for (neig = 0; neig < K; neig++) { // iterate over bound states
                l_curr = bound_states[neig];
                l = l_curr/2;
                
                COMPLEX SR[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX SL[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX U[4][4] = {{0}};
                
                PSI1[D] = SR[0][1]*CEXP(I*l*(T[1]+eps_t/2));
                PSI2[D] = SR[1][1]*CEXP(I*l*(T[1]+eps_t/2));
                n = D;
                do{
                    n--;
                    qn = q[n];
                    qnc = CONJ(qn);
                    ks = (-(CABS(qn)*CABS(qn))-(l*l)); //TODO:ks==0 case
                    k = CSQRT(ks);
                    ch = CCOSH(k*eps_t_n);
                    chi = ch/ks;
                    sh = CSINH(k*eps_t_n)/k;
                    u1 = l*sh*I;
                    ud1 = eps_t_n*l*l*chi*I;
                    ud2 = l*(eps_t_n*ch-sh)/ks;
                    
                    U[0][0] = ch-u1;
                    U[0][1] = qn*sh;
                    U[1][0] = -qnc*sh;
                    U[1][1] = ch + u1;
                    U[2][0] = ud1-(l*eps_t_n+I+(l*l*I)/ks)*sh;
                    U[2][1] = -qn*ud2;
                    U[3][0] = qnc*ud2;
                    U[3][1] = -ud1-(l*eps_t_n-I-(l*l*I)/ks)*sh;
                    U[2][2] = ch-u1;
                    U[2][3] = qn*sh;
                    U[3][2] = -qnc*sh;
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
                    
                    PSI1[n] = SR[0][1]*CEXP(I*l*(T[1]+eps_t/2));
                    PSI2[n] = SR[1][1]*CEXP(I*l*(T[1]+eps_t/2));
                    
                } while (n > 0);
                
                PHI1[0] = SL[0][0]*CEXP(-I*l*(T[0]-eps_t/2));
                PHI2[0] = SL[1][0]*CEXP(-I*l*(T[0]-eps_t/2));
                for (n = 0; n < D; n++){
                    qn = q[n];
                    qnc = CONJ(qn);
                    ks = (-(CABS(qn)*CABS(qn))-(l*l));
                    k = CSQRT(ks);
                    ch = CCOSH(k*eps_t);
                    chi = ch/ks;
                    sh = CSINH(k*eps_t)/k;
                    u1 = l*sh*I;
                    ud1 = eps_t*l*l*chi*I;
                    ud2 = l*(eps_t*ch-sh)/ks;
                    U[0][0] = ch-u1;
                    U[0][1] = qn*sh;
                    U[1][0] = -qnc*sh;
                    U[1][1] = ch+u1;
                    U[2][0] = ud1-(l*eps_t+I+(l*l*I)/ks)*sh;
                    U[2][1] = -qn*ud2;
                    U[3][0] = qnc*ud2;
                    U[3][1] = -ud1-(l*eps_t-I-(l*l*I)/ks)*sh;
                    U[2][2] = ch-u1;
                    U[2][3] = qn*sh;
                    U[3][2] = -qnc*sh;
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
                    PHI1[n+1] = SL[0][0]*CEXP(-I*l*(T[0]-eps_t/2));
                    PHI2[n+1] = SL[1][0]*CEXP(-I*l*(T[0]-eps_t/2));
                }
                
                a_vals[neig] = SL[0][0]*CEXP(I*l_curr*(-T[0]+T[1]+eps_t));
                aprime_vals[neig] = 0.5*(SL[2][0]+I*(T[1]+eps_t-T[0])*SL[0][0])*CEXP(I*l_curr*(-T[0]+T[1]+eps_t));
                
                // Calculation of b assuming a=0
                // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
                // computation point
                REAL error_metric = INFINITY, tmp = INFINITY;
                for (n = 0; n <= D; n++){
                    tmp = CABS((0.5*LOG(CABS((PHI2[n]/PSI2[n])/(PHI1[n]/PSI1[n])))));
                    if (tmp < error_metric){
                        b[neig] = PHI1[n]/PSI1[n];
                        error_metric = tmp;
                    }
                }
            }
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(discretization);
    }
    
    release_mem:
        free(PHI1);
        free(PSI1);
        free(PHI2);
        free(PSI2);
        return ret_code;
}
