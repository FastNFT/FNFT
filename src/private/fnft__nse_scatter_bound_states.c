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
* Shrinivas Chimmalgi (TU Delft) 2017.
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
        REAL const *const T,  UINT *trunc_index_ptr, UINT K,
        COMPLEX *bound_states, COMPLEX *a_vals,
        COMPLEX *aprime_vals, COMPLEX *b,
        nse_discretization_t discretization)
{
     
    INT ret_code = SUCCESS;
    REAL norm_left, norm_right;
    UINT i0, i1, neig;
    UINT n, c1, c2, c3;
    COMPLEX l, qn, qnc, ks, k, sum=0, TM[4][4],ch,chi,sh,u1,ud1,ud2;
    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL)
        return E_INVALID_ARGUMENT(eps_t);
    if (trunc_index_ptr == NULL || *trunc_index_ptr>D )
        return E_INVALID_ARGUMENT(trunc_index_ptr);
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
    
    
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    
    if (*trunc_index_ptr == D){
        // Heuristic for where to split the potential: Find the point where the
        // L1 norm of the left half is equal to that of the right half
        i0 = 0;
        i1 = D-1;
        norm_left = 0.0;
        norm_right = 0.0;
        while (i0 < i1) {
            if (norm_left < norm_right) {
                i0++;
                norm_left += eps_t*CABS(q[i0]);
            } else {
                i1--;
                norm_right += eps_t*CABS(q[i1]);
            }
        }
        // i0 is now the index where we will split
	*trunc_index_ptr = i0;
        
    }
    else
        i0 = *trunc_index_ptr;
    
    switch (discretization) {
        
        case nse_discretization_BO: // forward-backward bofetta-osborne scheme
            
            for (neig = 0; neig < K; neig++) { // iterate over bound states
                l = bound_states[neig];
                
                COMPLEX SR[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX SL[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
                COMPLEX U[4][4] = {{0}};
                
                for (n = D-1; n >= i0; n--){
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
                    U[1][1] = ch + u1;
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
                                sum = sum + SR[c1][c3]*U[c3][c2];
                            }
                            TM[c1][c2] = sum;
                            sum = 0;
                        }
                    }
                    for (c1 = 0; c1 < 4; c1++) {
                        for (c2 = 0; c2 < 4; c2++)
                            SR[c1][c2] = TM[c1][c2];
                    }
                }

                // Note that n is unsigned. A normal for (n=i0-1; n>=0; n--)
                // would therefore not stop when n==0.
                if (i0 > 1) {
                    n = i0;
                    do {
                        n--;
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
                                    sum = sum + SL[c1][c3]*U[c3][c2];
                                }
                                TM[c1][c2] = sum;
                                sum = 0;
                            }
                        }
                        for (c1 = 0; c1 < 4; c1++) {
                            for (c2 = 0; c2 < 4; c2++)
                                SL[c1][c2] = TM[c1][c2];
                        }
                    } while (n > 0);
               }

                // Compute the total transfer matrix (TM) from SL and SR
                for (c1 = 0; c1 < 4; c1++) {
                    for (c2 = 0; c2 < 4; c2++) {
                        for (c3 = 0; c3 < 4; c3++) {
                            sum = sum + SR[c1][c3]*SL[c3][c2];
                        }
                        TM[c1][c2] = sum;
                        sum = 0;
                    }
                }
                a_vals[neig] = TM[0][0]*CEXP(-I*l*(T[0]-eps_t/2))*CEXP(I*l*(T[1]+eps_t/2));
                aprime_vals[neig] = (TM[2][0]+I*(T[1]+eps_t/2)*(TM[0][0]+TM[2][2]))*CEXP(-I*l*(T[0]-eps_t/2))*CEXP(I*l*(T[1]+eps_t/2));
                if (CIMAG(l) == 0)
                    b[neig] = TM[1][0]*CEXP(-I*l*(T[0]-eps_t/2))*CEXP(-I*l*(T[1]+eps_t/2));
                else
                    b[neig] = SL[1][0]/SR[0][0];
            }
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(discretization);
    }
    return ret_code;
}
