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
#include "fnft__akns_scatter.h"
#include <stdio.h>

/**
 * Returns [S11 S12 S21 S22 S11' S12' S21' S22'] in result 
 * where S = [S11, S12; S21, S22] is the scattering matrix.
 * computed using the chosen scheme.
 * Result should be preallocated with size 8 * K
 * Default scheme should be set as BO.
 */
INT akns_scatter_matrix(const UINT D, COMPLEX const * const q, 
    COMPLEX const * const r, const REAL eps_t,
    const UINT K, COMPLEX const * const lambda,
    COMPLEX * const result, discretization_t discretization)
{
     
    INT ret_code = SUCCESS;
    UINT  neig;
    INT n, c1, c2, c3;
    COMPLEX l, qn, rn, ks, k, sum=0, TM[4][4], ch,
                         chi, sh, u1, ud1, ud2;
    
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
    
    switch (discretization) {
        
        case discretization_BO: // Bofetta-Osborne scheme
            
            for (neig = 0; neig < K; neig++) { // iterate over lambda
                l = lambda[neig];
                
                COMPLEX T[4][4] =
                    { {1,0,0,0}, {0,1,0,0},{0,0,1,0},{0,0,0,1} };
                COMPLEX U[4][4] = {{ 0 }};
               
                for (n = D-1; n >= 0; n--){
                    qn = q[n];
                    rn = r[n];
                    ks = ((q[n]*r[n])-(l*l));
                    k = CSQRT(ks);
                    ch = CCOSH(k*eps_t);
                    chi = ch/ks;
                    sh = CSINH(k*eps_t)/k;
                    u1 = l*sh*I;
                    ud1 = eps_t*l*l*chi*I;
                    ud2 = l*(eps_t*ch-sh)/ks;
                    
		    if (ks != 0){
                    U[0][0] = ch-u1;
                    U[0][1] = qn*sh;
                    U[1][0] = rn*sh;
                    U[1][1] = ch + u1;
                    U[2][0] = ud1-(l*eps_t+I+(l*l*I)/ks)*sh;
                    U[2][1] = -qn*ud2;
                    U[3][0] = -rn*ud2;
                    U[3][1] = -ud1-(l*eps_t-I-(l*l*I)/ks)*sh;
                    U[2][2] = ch-u1;
                    U[2][3] = qn*sh;
                    U[3][2] = rn*sh;
                    U[3][3] = ch+u1;}
                    else{
                    U[0][0] = 1;
                    U[0][1] = 0;
                    U[1][0] = 0;
                    U[1][1] = 1;
                    U[2][0] = 1;
                    U[2][1] = 0;
                    U[3][0] = 0;
                    U[3][1] = 1;
                    U[2][2] = 1;
                    U[2][3] = 0;
                    U[3][2] = 0;
                    U[3][3] = 1;}

                    for (c1 = 0; c1 < 4; c1++) {
                        for (c2 = 0; c2 < 4; c2++) {
                            for (c3 = 0; c3 < 4; c3++) {
                                sum = sum + T[c1][c3]*U[c3][c2];
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
               
		        result[neig*8] = T[0][0];
		        result[neig*8 + 1] = T[0][1];
		        result[neig*8 + 2] = T[1][0];
		        result[neig*8 + 3] = T[1][1];
		        result[neig*8 + 4] = T[2][0];
		        result[neig*8 + 5] = T[2][1];
		        result[neig*8 + 6] = T[3][0];
		        result[neig*8 + 7] = T[3][1];
            }
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(discretization);
    }
    return ret_code;
}
