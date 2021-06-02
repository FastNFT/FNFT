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
 * Lianne de Vries (TU Delft) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"


static INT manakov_fscatter_test_FTES4_4B()
{
    UINT i, j, D = 512, deg, nz = 5; 
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_FTES4_4B;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[512], q2[512];   // [D]
    COMPLEX result[45];
    const REAL err_bnd = 2*1000*EPSILON;
    UINT kappa = 1;
    /* %% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        if n == 1
    Q1 = [0                                 q1(2)-q1(1) q2(2)-q2(1);
          -kappa*(conj(q1(2))-conj(q1(1)))  0           0;
          -kappa*(conj(q2(2))-conj(q2(1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(2)-2*q1(1)+q1(1)       q2(2)-2*q2(1)+q2(1);
          -kappa*(conj(q1(2))-2*conj(q1(1))+conj(q1(1)))    0 0;
          -kappa*(conj(q2(2))-2*conj(q2(1))+conj(q2(1)))    0 0]/(eps_t^2);
        elseif n == D
    Q1 = [0                                 q1(n)-q1(n-1) q2(n)-q2(n-1);
          -kappa*(conj(q1(n))-conj(q1(n-1)))  0           0;
          -kappa*(conj(q2(n))-conj(q2(n-1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(n)-2*q1(n)+q1(n-1)       q2(n)-2*q2(n)+q2(n-1);
          -kappa*(conj(q1(n))-2*conj(q1(n))+conj(q1(n-1)))    0 0;
          -kappa*(conj(q2(n))-2*conj(q2(n))+conj(q2(n-1)))    0 0]/(eps_t^2);
        else
    Q1 = [0                                 q1(n+1)-q1(n-1) q2(n+1)-q2(n-1);
          -kappa*(conj(q1(n+1))-conj(q1(n-1)))  0           0;
          -kappa*(conj(q2(n+1))-conj(q2(n-1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(n+1)-2*q1(n)+q1(n-1)       q2(n+1)-2*q2(n)+q2(n-1);
          -kappa*(conj(q1(n+1))-2*conj(q1(n))+conj(q1(n-1)))    0 0;
          -kappa*(conj(q2(n+1))-2*conj(q2(n))+conj(q2(n-1)))    0 0]/(eps_t^2);
        end
  
E1 = expm(eps_t^2*Q1/12 + eps_t^3*Q2/48);
E2 = expm(-eps_t^2*Q1/12 + eps_t^3*Q2/48);
        eA_4 = [z^-1 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n); -kappa*conj(q1(n)), 0, 0; -kappa*conj(q2(n)), 0, 0];
        eB_1_2 = expm(B*eps_t/2);
        eB_1_4 = expm(B*eps_t/4);
        U = E1*((4/3)*eB_1_4*eA_4*eA_4*eB_1_2*eA_4*eA_4*eB_1_4+...
            -(1/3)*eB_1_2*eA_4*eA_4*eA_4*eA_4*eB_1_2)*E2;
        S = U*S;
    end
    S=S*z^(D*4);
    result_exact(i) = S(1,1);
    result_exact(5+i) = S(1,2);
    result_exact(10+i) = S(1,3);
    result_exact(15+i) = S(2,1);
    result_exact(20+i) = S(2,2);
    result_exact(25+i) = S(2,3);
    result_exact(30+i) = S(3,1);
    result_exact(35+i) = S(3,2);
    result_exact(40+i) = S(3,3);
end
*/

    COMPLEX result_exact[45] = {    
        0.535312216834847 + 0.000000000000000*I,
 -0.991215593863820 + 0.000000000000102*I,
  1.035523329343331 - 0.432580480596043*I,
  1.135536773929817 + 0.352064695736202*I,
  1.275152492352463 + 0.156701914027315*I,
  0.334921105634324 + 0.000000000000000*I,
 -0.367596221582111 + 0.000000000000083*I,
  0.000000000000053 - 0.000000000000233*I,
  0.000000000058863 - 0.000000000101954*I,
 -0.000000214024617 - 0.000000658700039*I,
  0.775415168479466 + 0.000000000000000*I,
 -0.851065165184663 + 0.000000000000192*I,
  0.000000000000123 - 0.000000000000540*I,
  0.000000000136281 - 0.000000000236047*I,
 -0.000000495513515 - 0.000001525033787*I,
 -0.334921105634325 + 0.000000000000000*I,
  0.367596221582107 - 0.000000000000083*I,
 -0.000000000000053 + 0.000000000000233*I,
 -0.000000000058863 + 0.000000000101954*I,
  0.000000214024617 + 0.000000658700039*I,
  0.926938543333828 + 0.000000000000000*I,
  0.686927186178230 - 0.000000000000253*I,
 -0.876491107795999 - 0.497584963899090*I,
 -0.558593048573332 - 0.856803403316768*I,
 -0.829534481336236 + 0.633146036016526*I,
 -0.169153155107543 + 0.000000000000000*I,
 -0.724831623304554 + 0.000000000000258*I,
  0.056671335898622 - 0.147482183462203*I,
 -0.135655645062242 + 0.021350935863545*I,
 -0.047502442638281 + 0.105020075360860*I,
 -0.775415168479466 + 0.000000000000000*I,
  0.851065165184666 - 0.000000000000192*I,
 -0.000000000000123 + 0.000000000000540*I,
 -0.000000000136282 + 0.000000000236047*I,
  0.000000495513515 + 0.000001525033787*I,
 -0.169153155107542 + 0.000000000000000*I,
 -0.724831623304554 + 0.000000000000258*I,
  0.056671335898622 - 0.147482183462203*I,
 -0.135655645062242 + 0.021350935863545*I,
 -0.047502442638281 + 0.105020075360860*I,
  0.608373673501010 + 0.000000000000000*I,
 -0.678142780042149 + 0.000000000000234*I,
 -0.769762405441516 - 0.775337055176522*I,
 -0.814072308678563 - 0.816593345751752*I,
 -0.918995475700696 + 0.830929557204138*I
    };
        
        i = manakov_fscatter_numel(D, manakov_discretization);
        if (i == 0) { 
            ret_code = E_INVALID_ARGUMENT(manakov_discretization);
            goto leave_fun;
        }
        

        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<D; i++){
            q1[i] = 0.92*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
            q2[i] = 2.13*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
        }
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, manakov_discretization); 
    
        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }

        
        for (i=0; i<9; i++){
            for (j=0; j<nz; j++){
                result[i*nz+j] = z[j];
            }
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }

#define DEBUG
#ifdef DEBUG
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }
#endif  // DEBUG

        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        // with normalization
        W_ptr = &W;
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, W_ptr, manakov_discretization);

        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }
        if (W == 0){
            return E_TEST_FAILED;
            goto leave_fun;
        }
        scl = POW(2.0, W);
        for (i=0; i<9*(deg+1); i++)
            transfer_matrix[i] *= scl;
        
        for (i=0; i<9; i++){   
            for (j=0; j<nz; j++)
                result[i*nz+j] = z[j];
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
        
#ifdef DEBUG
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }
#endif
        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    if (manakov_fscatter_test_FTES4_4B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
