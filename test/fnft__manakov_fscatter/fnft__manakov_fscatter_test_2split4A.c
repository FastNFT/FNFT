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


static INT manakov_fscatter_test_2split4A()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_2SPLIT4A;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[512], q2[512];   // [D]
    COMPLEX result[45];
    const REAL err_bnd = 3*1000*EPSILON;
    UINT kappa = 1;
    /* %% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
r1 = -kappa*conj(q1);
r2 = -kappa*conj(q2);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);


%%
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 20 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_4 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB = expm(B*eps_t);
        eB_2 = expm(B*eps_t/2);
        U = (4/3)*eA_4*eB_2*eA_4*eA_4*eB_2*eA_4+...
            -(1/3)*eA_4*eA_4*eB*eA_4*eA_4;
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
        0.535312216834838 + 0.000000000000000*I,
 -0.991215593863810 + 0.000000000000108*I,
  2.121295636716989 - 0.506137873379882*I,
  1.953704106189357 + 0.017036261278890*I,
  1.119847103859786 + 0.304166814790438*I,
  0.334921105634321 + 0.000000000000000*I,
  0.367596221582106 - 0.000000000000083*I,
 -0.000000000000006 + 0.000000000000028*I,
 -0.000000000001609 + 0.000000000002788*I,
 -0.000000117569248 - 0.000000361840938*I,
  0.775415168479460 + 0.000000000000000*I,
  0.851065165184651 - 0.000000000000194*I,
 -0.000000000000015 + 0.000000000000065*I,
 -0.000000000003726 + 0.000000000006455*I,
 -0.000000272198367 - 0.000000837740433*I,
 -0.334921105634321 + 0.000000000000000*I,
 -0.367596221582103 + 0.000000000000083*I,
  0.000000000000006 - 0.000000000000028*I,
  0.000000000001610 - 0.000000000002788*I,
  0.000000117569248 + 0.000000361840938*I,
  0.926938543333818 + 0.000000000000000*I,
  0.686927186178219 - 0.000000000000248*I,
 -1.025279981233664 - 0.582074319927440*I,
 -0.577293709014888 - 0.994544915085350*I,
 -0.796151672832402 + 0.637550811775035*I,
 -0.169153155107540 + 0.000000000000000*I,
 -0.724831623304542 + 0.000000000000255*I,
 -0.287807251516808 - 0.343093409919144*I,
 -0.178951739345421 - 0.297550607687590*I,
  0.029786016180405 + 0.115218088801770*I,
 -0.775415168479461 + 0.000000000000000*I,
 -0.851065165184653 + 0.000000000000194*I,
  0.000000000000015 - 0.000000000000065*I,
  0.000000000003726 - 0.000000000006455*I,
  0.000000272198367 + 0.000000837740433*I,
 -0.169153155107540 + 0.000000000000000*I,
 -0.724831623304541 + 0.000000000000255*I,
 -0.287807251516807 - 0.343093409919144*I,
 -0.178951739345421 - 0.297550607687590*I,
  0.029786016180405 + 0.115218088801770*I,
  0.608373673501009 + 0.000000000000000*I,
 -0.678142780042136 + 0.000000000000233*I,
 -1.567305221957693 - 1.228219568604027*I,
 -0.914312179138530 - 1.554919745495158*I,
 -0.740055891696565 + 0.854540175279291*I
    };
        
        i = manakov_fscatter_numel(D, manakov_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
            ret_code = E_INVALID_ARGUMENT(manakov_discretization);
            goto leave_fun;
        }
        

        transfer_matrix = malloc(2*i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<D; i++){
            q1[i] = 0.92*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
            q2[i] = 2.13*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
        }
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, manakov_discretization);  // with kappa =1
    
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
        ret_code = manakov_fscatter(D, q1, q2, 1, eps_t, transfer_matrix, &deg, W_ptr, manakov_discretization); // with kappa = 1

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
        
        for (i=0; i<9; i++){    // replaced 4 by 9
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
    if (manakov_fscatter_test_2split4A() != SUCCESS)
        return EXIT_FAILURE;
    
    printf("SUCCES");
    return EXIT_SUCCESS;
}
