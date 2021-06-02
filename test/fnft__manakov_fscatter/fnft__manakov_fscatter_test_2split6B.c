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


static INT manakov_fscatter_test_2split6B()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_2SPLIT6B;
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
        eA_6 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB = expm(B*eps_t);
        eB_2 = expm(B*eps_t/2);
        eB_3 = expm(B*eps_t/3);
        eB_4 = expm(B*eps_t/4);
        eB_6 = expm(B*eps_t/6);
        U = (81/40) * eB_6*(eA_6^2*eB_3)^2*eA_6^2*eB_6+...
            -(16/15)*eB_4*eA_6^3*eB_2*eA_6^3*eB_4+...
            (1/24)*eB_2*eA_6^6*eB_2;
        S = U*S;
    end
    S=S*z^(D*6);
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

    COMPLEX result_exact[45] = {    //(for D = 512)
       -0.061878294244769 + 0.000000000000000*I,
  1.603083695263923 - 1.655106374732997*I,
-15.763032123115192 +18.881417962551215*I,
-12.385635306512643 + 0.151352742938106*I,
  0.759315399482681 + 0.767938690936159*I,
 -0.395758726609693 + 0.000000000000000*I,
  0.000000000000030 - 0.000000000000000*I,
 -0.000309183062186 + 0.000148894715312*I,
  4.794742059479877 - 0.000000000008034*I,
  0.000000006670901 - 0.000000020530922*I,
 -0.916267486607223 + 0.000000000000000*I,
  0.000000000000069 - 0.000000000000000*I,
 -0.000715826002682 + 0.000344723634361*I,
 11.100870202926329 - 0.000000000018604*I,
  0.000000015444586 - 0.000000047533549*I,
  0.395758726609694 + 0.000000000000000*I,
 -0.000000000000030 + 0.000000000000000*I,
  0.000309183062194 - 0.000148894715316*I,
 -4.794742059479886 + 0.000000000008034*I,
 -0.000000006670901 + 0.000000020530922*I,
  0.833044082951156 + 0.000000000000000*I,
  1.094821027933055 + 0.260227376436320*I,
 -3.340778841544879 - 0.572170136534140*I,
 -1.104583011059334 - 0.023796734642608*I,
 -0.849371673652169 - 0.467860953630042*I,
 -0.386539242732513 + 0.000000000000000*I,
  0.219531292932126 + 0.602482947620621*I,
 -9.178143707009369 + 0.485411594390734*I,
 -4.872567188649286 - 0.055094613893776*I,
 -0.093429855282566 + 0.277650821903548*I,
  0.916267486607224 + 0.000000000000000*I,
 -0.000000000000069 + 0.000000000000000*I,
  0.000715826002699 - 0.000344723634365*I,
-11.100870202926252 + 0.000000000018604*I,
 -0.000000015444586 + 0.000000047533549*I,
 -0.386539242732512 + 0.000000000000000*I,
  0.219531292932125 + 0.602482947620619*I,
 -9.178143707009387 + 0.485411594390721*I,
 -4.872567188649286 - 0.055094613893776*I,
 -0.093429855282566 + 0.277650821903548*I,
  0.105077622804010 + 0.000000000000000*I,
  1.508262667331915 + 1.394878998294802*I,
-20.625908128500338 + 0.342001882806074*I,
-10.281052295459270 - 0.127556008257117*I,
 -1.025327420192444 + 0.055036759288440*I

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
            q1[i] = 1.5*0.92*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
            q2[i] = 1.5*2.13*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
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
    if (manakov_fscatter_test_2split6B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
