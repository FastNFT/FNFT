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
 * Lianne de Vries (TU Delft student) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"


static INT manakov_fscatter_test_2split4B()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_2SPLIT4B;
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

for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 20 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_2 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB_2 = expm(B*eps_t/2);
        eB_4 = expm(B*eps_t/4);
        U = ((4/3)*eB_4*(eA_2*eA_2)*(eB_2)*(eA_2*eA_2)*eB_4-...
        (1/3)*(eB_2)*eA_2*eA_2*eA_2*eA_2*(eB_2));
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

    COMPLEX result_exact[45] = {    //(for D = 512)
        0.535312216834846 + 0.000000000000000*I,
  0.966872687808956 - 0.482269048895960*I,
  1.187614580111148 + 0.297208747675138*I,
  1.135602323091052 - 0.351853206210830*I,
  0.872290767293849 + 0.562341955978110*I,
  0.334921105634325 + 0.000000000000000*I,
 -0.000000000000002 + 0.000000000000000*I,
 -0.000000001878687 - 0.000000002355799*I,
  0.000000000067879 + 0.000000000117571*I,
 -0.000000000000042 - 0.000000000000030*I,
  0.775415168479466 + 0.000000000000000*I,
 -0.000000000000005 + 0.000000000000000*I,
 -0.000000004349569 - 0.000000005454187*I,
  0.000000000157156 + 0.000000000272201*I,
 -0.000000000000096 - 0.000000000000070*I,
 -0.334921105634324 + 0.000000000000000*I,
  0.000000000000002 + 0.000000000000000*I,
  0.000000001878687 + 0.000000002355799*I,
 -0.000000000067879 - 0.000000000117570*I,
  0.000000000000042 + 0.000000000000030*I,
  0.926938543333827 + 0.000000000000000*I,
  0.994791492757514 + 0.075825705976718*I,
 -0.183527238589528 + 1.014084640406814*I,
 -0.558569404664136 + 0.856828954589562*I,
  0.386900111208145 + 0.904638093349349*I,
 -0.169153155107542 + 0.000000000000000*I,
 -0.012058826550506 + 0.175552993185681*I,
  0.090278881664109 + 0.090656338173164*I,
 -0.135600904278556 - 0.021291779115683*I,
  0.180316346581180 - 0.107468740081207*I,
 -0.775415168479465 + 0.000000000000000*I,
  0.000000000000005 + 0.000000000000000*I,
  0.000000004349569 + 0.000000005454187*I,
 -0.000000000157156 - 0.000000000272201*I,
  0.000000000000096 + 0.000000000000070*I,
 -0.169153155107542 + 0.000000000000000*I,
 -0.012058826550506 + 0.175552993185681*I,
  0.090278881664109 + 0.090656338173165*I,
 -0.135600904278557 - 0.021291779115683*I,
  0.180316346581180 - 0.107468740081207*I,
  0.608373673501007 + 0.000000000000000*I,
  0.972081195051541 + 0.406443342918829*I,
 -0.013505697059928 + 1.184817042952263*I,
 -0.813945571861527 + 0.816730306484617*I,
  0.726488535916168 + 0.702243020237650*I

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
    if (manakov_fscatter_test_2split4B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
