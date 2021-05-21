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

// For now: using this file to check if misc_matrix_mult works as I think it does


static INT manakov_fscatter_test_FTES4_suzuki()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_FTES4_suzuki;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[512], q2[512];   // [D]
    COMPLEX result[45];
    const REAL err_bnd = 3*1000*EPSILON;
    UINT kappa = 1;
    // TODO: adjust matlab and result
    /* %% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
c1 = 0.5-sqrt(3)/6;
c2 = 0.5+sqrt(3)/6;
a1  = 0.25 + sqrt(3)/6;
a2  = 0.25 - sqrt(3)/6;

% getting the non-equidistant samples
[q1_c1, q2_c1] =  bandlimited_interpolation(eps_t,...
                                    [q1; q2],c1*eps_t);
[q1_c2, q2_c2] = bandlimited_interpolation(eps_t,...
                                    [q1; q2],c2*eps_t);

%%
% Interlacing the samples:
qeff = zeros(2,2*D);
reff = zeros(2,2*D);

qeff(1,2:2:end)=a2*q1_c1+a1*q1_c2;
qeff(2,2:2:end)=a2*q2_c1+a1*q2_c2;
reff(1,2:2:end)=-kappa*(a2*conj(q1_c1)+a1*conj(q1_c2));
reff(2,2:2:end)=-kappa*(a2*conj(q2_c1)+a1*conj(q2_c2));

qeff(1,1:2:end)=a1*q1_c1+a2*q1_c2;
qeff(2,1:2:end)=a1*q2_c1+a2*q2_c2;
reff(1,1:2:end)=-kappa*(a1*conj(q1_c1)+a2*conj(q1_c2));
reff(2,1:2:end)=-kappa*(a1*conj(q2_c1)+a2*conj(q2_c2));

for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:2*D
        % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_4 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, qeff(1,n), qeff(2,n);...
             reff(1,n), 0, 0;
             reff(2,n), 0, 0];
        eB = expm(B*eps_t);
        eB_1_2 = expm(B*eps_t/2);
        U = (4/3)*eA_4*eB_1_2*eA_4*eA_4*eB_1_2*eA_4+...
            -(1/3)*eA_4*eA_4*eB*eA_4*eA_4;
        S = U*S;
    end
    S=S*z^(2*D*4);
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

%%
function [q1s, q2s] = bandlimited_interpolation(eps_t, qn, ts)
% implements bandlimited interpolation (interpolation using the time shift
% property of the Fourier Transform)
% 
% inputs:
% eps_t sampling period
% qn    values of the samples of q taken at the points in vector tn
% ts    desired shift of the samples
% 
% output: qs    the shifted vector of samples

Qn = [fft(qn(1,:)); fft(qn(2,:))];
N = length(qn(1,:));
Np = floor(N/2);
Nn = -floor((N-1)/2);
Qn = [Qn(1,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t));...
      Qn(2,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t))];
q1s = ifft(Qn(1,:));
q2s = ifft(Qn(2,:));
end
*/

    COMPLEX result_exact[45] = {    //(for D = 512)
        0.535312216834849 + 0.000000000000000*I,
  0.965987341039084 - 0.258589359704154*I,
 -0.403732646652832 + 0.914876241605129*I,
 -0.136420616209859 + 0.618301251926966*I,
 -0.038523111583892 + 0.999257709439486*I,
  0.334921105634327 + 0.000000000000000*I,
 -0.000000042911225 + 0.000000000000000*I,
 -0.000471269039502 + 0.000000000000000*I,
  0.153454862586720 - 0.265791618669811*I,
  0.000000000001424 + 0.000000000001034*I,
  0.775415168479470 + 0.000000000000000*I,
 -0.000000099348816 + 0.000000000000000*I,
 -0.001091090276238 + 0.000000000000000*I,
  0.355281366640996 - 0.615365378007285*I,
  0.000000000003296 + 0.000000000002395*I,
 -0.334921105634328 + 0.000000000000000*I,
  0.000000042911226 + 0.000000000000000*I,
  0.000471269039502 - 0.000000000000000*I,
 -0.153454862586720 + 0.265791618669810*I,
 -0.000000000001424 - 0.000000000001034*I,
  0.926938543333835 + 0.000000000000000*I,
  0.994652292358891 + 0.040657224016971*I,
 -0.251012304236169 - 0.965486169388773*I,
 -0.494851489855136 + 0.797045044529699*I,
  0.990322091347419 - 0.054310171969179*I,
 -0.169153155107543 + 0.000000000000000*I,
 -0.012381105734287 + 0.094130312127464*I,
 -0.065963715973562 + 0.021859687117875*I,
  0.011919920222042 - 0.159704527408744*I,
 -0.022406462424323 - 0.125739854668027*I,
 -0.775415168479471 + 0.000000000000000*I,
  0.000000099348815 - 0.000000000000000*I,
  0.001091090276238 - 0.000000000000000*I,
 -0.355281366640992 + 0.615365378007278*I,
 -0.000000000003296 - 0.000000000002395*I,
 -0.169153155107542 + 0.000000000000000*I,
 -0.012381105734287 + 0.094130312127464*I,
 -0.065963715973562 + 0.021859687117875*I,
  0.011919920222042 - 0.159704527408743*I,
 -0.022406462424323 - 0.125739854668027*I,
  0.608373673501015 + 0.000000000000000*I,
  0.971335048680390 + 0.217932135685897*I,
 -0.375241276373556 - 0.924317984397696*I,
 -0.472402793396961 + 0.496274704459300*I,
  0.948124168517778 - 0.291115098307387*I
    };
        
        i = manakov_fscatter_numel(D, manakov_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
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
        misc_print_buf(4,q1,"q1");
        misc_print_buf(4,q2,"q2");
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, manakov_discretization);  // with kappa =1
        misc_print_buf(540,transfer_matrix,"TM_C");
    
        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }

        
        printf("deg=%d\n",deg);
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
        misc_print_buf(45,result,"result_C");
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

    // Trying out misc_matrix_mult
    COMPLEX A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    COMPLEX B[9] = {7, 8, 9, 1, 2, 3, 4, 5, 6};
    COMPLEX C[9];

    fnft__misc_matrix_mult(3,3,3,&A[0],&B[0],&C[0]);
    for (UINT i=0; i<9; i++){
        printf("%f+i%f\n",creal(C[i]),cimag(C[i]));
    }




    if (manakov_fscatter_test_FTES4_suzuki() != SUCCESS)
        return EXIT_FAILURE;
    
    printf("SUCCES");
    return EXIT_SUCCESS;
}
