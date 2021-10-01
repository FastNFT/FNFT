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
    COMPLEX q1[512], q2[512];
    COMPLEX result[45];
    const REAL err_bnd = 3*1000*EPSILON;
    UINT kappa = 1;
    /*
%% Matlab file to get the polynomial coefficients for FTES4
% also to generate results for the test file

syms q1 q2 l kappa h z...

A = [-1i*l,0,0;0,1i*l,0;0,0,1i*l];
B = [0,q1,q2;-kappa*conj(q1),0,0;-kappa*conj(q2),0,0];
% Let AE = expm(A*h) and BE = expm(B*h)
% Set z = exp(j lambda h/3) such that
AE = [1/z^3,0,0;0,z^3,0;0,0,z^3];
BE = sym('BE',[3,3]);
% CE=expm((A+B)*h)
% Let CE_approx be approximation of CE after application of splitting
% scheme
% AE_1_3 = expm(A*h/3)
AE_1_3 = [1/z,0,0;0,z,0;0,0,z];
AE_m1_3 = [z, 0, 0; 0, 1/z, 0; 0, 0, 1/z];
BE7_48_ = sym('BE7_48_',[3,3]);
BE3_8_ = sym('BE3_8_',[3,3]);
BEm1_48_ = sym('BEm1_48_',[3,3]);
E1 = sym('E1',[3,3]);
E2 = sym('E2',[3,3]);

CE_approx = E1*(BE7_48_*AE_1_3*BE3_8_*AE_m1_3*BEm1_48_*...
    AE*BEm1_48_*AE_m1_3*BE3_8_*AE_1_3*BE7_48_)*E2;

% Dividing throughout by z^-7 (z^-7 is the lowest power of z found in the
% polynomial)
CE_approx_pos = expand(CE_approx*z^7);

% We can now look at coefficients of individual polynomials. Here c gives
% the coefficients of the element specified by (i,j) in CE_approx_pos(i,j)
% and t gives which power of z the coefficients belong to
[c11,t11] = coeffs(CE_approx_pos(1,1),z);
[c12,t12] = coeffs(CE_approx_pos(1,2),z);
[c13,t13] = coeffs(CE_approx_pos(1,3),z);

[c21,t21] = coeffs(CE_approx_pos(2,1),z);
[c22,t22] = coeffs(CE_approx_pos(2,2),z);
[c23,t23] = coeffs(CE_approx_pos(2,3),z);

[c31,t31] = coeffs(CE_approx_pos(3,1),z);
[c32,t32] = coeffs(CE_approx_pos(3,2),z);
[c33,t33] = coeffs(CE_approx_pos(3,3),z);

%% Computing values for the test file
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
        AE_1_3 = [1/z 0 0; 0 z 0; 0 0 z];
        AE_m1_3 = [z 0 0; 0 1/z 0; 0 0 1/z];
        AE = [1/z^3 0 0; 0 z^3 0; 0 0 z^3];
        B = [0, q1(n), q2(n); -kappa*conj(q1(n)), 0, 0; -kappa*conj(q2(n)), 0, 0];
        BE7_48_ = expm(7*B*eps_t/48);
        BE3_8_ = expm(3*B*eps_t/8);
        BEm1_48_ = expm(-B*eps_t/48);
        U = E1*(BE7_48_*AE_1_3*BE3_8_*AE_m1_3*BEm1_48_*...
        AE*BEm1_48_*AE_m1_3*BE3_8_*AE_1_3*BE7_48_)*E2;
%        U = (BE7_48_*AE_1_3*BE3_8_*AE_m1_3*BEm1_48_*...
%        AE*BEm1_48_*AE_m1_3*BE3_8_*AE_1_3*BE7_48_);
        S = U*S;
    end
    S=S*z^(D*7);
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

#ifdef DEBUG
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
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
#endif
        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{

    if (manakov_fscatter_test_FTES4_suzuki() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
