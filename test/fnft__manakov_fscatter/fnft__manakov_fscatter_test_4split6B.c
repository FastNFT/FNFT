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


static INT manakov_fscatter_test_4split6B()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_4SPLIT6B;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[512], q2[512];   // [D]
    COMPLEX result[45];
    const REAL err_bnd = 4*1000*EPSILON;
    UINT kappa = 1;
    /* %% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 1.5*0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 1.5*2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
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
        % Eq. 22 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_6 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, qeff(1,n), qeff(2,n);...
             reff(1,n), 0, 0;
             reff(2,n), 0, 0];
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
    S=S*z^(2*D*6);
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
       -0.426881661014678 + 0.000000000000000*I,
  1.589056813406154 - 1.423091918606036*I,
-16.628859139917683 +14.026459728346595*I,
 -2.278074294410628 + 0.013845125922587*I,
  0.809390460113850 + 0.695158762739755*I,
  0.358574719043732 + 0.000000000000000*I,
 -0.023629490528572 + 0.000000000000016*I,
 -0.000000003609930 - 0.000000000000005*I,
 -5.964082113257203 + 0.000000000019956*I,
  0.000000000000014 - 0.000000000000043*I,
  0.830178425612130 + 0.000000000000000*I,
 -0.054707407419411 + 0.000000000000037*I,
 -0.000000008357790 - 0.000000000000009*I,
-13.808146631780529 + 0.000000000046220*I,
  0.000000000000033 - 0.000000000000099*I,
 -0.358574719043732 + 0.000000000000000*I,
  0.023629490528572 - 0.000000000000016*I,
 -0.000000000803307 - 0.000000003519426*I,
  5.964082113257278 - 0.000000000020019*I,
  0.000000000000036 - 0.000000000000027*I,
  0.775655706001278 + 0.000000000000000*I,
  1.092615623663557 + 0.223748444244153*I,
 -1.755797952576416 + 2.218045100514900*I,
  0.484598279347938 - 0.002176827335564*I,
  0.403704377612041 + 0.888779488776262*I,
 -0.519405811105474 + 0.000000000000000*I,
  0.214425302612647 + 0.518026289394804*I,
 -3.549869619194624 + 7.392426648959383*I,
 -1.193267027163931 - 0.005039828491079*I,
  0.219221876407042 - 0.144184857191663*I,
 -0.830178425612129 + 0.000000000000000*I,
  0.054707407419411 - 0.000000000000037*I,
 -0.000000001859816 - 0.000000008148235*I,
 13.808146631780396 - 0.000000000046365*I,
  0.000000000000085 - 0.000000000000061*I,
 -0.519405811105477 + 0.000000000000000*I,
  0.214425302612647 + 0.518026289394804*I,
 -3.549869619194631 + 7.392426648959386*I,
 -1.193267027163929 - 0.005039828491079*I,
  0.219221876407042 - 0.144184857191663*I,
 -0.202537367016067 + 0.000000000000000*I,
  1.496441189744397 + 1.199343474358212*I,
 -8.441240813179398 +16.140146829431195*I,
 -1.762672573759059 - 0.011668298577827*I,
  0.816563295187157 + 0.617237227362817*I
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
            q1[i] = 2*0.92*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
            q2[i] = 2*2.13*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
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
                printf("w = %d\n",W);
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
    if (manakov_fscatter_test_4split6B() != SUCCESS)
        return EXIT_FAILURE;
    
    printf("SUCCES");
    return EXIT_SUCCESS;
}
