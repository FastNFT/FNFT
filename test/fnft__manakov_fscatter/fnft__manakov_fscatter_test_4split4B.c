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


static INT manakov_fscatter_test_4split4B()
{
    UINT i, j, D = 512, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_4SPLIT4B;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[512], q2[512];   // [D]
    COMPLEX result[45];
    const REAL err_bnd = 100*EPSILON;
    UINT kappa = 1;
    /* matlab code to compute results:
    %% Exact answer for test file:
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
        eB_1_4 = expm(B*eps_t/4);
        U = ((4/3)*eB_1_4*eA_4*eA_4*eB_1_2*eA_4*eA_4*eB_1_4)+...
            -((1/3)*eB_1_2*eA_4*eA_4*eA_4*eA_4*eB_1_2);
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

    COMPLEX result_exact[45] = {    //(for D = 4)
          0.535312216834830 + 0.000000000000000*I,
 -0.876100554400426 + 0.000000000000181*I,
  1.038800200364441 - 0.209527447358909*I,
  1.078710463977776 + 0.164380951016611*I,
  1.133018291154739 + 0.069406226135903*I,
  0.334921105634316 + 0.000000000000000*I,
 -0.305961558188488 + 0.000000000000139*I,
  0.000002038933887 + 0.000000000000000*I,
  0.000000000000224 - 0.000000000000388*I,
 -0.000000000000059 - 0.000000000000182*I,
  0.775415168479455 + 0.000000000000000*I,
 -0.708367520588573 + 0.000000000000321*I,
  0.000004720575195 + 0.000000000000000*I,
  0.000000000000518 - 0.000000000000898*I,
 -0.000000000000137 - 0.000000000000422*I,
 -0.334921105634316 + 0.000000000000000*I,
  0.305961558188484 - 0.000000000000139*I,
 -0.000001271254485 - 0.000001594102704*I,
  0.000000000000447 + 0.000000000000000*I,
 -0.000000000000155 - 0.000000000000113*I,
  0.926938543333805 + 0.000000000000000*I,
  0.705026376154915 - 0.000000000000517*I,
  0.601537200914534 + 0.807140851012290*I,
 -0.483805173381847 + 0.889665382178414*I,
  0.305101361612649 - 0.974319130930809*I,
 -0.169153155107540 + 0.000000000000000*I,
 -0.682928063901840 + 0.000000000000489*I,
 -0.050825043491043 + 0.058596690217650*I,
  0.037494544228042 + 0.054731689102213*I,
 -0.009065541069525 - 0.053858009970816*I,
 -0.775415168479451 + 0.000000000000000*I,
  0.708367520588559 - 0.000000000000322*I,
 -0.000002943230493 - 0.000003690694303*I,
  0.000000000001036 + 0.000000000000000*I,
 -0.000000000000359 - 0.000000000000260*I,
 -0.169153155107539 + 0.000000000000000*I,
 -0.682928063901840 + 0.000000000000489*I,
 -0.050825043491043 + 0.058596690217650*I,
  0.037494544228042 + 0.054731689102213*I,
 -0.009065541069525 - 0.053858009970816*I,
  0.608373673501002 + 0.000000000000000*I,
 -0.581126930555507 + 0.000000000000403*I,
  0.505818777254755 + 0.917495558732356*I,
 -0.413191979120976 + 0.992741162251760*I,
  0.288028296029328 - 1.075749517640996*I
    };
        
        i = manakov_fscatter_numel(2*D, akns_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto leave_fun;
        }
        

        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<D; i++){
//            q1[i] = (0.41*COS(2*i+1) + 0.59*I*SIN(0.28*(2*i+1)))*50;
  //          q2[i] = (0.41*COS(2*(i+1)) + 0.59*I*SIN(0.28*(2*(i+1))))*50;
//            q1[i] = (0.41*cos(0.012*(2*i+1))+0.59*I*sin(0.08*(2*i+1)))*50;
  //          q2[i] = (0.41*cos(0.012*(2*(i+1)))+0.59*I*sin(0.08*(2*(i+1))))*50;
            q1[i] = 0.92*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
            q2[i] = 2.13*1/COSH(i*eps_t-(D-1)*eps_t*0.5);
        }
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, akns_discretization);  // with kappa =1
      misc_print_buf(deg+1,transfer_matrix,"p1");
    
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

        if (misc_rel_err(9*nz, result, result_exact) > 100*err_bnd)
            return E_TEST_FAILED;
        
        
        // with normalization
        W_ptr = &W;
        ret_code = manakov_fscatter(D, q1, q2, 1, eps_t, transfer_matrix, &deg, W_ptr, akns_discretization); // with kappa = 1

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
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),err_bnd);
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
    if (manakov_fscatter_test_4split4B() != SUCCESS)
        return EXIT_FAILURE;
    
    printf("SUCCES");
    return EXIT_SUCCESS;
}
