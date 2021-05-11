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

#include "fnft__manakov_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

// TODO: update result_exact and matlab code
INT manakov_scatter_matrix_test_focusing_CF4_2()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    const INT kappa = 1;
    COMPLEX q1[8], q2[8];
    COMPLEX result[18];
    COMPLEX lam[2] = {2, 1+0.5*I};
    // TODO: update matlab code and result
    COMPLEX result_exact[18] = { -0.588197535835346 -     0.613751291319182*I,
         0.0161948933586321 +     0.304375971036892*I,
         -0.134233635547494 +     0.407930552929563*I,
        0.00148492925027899 +      0.31666765174209*I,
         -0.496622965092667 +     0.801077304567259*I,
        -0.0174131387291388 -     0.105152073563462*I,
           0.12239188575443 +     0.402582452260418*I,
         -0.063319103856083 -     0.121344416522301*I,
         -0.579202725562012 +     0.684644449185652*I,
          0.514141690388174 -      1.15666425151093*I,
         0.0328270944528962 +     0.368574677089231*I,
        -0.0151729229597683 +     0.594636336925263*I,
        0.00917285963472292 +     0.391934158711482*I,
          0.258469450645832 +     0.461930501729423*I,
        -0.0703522193258069 -    0.0707967270577841*I,
          0.155348140750472 +     0.668179952413959*I,
          -0.11350259635026 -    0.0724991852912159*I,
          0.117556373723012 +     0.411012176721633*I };
    /* Matlab code to generate result_exact:
clear
eps_t = 0.13;
    kappa = +1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact = []; lam = [2, 1+0.5*i];
    
    %% getting some values needed for the num. method
a1  = 0.25 + sqrt(3)/6;
a2  = 0.25 - sqrt(3)/6;

c1 = 0.5-sqrt(3)/6;
c2 = 0.5+sqrt(3)/6;

% getting the non-equidistant samples
[q1_c1, q2_c1] =  bandlimited_interpolation_CF24(eps_t,...
                                    [q1; q2],c1*eps_t);
[q1_c2, q2_c2] = bandlimited_interpolation_CF24(eps_t,...
                                    [q1; q2],c2*eps_t);

% Getting the interlaced "samples".
Neff = length(q1_c1)+length(q1_c2);
qeff = zeros(2,Neff);
reff = zeros(2,Neff);

qeff(1,2:2:end)=a2*q1_c1+a1*q1_c2;
qeff(2,2:2:end)=a2*q2_c1+a1*q2_c2;
reff(1,2:2:end)=-kappa*(a2*conj(q1_c1)+a1*conj(q1_c2));
reff(2,2:2:end)=-kappa*(a2*conj(q2_c1)+a1*conj(q2_c2));

qeff(1,1:2:end)=a1*q1_c1+a2*q1_c2;
qeff(2,1:2:end)=a1*q2_c1+a2*q2_c2;
reff(1,1:2:end)=-kappa*(a1*conj(q1_c1)+a2*conj(q1_c2));
reff(2,1:2:end)=-kappa*(a1*conj(q2_c1)+a2*conj(q2_c2));
    
    for i = 1:length(lam)
        lambda = lam(i)/2;
        S=eye(3);
        for n=1:2*D         % 2*D because we are effectively doing BO with twice the number of samples
                P = [-1i*lambda,           qeff(1,n),     qeff(2,n);
                reff(1,n),   1i*lambda, 0;
                reff(2,n),   0,      1i*lambda];

            TM = expm(eps_t*P);
            S=TM*S;
        end
        result_exact = [result_exact,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact.'
    
    
    function [q1s, q2s] = bandlimited_interpolation_CF24(eps_t, qn, ts)
% implements bandlimited interpolation (interpolation using the time shift
% property of the Fourier Transform)
% 
% inputs:
% eps_t timestep size
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
    for (i=0; i<D; i++) {
        q1[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));
        q2[i] = 0.21*cos(i+1) + 1.05*I*sin(0.2*(i+1));
    }

    ret_code = manakov_scatter_matrix(D, q1, q2, eps_t, 2, lam, kappa, result, manakov_discretization_CF4_2);
    for (i=0; i<18; i++){
        printf("result[%d] = %f + i%f\n", i, creal(result[i]), cimag(result[i]));
    }

    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(18, result, result_exact) > 10*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (manakov_scatter_matrix_test_focusing_CF4_2() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
