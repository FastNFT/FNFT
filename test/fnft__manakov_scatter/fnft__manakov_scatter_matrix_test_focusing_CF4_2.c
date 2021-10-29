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

#include "fnft__manakov_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT manakov_scatter_matrix_test_focusing_CF4_2()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    const INT kappa = 1;
    COMPLEX q1[8], q2[8];
    COMPLEX result[18];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[18] = { 0.235870243814537 -      0.67240741647598*I,
        -0.0160013006792082 +     0.343757091014225*I,
         -0.139417083479491 +     0.595291014495476*I,
         0.0316718715200687 +      0.35616288784048*I,
          0.454838382307009 +     0.806518284627255*I,
        -0.0876930926159252 -    0.0842918378566583*I,
          0.129153188764432 +     0.589659086305378*I,
         -0.130670744901392 -    0.0845605092414831*I,
          0.286900348409401 +     0.727381937193787*I,
           0.75926006766599 -     0.496251446738425*I,
        0.00411038462756639 +     0.344846116857216*I,
        -0.0496904558168458 +     0.631319952907485*I,
        0.00847708904602566 +     0.380366077602955*I,
          0.610198362265694 +     0.346117553467807*I,
        -0.0993807831886292 -    0.0536042664752771*I,
         0.0829944502053134 +     0.700732853349535*I,
         -0.138508583171897 -    0.0489598976143097*I,
          0.434538899204319 +      0.31191061443021*I };
    /* Matlab code to generate result_exact:
clear lam lambda result_exact_CF
clear
eps_t = 0.13;
    kappa = +1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact = []; lam = [2, 1+0.5*i]./2;       %Dividing by 2 because CF4_2 is essentially BO but with 2x the amnt of samples and lambda/2 to compensate. Getting the samples is taken care of outside manakov_scatter
    for i = 1:length(lam)
        lambda = lam(i);
        S=eye(3);
        for n=1:D
            P = [-1i*lambda,           q1(n),     q2(n);
                -kappa*conj(q1(n)),   1i*lambda, 0;
                -kappa*conj(q2(n)),   0,      1i*lambda];

            TM = expm(eps_t*P);
            TM_reshaped = reshape(TM.',[1,9]);
            S=TM*S;
            S_reshaped = reshape(S.',[1,9]);
        end
        result_exact = [result_exact,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact.'
    */
    for (i=0; i<D; i++) {
        q1[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));
        q2[i] = 0.21*cos(i+1) + 1.05*I*sin(0.2*(i+1));
    }

    ret_code = manakov_scatter_matrix(D, q1, q2, eps_t, 2, lam, kappa, result, manakov_discretization_CF4_2);
    #ifdef DEBUG
    printf("result\n");
        for (UINT j = 0; j<18; j++){
                printf("%f + i%f\n", creal(result[j]), cimag(result[j]));
            }
    #endif  // DEBUG

    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(18, result, result_exact) > 1000*EPSILON){
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

INT main()
{
    if (manakov_scatter_matrix_test_focusing_CF4_2() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
