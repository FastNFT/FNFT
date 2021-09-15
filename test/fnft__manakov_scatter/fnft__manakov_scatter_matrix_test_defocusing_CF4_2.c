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

INT manakov_scatter_matrix_test_defocusing_CF4_2()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    const INT kappa = -1;
    COMPLEX q1[8], q2[8];
    COMPLEX result[18];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[18] = { 0.819736220396989 -      1.06971941301528*I,
        -0.0457063915474174 +     0.450765460617743*I,
         -0.167493111329713 +     0.763501523610377*I,
        -0.0289195047857921 -     0.436833992421533*I,
            0.5651496162473 +     0.924024643501244*I,
          0.100514910636476 +    0.0913312481556059*I,
         -0.178713039247989 -     0.769849640204658*I,
           0.15234015174402 +    0.0941032536034179*I,
           0.76043089059641 +      1.00711902741191*I,
           1.54431774593891 -     0.805135781356946*I,
        -0.0136222082924246 +     0.458040841381709*I,
         -0.059143991192552 +     0.808570036690227*I,
       -0.00256684120826156 -     0.470316855748848*I,
          0.736631549065503 +     0.424083338792149*I,
          0.113125772574669 +    0.0579047421322912*I,
         -0.123508004771838 -     0.908270882204866*I,
          0.161467503729708 +    0.0540828673048715*I,
          0.939354119069586 +     0.457659292402481*I };
    /* Matlab code to generate result_exact:
clear lam lambda result_exact
clear
eps_t = 0.13;
    kappa = -1; D=8;
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
end
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
    if (misc_rel_err(18, result, result_exact) > 1000*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (manakov_scatter_matrix_test_defocusing_CF4_2() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
