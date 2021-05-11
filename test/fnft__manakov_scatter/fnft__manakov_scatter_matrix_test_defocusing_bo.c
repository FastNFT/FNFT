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
INT manakov_scatter_matrix_test_defocusing_bo()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    const INT kappa = -1;
    COMPLEX q1[8], q2[8];
    COMPLEX result[18];
    COMPLEX lam[2] = {2, 1+0.5*I};
    // TODO: update matlab code and result
    COMPLEX result_exact[18] = { -0.358090090099409 -      1.13917188728172*I,
        -0.0936860689277801 +       0.3557030727293*I,
         -0.255375942444119 +     0.474786766834309*I,
        -0.0748769150154918 -     0.347254551668194*I,
          -0.47398368264259 +     0.943119274834602*I,
         0.0299267209047495 +     0.105651168741426*I,
         -0.267845497420234 -     0.477501429166025*I,
         0.0729420978608351 +     0.125883707618941*I,
         -0.371424472524121 +      1.06800066903998*I,
           1.23905565295345 -      1.74828015604221*I,
         -0.014909213651051 +      0.43237563683829*I,
        -0.0666271402229915 +     0.713230897578002*I,
        -0.0594717745495415 -     0.476084941782815*I,
          0.350489581988108 +     0.557729604600123*I,
         0.0848592596301212 +    0.0656854715680317*I,
         -0.305372161239983 -     0.870389481046401*I,
          0.129916297678531 +    0.0667430689661592*I,
          0.516789279241878 +     0.610111903299273*I };
    /* Matlab code to generate result_exact:
    clear
    eps_t = 0.13;
    kappa = -1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    result_exact = []; lam = [2 , 1+0.5i];
    for l = lam
        S=eye(4);
        for n=D:-1:1
            ks=(-kappa*abs(q(n))^2-l^2);
            k=sqrt(ks);
            U=[cosh(k*eps_t)-1i*l*sinh(k*eps_t)/k,q(n)*sinh(k*eps_t)/k;...
                -kappa*conj(q(n))*sinh(k*eps_t)/k,cosh(k*eps_t)+1i*l*sinh(k*eps_t)/k];
            Udash=[1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t+1i+1i*l^2/ks)*sinh(k*eps_t)/k,...
                -q(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks;...
                kappa*conj(q(n))*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks,...
                -1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t-1i-1i*l^2/ks)*sinh(k*eps_t)/k];
            T=[U,zeros(2);Udash,U];
            S=S*T;
        end
        result_exact = [result_exact,[S(1,1) S(1,2) S(2,1) S(2,2) S(3,1) S(3,2) S(4,1) S(4,2)]];
    end
    format long g; result_exact.'
    */
    for (i=0; i<D; i++) {
        q1[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));
        q2[i] = 0.21*cos(i+1) + 1.05*I*sin(0.2*(i+1));
    }

    ret_code = manakov_scatter_matrix(D, q1, q2, eps_t, 2, lam, kappa, result, manakov_discretization_BO);
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
    if (manakov_scatter_matrix_test_defocusing_bo() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
