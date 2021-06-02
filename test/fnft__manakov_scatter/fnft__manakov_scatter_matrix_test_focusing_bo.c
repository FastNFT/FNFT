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
INT manakov_scatter_matrix_test_focusing_bo()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    const INT kappa = 1;
    COMPLEX q1[8], q2[8];
    COMPLEX result[18];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[18] = { -0.590289500527102 -     0.632629148962372*I,
        -0.0511893048465224 +     0.279426275135839*I,
         -0.209603512694856 +     0.355957363668216*I,
         0.0687240524043816 +     0.286905729978922*I,
         -0.497489772558467 +      0.80970142687034*I,
        -0.0226543315497955 -    0.0966251184334277*I,
           0.19819895871805 +     0.353578767660151*I,
        -0.0599351238096392 -     0.112231538752811*I,
         -0.580126285738488 +     0.694956533217276*I,
          0.514476062248241 -      1.17775938388385*I,
        0.00290970524770017 +     0.331371924186705*I,
        -0.0566170146699455 +     0.553800705523123*I,
         0.0540346314461356 +     0.388596516118113*I,
          0.258219292512569 +     0.471820101061512*I,
        -0.0734657706676645 -    0.0606847546267588*I,
          0.231635690137041 +     0.671583582363646*I,
         -0.110755098081325 -    0.0603997043245748*I,
          0.116003702766256 +       0.4211835584256*I };
    /* Matlab code to generate result_exact:
    clear
    eps_t = 0.13;
    kappa = +1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
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

    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(18, result, result_exact) > 10*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (manakov_scatter_matrix_test_focusing_bo() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
