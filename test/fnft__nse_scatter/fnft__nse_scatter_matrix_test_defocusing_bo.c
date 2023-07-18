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
* Sander Wahls (TU Delft) 2017-2018, 2020.
* Shrinivas Chimmalgi (TU Delft) 2017-2018,2020.
* Sander Wahls (KIT) 2023.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT nse_scatter_matrix_test_defocusing_bo(INT normalization_flag)
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8], r[8];
    COMPLEX result[2*8];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[16] = { \

         -0.475166844884886 -      0.94058850193694*I, \
        -0.0737839235972917 +     0.324108304540972*I, \
        -0.0737839235972917 -     0.324108304540972*I, \
         -0.475166844884886 +      0.94058850193694*I, \
         -0.953568191102039 +     0.514608486676543*I, \
        -0.0296784708985577 -     0.102190009894233*I, \
        -0.0296784708985577 +     0.102190009894233*I, \
         -0.953568191102039 -     0.514608486676543*I, \
          0.922265187166008 -      1.53250405713011*I, \
       -0.00545204324801093 +     0.390264205099379*I, \
        -0.0587177918688801 -     0.442673550385969*I, \
          0.348002709127808 +       0.5560370203054*I, \
          -1.56273289561246 -     0.921663361441779*I, \
        -0.0131572520813856 -    0.0545206565775057*I, \
        -0.0972189335180089 +     0.061610395802929*I, \
         -0.560360759744643 +     0.327099280806848*I };
    /* Matlab code to generate result_exact:
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
    INT W[8];
    INT * const W_ptr = (normalization_flag) ? W : NULL;
    for (i=0; i<D; i++) {
        q[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));
        r[i] = CONJ(q[i]);
    }

    ret_code = nse_scatter_matrix(D, q, r, eps_t, -1, 2, lam, result, W_ptr, nse_discretization_BO, 1);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (normalization_flag) {
        REAL scl = POW(2, W[0]);
        for (i=0; i<8; i++)
            result[i] *= scl;
        scl = POW(2, W[1]);
        for (i=8; i<16; i++)
            result[i] *= scl;
    }
    if (misc_rel_err(16, result, result_exact) > 10*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    INT ret_code = nse_scatter_matrix_test_defocusing_bo(0);
    CHECK_RETCODE(ret_code, failure);

    ret_code = nse_scatter_matrix_test_defocusing_bo(1);
    CHECK_RETCODE(ret_code, failure);

    return EXIT_SUCCESS;

failure:
    return EXIT_SUCCESS;
}
