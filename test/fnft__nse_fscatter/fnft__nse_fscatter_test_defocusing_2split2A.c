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
* Sander Wahls (TU Delft) 2017-2018.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__nse_fscatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static INT nse_fscatter_test_defocusing_2split2A()
{
    UINT i, D = 8, deg;
    INT W;
    REAL scl;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8];
    COMPLEX result[4*2*8];
    COMPLEX result_exact[36] = {
                          0 +                     0*I,
       0.000639850511863557 +   0.00139867931146051*I,
        0.00401458751274878 +  0.000186256997598998*I,
        0.00659872422074629 -   0.00370818909843262*I,
        0.00687092135744545 -   0.00763968758325248*I,
        0.00799923534398637 -   0.00861693168786061*I,
         0.0141285634061023 -   0.00601148920849778*I,
         0.0241474208616003 -   0.00210777486264419*I,
           1.01542228581741 +                     0*I,
       -0.00767760770342658 +    0.0445527501148382*I,
          0.039391619284806 +    0.0578159182259416*I,
         0.0505007832197083 +    0.0656803906040225*I,
         0.0148098903060824 +     0.067415207521504*I,
        -0.0346462389904789 +    0.0628353020665661*I,
        -0.0519037555064953 +    0.0525189049051206*I,
        -0.0212160144049619 +    0.0376066095495185*I,
         0.0285180097182471 +    0.0194975295172145*I,
                          0 +                     0*I,
                          0 +                     0*I,
         0.0285180097182471 -    0.0194975295172145*I,
        -0.0212160144049619 -    0.0376066095495185*I,
        -0.0519037555064953 -    0.0525189049051206*I,
        -0.0346462389904789 -    0.0628353020665661*I,
         0.0148098903060824 -     0.067415207521504*I,
         0.0505007832197083 -    0.0656803906040225*I,
          0.039391619284806 -    0.0578159182259416*I,
       -0.00767760770342658 -    0.0445527501148382*I,
           1.01542228581741 +                     0*I,
         0.0241474208616003 +   0.00210777486264419*I,
         0.0141285634061023 +   0.00601148920849778*I,
        0.00799923534398637 +   0.00861693168786061*I,
        0.00687092135744545 +   0.00763968758325248*I,
        0.00659872422074629 +   0.00370818909843262*I,
        0.00401458751274878 -  0.000186256997598998*I,
       0.000639850511863557 -   0.00139867931146051*I,
                          0 +                     0*I };
    /* Matlab code to generate result_exact:
    p11 = 1; p12 = 0; p21 = 0; p22 = 1; eps_t = 0.13;
    kappa = -1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    Delta = sqrt(-kappa)*eps_t*abs(q);
    rho = eps_t*sinh(Delta)./Delta;
    rho(Delta == 0) = eps_t;
    for i=1:D
        t11 = conv([0 cosh(Delta(i))],p11)+conv([rho(i)*q(i) 0], p21);
        t12 = conv([0 cosh(Delta(i))],p12)+conv([rho(i)*q(i) 0], p22);
        t21 = conv([0 -kappa*rho(i)*q(i)'],p11)+conv([cosh(Delta(i)) 0], p21);
        t22 = conv([0 -kappa*rho(i)*q(i)'],p12)+conv([cosh(Delta(i)) 0], p22);
        p11 = t11; p12 = t12; p21 = t21; p22 = t22;
    end
    format long g; result_exact = [p11 p12 p21 p22].'
    */
    for (i=0; i<D; i++)
        q[i] = 0.4*COS(i+1) + 0.5*I*SIN(0.3*(i+1));

    // without normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, NULL, nse_discretization_2SPLIT2A);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*DBL_EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, &W, nse_discretization_2SPLIT2A);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W != 0) // it should really be zero in this case
        return E_TEST_FAILED;
    scl = POW(2.0, W);
    for (i=0; i<4*(deg+1); i++)
        result[i] *= scl;
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (nse_fscatter_test_defocusing_2split2A() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
