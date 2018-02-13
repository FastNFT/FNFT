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

static INT nse_fscatter_test_focusing_2split2A()
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
        -0.0006218488666558 -   0.00135932866899402*I,
       -0.00389486465722303 -  0.000116415544549393*I,
       -0.00631276117670343 +   0.00370814567471042*I,
       -0.00644969591990114 +   0.00752481420549383*I,
       -0.00745312456702635 +   0.00845167335191335*I,
        -0.0134524889928082 +   0.00589367207333885*I,
        -0.0235703477995799 +    0.0020704635641732*I,
          0.984788462154351 +                     0*I,
       -0.00745584493130653 +    0.0432658725153708*I,
         0.0390132523297099 +    0.0544316882818322*I,
         0.0493595931916023 +    0.0609730121105069*I,
         0.0146919904705362 +    0.0623212173953222*I,
        -0.0333333805616417 +    0.0583826525051374*I,
        -0.0510476257002306 +    0.0493474287305426*I,
        -0.0220503548310495 +    0.0358125415341774*I,
         0.0276790289397556 +    0.0189239252350563*I,
                          0 +                     0*I,
                          0 +                     0*I,
        -0.0276790289397556 +    0.0189239252350563*I,
         0.0220503548310495 +    0.0358125415341774*I,
         0.0510476257002306 +    0.0493474287305426*I,
         0.0333333805616417 +    0.0583826525051374*I,
        -0.0146919904705362 +    0.0623212173953222*I,
        -0.0493595931916023 +    0.0609730121105069*I,
        -0.0390132523297099 +    0.0544316882818322*I,
        0.00745584493130653 +    0.0432658725153708*I,
          0.984788462154351 +                     0*I,
        -0.0235703477995799 -    0.0020704635641732*I,
        -0.0134524889928082 -   0.00589367207333885*I,
       -0.00745312456702635 -   0.00845167335191335*I,
       -0.00644969591990114 -   0.00752481420549383*I,
       -0.00631276117670343 -   0.00370814567471042*I,
       -0.00389486465722303 +  0.000116415544549393*I,
        -0.0006218488666558 +   0.00135932866899402*I,
                          0 +                     0*I };
    /* Matlab code to generate result_exact:
    p11 = 1; p12 = 0; p21 = 0; p22 = 1; eps_t = 0.13;
    kappa = +1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
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
        q[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));

    // without normalization
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, NULL, nse_discretization_2SPLIT2A);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, &W, nse_discretization_2SPLIT2A);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W == 0)
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
    if (nse_fscatter_test_focusing_2split2A() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
