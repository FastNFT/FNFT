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

static INT nse_fscatter_test_defocusing_2split2_modal()
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
       0.000640551624374837 +   0.00140021190625625*I,
        0.00402138731457545 +  0.000187175671792273*I,
        0.00661517316353962 -   0.00371641703326692*I,
        0.00688808601726169 -   0.00766473549226468*I,
        0.00802055965788163 -   0.00865094187719025*I,
         0.0141722122264573 -   0.00603748291258251*I,
         0.0242264202084207 -   0.00212137695307512*I,
           1.01547082891319 +                     0*I,
       -0.00768305407680358 +    0.0445843551305502*I,
           0.03945396460556 +    0.0579111756213562*I,
         0.0506122555584725 +    0.0658294591284749*I,
         0.0148317044630103 +    0.0675237362559048*I,
        -0.0347048112757902 +    0.0629437931581923*I,
        -0.0519964361407694 +     0.052615049271045*I,
        -0.0212275567248036 +    0.0376321554770775*I,
         0.0285303839809887 +    0.0195059896992343*I,
                          0 +                     0*I,
                          0 +                     0*I,
         0.0285303839809887 -    0.0195059896992343*I,
        -0.0212275567248036 -    0.0376321554770775*I,
        -0.0519964361407694 -     0.052615049271045*I,
        -0.0347048112757902 -    0.0629437931581923*I,
         0.0148317044630103 -    0.0675237362559048*I,
         0.0506122555584725 -    0.0658294591284749*I,
           0.03945396460556 -    0.0579111756213562*I,
       -0.00768305407680358 -    0.0445843551305502*I,
           1.01547082891319 +                     0*I,
         0.0242264202084207 +   0.00212137695307512*I,
         0.0141722122264573 +   0.00603748291258251*I,
        0.00802055965788163 +   0.00865094187719025*I,
        0.00688808601726169 +   0.00766473549226468*I,
        0.00661517316353962 +   0.00371641703326692*I,
        0.00402138731457545 -  0.000187175671792273*I,
       0.000640551624374837 -   0.00140021190625625*I,
                          0 +                     0*I };
    /* Matlab code to generate result_exact:
    p11 = 1; p12 = 0; p21 = 0; p22 = 1; eps_t = 0.13;
    kappa = -1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    Q = eps_t * q; for i=1:D
        scl = 1/sqrt(1 + kappa*abs(Q(i))^2);
        t11 = conv([0 scl],p11)+conv([scl*Q(i) 0], p21);
        t12 = conv([0 scl],p12)+conv([scl*Q(i) 0], p22);
        t21 = conv([0 -kappa*scl*Q(i)'],p11)+conv([scl 0], p21);
        t22 = conv([0 -kappa*scl*Q(i)'],p12)+conv([scl 0], p22);
        p11 = t11; p12 = t12; p21 = t21; p22 = t22;
    end
    format long g; result_exact = [p11 p12 p21 p22].'
    */

    for (i=0; i<D; i++)
        q[i] = 0.4*COS(i+1) + 0.5*I*SIN(0.3*(i+1));

    // without normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, NULL, nse_discretization_2SPLIT2_MODAL);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, &W, nse_discretization_2SPLIT2_MODAL);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W != 0) // it really should be zero in this case!
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
    if (nse_fscatter_test_defocusing_2split2_modal() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
