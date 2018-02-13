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

static INT nse_fscatter_test_focusing_2split2_modal()
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
      -0.000621226975062701 -   0.00135796924692729*I,
       -0.00388866907570591 -  0.000115823138793798*I,
       -0.00629791886669818 +   0.00370003110637529*I,
       -0.00643495731890653 +    0.0075008058387984*I,
       -0.00743509432141662 +    0.0084192796324154*I,
        -0.0134132704091659 +   0.00586897194982142*I,
        -0.0234955391103155 +   0.00205737763730251*I,
          0.984835331462689 +                     0*I,
        -0.0074512658492344 +    0.0432393003450808*I,
         0.0389531486115869 +    0.0543514135900098*I,
         0.0492546645712511 +    0.0608472810955837*I,
         0.0146703321016289 +    0.0622380374023295*I,
        -0.0332825430879091 +    0.0582957550972926*I,
        -0.0509603967571021 +    0.0492653189605623*I,
        -0.0220358075771321 +     0.035793541521466*I,
         0.0276696576254646 +    0.0189175181442808*I,
                          0 +                     0*I,
                          0 +                     0*I,
        -0.0276696576254646 +    0.0189175181442808*I,
         0.0220358075771321 +     0.035793541521466*I,
         0.0509603967571021 +    0.0492653189605623*I,
         0.0332825430879091 +    0.0582957550972926*I,
        -0.0146703321016289 +    0.0622380374023295*I,
        -0.0492546645712511 +    0.0608472810955837*I,
        -0.0389531486115869 +    0.0543514135900098*I,
         0.0074512658492344 +    0.0432393003450808*I,
          0.984835331462689 +                     0*I,
        -0.0234955391103155 -   0.00205737763730251*I,
        -0.0134132704091659 -   0.00586897194982142*I,
       -0.00743509432141662 -    0.0084192796324154*I,
       -0.00643495731890653 -    0.0075008058387984*I,
       -0.00629791886669818 -   0.00370003110637529*I,
       -0.00388866907570591 +  0.000115823138793798*I,
      -0.000621226975062701 +   0.00135796924692729*I,
                          0 +                     0*I };
    /* Matlab code to generate result_exact:
    p11 = 1; p12 = 0; p21 = 0; p22 = 1; eps_t = 0.13;
    kappa = +1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
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
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, NULL, nse_discretization_2SPLIT2_MODAL);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*DBL_EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, &W, nse_discretization_2SPLIT2_MODAL);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W == 0)
        return E_TEST_FAILED;
    scl = POW(2.0, W);
    for (i=0; i<4*(deg+1); i++)
        result[i] *= scl;
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*DBL_EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (nse_fscatter_test_focusing_2split2_modal() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
