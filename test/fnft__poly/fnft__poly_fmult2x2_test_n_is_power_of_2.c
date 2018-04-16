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

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static INT poly_fmult2x2_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 1, n = 4;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:
    i = 0:7;
    p11 = sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
    p12 = sqrt(i+1).*(cos(i+0.1) + 1j*sin(-2*i+0.1));
    p21 = sqrt(i+1).*(cos(i+0.2) + 1j*sin(-2*i+0.2));
    p22 = sqrt(i+1).*(cos(i+0.3) + 1j*sin(-2*i+0.3));

    ind1 = 1:2; ind2 = 3:4;
    s11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2));
    s12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2));
    s21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2));
    s22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2));

    ind1 = 5:6; ind2 = 7:8;
    t11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2));
    t12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2));
    t21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2));
    t22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2));

    r11 = conv(s11,t11) + conv(s12,t21);
    r12 = conv(s11,t12) + conv(s12,t22);
    r21 = conv(s21,t11) + conv(s22,t21);
    r22 = conv(s21,t12) + conv(s22,t22);

    format long g; result_exact = [r11 r12 r21 r22].'
    */
    const UINT memsize = poly_fmult2x2_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    COMPLEX result_exact[20] = {
        60.6824426714241 + I*64.8661118935554,
        192.332936227757 - I*1.10233198818688,
        162.110387625956 - I*111.127467380607,
       -89.7931323551994 - I*80.292426677456,
        10.9405573209266 + I*121.130913146948,
        60.232171606396 + I*70.7074200214283,
        188.004055601206 + I*7.96812436286603,
        156.207023899884 - I*102.153018093302,
       -100.627468409158 - I*80.5748163736639,
        14.8569933196036 + I*113.838094358267,
        46.5859268828598 + I*75.6288485779254,
        171.341929283725 - I*10.2410011372218,
        155.475099221666 - I*110.340545813463, 
       -117.458524441818 - I*87.5204294848652,
        31.9873641358437 + I*114.857725853133, 
          44.9841424843245 + I*81.2642643937605, 
          166.826592700681 - I*5.16720885749995, 
          155.130568557566 - I*102.464890587643, 
          -127.981613283912 - I*83.7735775388555, 
          34.4728506991058 + I*107.132856593172};
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/4 - 1;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
        p[i+n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.1) + I*SIN(-2.0*i+0.1));
        p[i+2*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.2) + I*SIN(-2.0*i+0.2));
        p[i+3*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.3) + I*SIN(-2.0*i+0.3));
    }
   
    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult2x2(&deg, n, p, result, W_ptr);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (deg != deg_exact)
        return E_TEST_FAILED;
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<4*(deg+1); i++)
            result[i] *= scl;
    }
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult2x2_test_n_is_power_of_2(0); // 2x2 w/o normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult2x2_test_n_is_power_of_2(1); // 2x2 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
