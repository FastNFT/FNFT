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

static INT poly_fmult_test_D_is_power_of_2(INT normalize_flag)
{
    UINT deg = 2, n = 8, i;
    COMPLEX p[24];
    /* Matlab code to generate the "exact" result:
        i=0:23; p=sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
        r=p(1:3); for (k=4:3:24), r=conv(r,p(k:(k+2))); end
        format long g; r.'
    */
    COMPLEX result_exact[18] = { \
           10690.1093069785 - I*7285.45565246217,
           22605.6811426431 - I*19510.065631389, 
           57389.5004107115 - I*20763.2575587549,
           53847.0938223071 - I*18963.5560762862,
           88648.7268642591 - I*2309.98190237675,
           66515.8698354361 - I*3878.11603829108,
           45160.5175520559 + I*12083.5337878339,
           64021.1527510823 + I*34487.976136233, 
          -26614.7840829275 + I*18912.0673174305,
           624.669240812309 + I*18294.2454730283,
           7755.19752834515 + I*18447.5313133858,
          -33812.2918248361 + I*1257.46809193287,
          -711.103319555697 - I*2023.66741077336,
           9781.54461527385 + I*407.75394285625, 
           39.3498537600159 + I*420.245708764338,
          -508.134987301475 + I*156.212629759735,
           4.42978955469815 - I*0.0813363246183392 };
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
    }
    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult(&deg, n, p, W_ptr);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<deg+1; i++)
            p[i] *= scl;
    }
    if (misc_rel_err(deg+1, p, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;
    return SUCCESS;
}

static INT poly_fmult_test_D_is_no_power_of_2(INT normalize_flag)
{
    UINT deg = 2, n = 11, i;
    COMPLEX p[(deg+1)*16]; // 16 = next power of two of n=11
    /* Matlab code to generate the "exact" result:
        deg=2; n=11; N=(deg+1)*n; i=0:(N-1);
        p=sqrt(i+1).*(cos(i) + 1j*sin(-2*i)); r=p(1:(deg+1));
        for (k=(deg+2):(deg+1):N), r=conv(r,p(k:(k+deg))); end
        format long g; r.'
    */
    COMPLEX result_exact[23] = {
          -319308.705234574 -      166784.182995777*I, 
          -2743910.31810172 -      794705.736952064*I,
           -9075447.9140096 -      882018.812864758*I,
          -18191941.8165097 -      122049.064791405*I,
          -30115074.8590672 +      760936.799294374*I,
          -44814146.2207504 +      1917904.17878258*I,
          -59509427.3900375 +       1151280.7345219*I,
          -67870895.8105611 +       6758914.6822708*I,
          -65004270.8979727 +      2378651.90695769*I,
           -62350519.934101 -      4351181.14834246*I,
           -49802387.816833 +      6324801.48028804*I,
          -25632453.5014797 -      15474491.9397379*I,
          -19770246.1054521 -      11669704.6069006*I,
          -5178072.79922669 +      454130.609010228*I,
           7883730.23930349 -      23609117.6699819*I,
          -1012854.76439042 -      6787648.75896478*I,
           2010055.29299395 +      1306386.13720826*I,
           4346918.37729613 -      8535073.39464626*I,
          -2082296.12229051 -      649844.924066411*I,
          -1785325.93258564 +      2098730.17691877*I,
            50207.974399346 -      84587.5774434785*I,
           56394.1550202597 -      147664.217837454*I,
          -801.466849402119 +      1048.08952721767*I };        
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
    }
    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult(&deg, n, p, W_ptr);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<deg+1; i++)
            p[i] *= scl;
    }
    if (misc_rel_err(deg+1, p, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;
    return SUCCESS;
}

static INT poly_fmult2x2_test(INT normalize_flag)
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
    COMPLEX p[32], result[2*32];
    COMPLEX result_exact[20] = { \
        60.6824426714241 + I*64.8661118935554, \
        192.332936227757 - I*1.10233198818688, \
        162.110387625956 - I*111.127467380607, \
       -89.7931323551994 - I*80.292426677456, \
        10.9405573209266 + I*121.130913146948, \
        60.232171606396 + I*70.7074200214283, \
        188.004055601206 + I*7.96812436286603, \
        156.207023899884 - I*102.153018093302, \
       -100.627468409158 - I*80.5748163736639, \
        14.8569933196036 + I*113.838094358267, \
        46.5859268828598 + I*75.6288485779254, \
        171.341929283725 - I*10.2410011372218, \
        155.475099221666 - I*110.340545813463, \
       -117.458524441818 - I*87.5204294848652, \
        31.9873641358437 + I*114.857725853133, \
          44.9841424843245 + I*81.2642643937605, \
          166.826592700681 - I*5.16720885749995, \
          155.130568557566 - I*102.464890587643, \
          -127.981613283912 - I*83.7735775388555, \
          34.4728506991058 + I*107.132856593172};

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
        p[i+8] = SQRT(i+1.0)*(COS(i+0.1) + I*SIN(-2.0*i+0.1));
        p[i+16] = SQRT(i+1.0)*(COS(i+0.2) + I*SIN(-2.0*i+0.2));
        p[i+24] = SQRT(i+1.0)*(COS(i+0.3) + I*SIN(-2.0*i+0.3));
    }
   
    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult2x2(&deg, n, p, result, W_ptr);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<20; i++)
            result[i] *= scl;
    }
    if (misc_rel_err(20, result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult_test_D_is_no_power_of_2(0); // 1x1 w/o normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_D_is_no_power_of_2(1); // 1x1 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_D_is_power_of_2(0); // 1x1 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_D_is_power_of_2(1); // 1x1 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_D_is_no_power_of_2(1); // 1x1 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult2x2_test(0); // 2x2 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult2x2_test(1); // 2x2 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
