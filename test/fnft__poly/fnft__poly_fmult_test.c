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

/*INT test_poly_fmult2()
{
    const INT deg = 4;
    kiss_fft_cpx p1[5] = {1, 0, 0, 2, 3, 0, 0, 4, 5, 0};
    kiss_fft_cpx p2[5] = {0, 5, 6, 0, 0, 7, 8, 0, 0, 9};
    kiss_fft_cpx result[9];
    kiss_fft_cpx result_exact[9] = {0, 5, -4, 0, 0, 34, -8, 0, 0, 95, 8, 0, \
        0, 94, 4, 0, 0, 45};
    void *mem;
    UINT lenmem;

    lenmem = poly_fmult2_lenmen(deg);
    mem = malloc(lenmem);
    if (mem == NULL)
        return E_NOMEM;

    if (poly_fmult2(deg, p1, p2, result, mem, lenmem, 0) == FNFT_ERROR)
        return E_SUBROUTINE;
    if (rel_err(2*deg+1, result, result_exact) > 10*DBL_EPSILON)
        return FNFT_E_INCORRECT_RESULT;
    return SUCCESS;
}*/

static INT poly_fmult_test(INT normalize_flag)
{
    UINT deg = 2, n = 8, i;
    COMPLEX p[24];
    /* Matlab code to generate the "exact" result:
        i=0:23; p=sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
        r=p(1:3); for (k=4:3:24), r=conv(r,p(k:(k+2))); end
        format long g; r.'
    */
    COMPLEX result_exact[18] = { \
           10690.1093069785 - I*7285.45565246217, \
           22605.6811426431 - I*19510.065631389, \
           57389.5004107115 - I*20763.2575587549, \
           53847.0938223071 - I*18963.5560762862, \
           88648.7268642591 - I*2309.98190237675, \
           66515.8698354361 - I*3878.11603829108, \
           45160.5175520559 + I*12083.5337878339, \
           64021.1527510823 + I*34487.976136233, \
          -26614.7840829275 + I*18912.0673174305, \
           624.669240812309 + I*18294.2454730283, \
           7755.19752834515 + I*18447.5313133858, \
          -33812.2918248361 + I*1257.46809193287, \
          -711.103319555697 - I*2023.66741077336, \
           9781.54461527385 + I*407.75394285625, \
           39.3498537600159 + I*420.245708764338, \
          -508.134987301475 + I*156.212629759735, \
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
    COMPLEX p[32], result[24];
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
    /*printf("poly_fmult2: %s\n", test_poly_fmult2()==FNFT_ERROR \
        ? "test failed" : "test passed");*/

    ret_code = poly_fmult_test(0); // 1x1 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test(1); // 1x1 without normalization
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
