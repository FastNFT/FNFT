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

static INT poly_fmult2x2_test_n_is_no_power_of_2(INT normalize_flag)
{
    UINT deg = 1, n = 5;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:
    i = 0:9;
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

    u11 = conv(s11,t11) + conv(s12,t21);
    u12 = conv(s11,t12) + conv(s12,t22);
    u21 = conv(s21,t11) + conv(s22,t21);
    u22 = conv(s21,t12) + conv(s22,t22);

    ind = 9:10;
    r11 = conv(u11,p11(ind)) + conv(u12,p21(ind));
    r12 = conv(u11,p12(ind)) + conv(u12,p22(ind));
    r21 = conv(u21,p11(ind)) + conv(u22,p21(ind));
    r22 = conv(u21,p12(ind)) + conv(u22,p22(ind));

    format long g; result_exact = [r11 r12 r21 r22].'
    */
    const UINT memsize = poly_fmult2x2_numel(deg, n);
    COMPLEX p[memsize];
    COMPLEX result[memsize];
    COMPLEX result_exact[24] = {
          -163.292988790261 -      31.2370829948429*I, 
          -984.955455390257 +      114.750671387418*I,
          -1258.50583714315 +      1286.47060781892*I,
          -171.431868738405 +      1460.18133348016*I,
           824.665563851057 -      164.619215022034*I,
          -677.169924780033 -      633.214936468468*I,
          -157.736448644796 -      105.591543065528*I,
          -1124.70654776703 +      9.72919735381461*I,
          -1445.87975819902 +      1317.90341910231*I,
          -155.534713048327 +      1635.14192990041*I,
           929.979021480405 -        257.7534144306*I,
          -719.402982559702 -      649.205484178286*I,
          -153.832150334564 -      63.0557298389851*I,
          -909.619744624336 -      30.7107489877624*I,
          -1071.61910789952 +      1239.31431276425*I,
          -102.976412016802 +      1415.13842812132*I,
           992.760691268736 -      241.500041506092*I,
          -764.995532068071 -      490.694096536849*I,
           -133.61442022569 -      134.824136955459*I,
          -1045.09531857318 -      123.739126689896*I,
          -1249.33331746925 +      1267.46980302693*I,
          -73.7197104407829 +      1607.82704486421*I,
           1089.20607658251 -      351.290776172497*I,
          -808.546624609249 -      498.637460827699*I };
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/4 - 1;

    for (i=0; i<n*(deg+1); i++) {
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
 
    ret_code = poly_fmult2x2_test_n_is_no_power_of_2(0); // 2x2 w/o normaliz.
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult2x2_test_n_is_no_power_of_2(1); // 2x2 with normaliz.
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
