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
 * Sander Wahls (TU Delft) 2017-2018 (fmult2x2 test files).
 * Lianne de Vries (TU Delft) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include <stdio.h>

static INT poly_fmult3x3_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 1, n = 4;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:
    % Define the polynomial coefficients:
% p11 = [p11_1(1) p11_1(0) p11_2(1) ... p11_4(0)]
% where the subscript indicates the number of the matrix and (.) the
% coefficient. So p11_4(1) is the coefficient in front of z^1 of the upper
% left element of matrix 4. In C, all these pij are collected in one
% long vector: p = [p11 p12 ... p33]. Pointers *pij indicate the memory
% address where each pij starts
i = 0:7;
    p11 = sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
    p12 = sqrt(i+1).*(cos(i+0.1) + 1j*sin(-2*i+0.1));
    p13 = sqrt(i+1).*(cos(i+0.2) + 1j*sin(-2*i+0.2));
    p21 = sqrt(i+1).*(cos(i+0.3) + 1j*sin(-2*i+0.3));
    p22 = sqrt(i+1).*(cos(i+0.4) + 1j*sin(-2*i+0.4));
    p23 = sqrt(i+1).*(cos(i+0.5) + 1j*sin(-2*i+0.5));
    p31 = sqrt(i+1).*(cos(i+0.6) + 1j*sin(-2*i+0.6));
    p32 = sqrt(i+1).*(cos(i+0.7) + 1j*sin(-2*i+0.7));
    p33 = sqrt(i+1).*(cos(i+0.8) + 1j*sin(-2*i+0.8));
    
% First multiply matrices 1 and 2:
    ind1 = 1:2; ind2 = 3:4;
    s11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    s12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    s13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    s21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    s22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    s23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    s31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    s32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    s33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    ind1 = 5:6; ind2 = 7:8;
    t11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    t12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    t13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    t21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    t22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    t23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    t31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    t32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    t33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    r11 = conv(s11,t11) + conv(s12,t21) + conv(s13,t31);
    r12 = conv(s11,t12) + conv(s12,t22) + conv(s13,t32);
    r13 = conv(s11,t13) + conv(s12,t23) + conv(s13,t33);
    r21 = conv(s21,t11) + conv(s22,t21) + conv(s23,t31);
    r22 = conv(s21,t12) + conv(s22,t22) + conv(s23,t32);
    r23 = conv(s21,t13) + conv(s22,t23) + conv(s23,t33);
    r31 = conv(s31,t11) + conv(s32,t21) + conv(s33,t31);
    r32 = conv(s31,t12) + conv(s32,t22) + conv(s33,t32);
    r33 = conv(s31,t13) + conv(s32,t23) + conv(s33,t33);

    format long g; result_exact = [r11 r12 r13 r21 r22 r23 r31 r32 r33]'
    */

    const UINT memsize = poly_fmult3x3_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    COMPLEX result_exact[45] = {
        27.1816382503119 +      I*257.127029644935,
           495.294918273546 +       I*370.79414146898,
           449.985211066113 -      I*82.3177853061955,
          -137.718854031403 +      I*176.402462123245,
           295.374188241712 +      I*230.810890384916,
           17.4941046840742 +      I*265.055038507231,
           465.358111926116 +      I*385.041795839986,
           401.320501169614 -      I*70.5740374634787,
          -182.354169142207 +       I*170.05221382589,
           280.765918587179 +      I*197.479740256079,
           7.63177580661544 +       I*270.33470504031,
           430.771601151261 +      I*395.442239864855,
           348.645929484348 -      I*58.1251371671016,
          -225.167461673215 +      I*162.002860019777,
           263.352328682996 +      I*162.175437840726,
          -50.0186220841524 +      I*253.675557113414,
           406.358796093538 +      I*277.247044099463,
            432.48534622583 -      I*176.148770488459,
          -220.301480277088 +      I*160.839194819966,
           335.308613960652 +       I*135.92414782824,
           -61.616363207007 +      I*258.386611344143,
             383.7260605612 +      I*277.811964436761,
           410.918340391173 -       I*165.97365998068,
          -253.323920091272 +      I*176.555618777935,
           312.916595162616 +      I*108.861808336948,
           -72.598453996359 +       I*260.51595196558,
           357.259261074707 +      I*275.601079457834,
           385.245574330871 -      I*154.140195525971,
          -283.815231033659 +      I*190.507957354603,
            287.39801718219 +      I*80.7117576416826,
          -122.750867875861 +      I*227.564002574487,
            281.12385289707 +      I*158.934293991886,
             376.3528534572 -      I*254.244910618278,
          -283.205231402771 +      I*130.908641262946,
           345.290919827773 +      I*28.8957059623247,
          -135.222824881815 +      I*228.637277729921,
           267.816903038957 +      I*145.766017644221,
           383.810068083634 -      I*246.547349763061,
          -301.664999720857 +      I*167.287836131531,
           317.115364236406 +      I*10.5195752968896,
          -146.343680120272 +      I*227.426084783715,
           251.834015214176 +      I*131.141295432093,
           387.432379380516 -      I*236.386369286729,
          -317.110631078923 +      I*201.995546239497,
           285.771296749992 -      I*7.96166348760254};
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/9 - 1;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
        p[i+n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.1) + I*SIN(-2.0*i+0.1));
        p[i+2*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.2) + I*SIN(-2.0*i+0.2));
        p[i+3*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.3) + I*SIN(-2.0*i+0.3));
        p[i+4*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.4) + I*SIN(-2.0*i+0.4));
        p[i+5*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.5) + I*SIN(-2.0*i+0.5));
        p[i+6*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.6) + I*SIN(-2.0*i+0.6));
        p[i+7*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.7) + I*SIN(-2.0*i+0.7));
        p[i+8*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.8) + I*SIN(-2.0*i+0.8));
    }

    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult3x3(&deg, n, p, result, W_ptr);
    if (ret_code != SUCCESS){
        return E_SUBROUTINE(ret_code);
    }
    if (deg != deg_exact){
        return E_TEST_FAILED;
    }
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<9*(deg+1); i++)
            result[i] *= scl;
    }
    #define DEBUG
    #ifdef DEBUG
        printf("error = %2.1e < %2.1e\n",misc_rel_err(45, result, result_exact),100*EPSILON);
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }
    #endif
    if (!(misc_rel_err(9*(deg+1), result, result_exact) <= 100*EPSILON)){
        REAL err;
        err = misc_rel_err(9*(deg+1), result, result_exact);
        return E_TEST_FAILED;
    }

    return SUCCESS;
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult3x3_test_n_is_power_of_2(0); // 3x3 w/o normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    
    ret_code = poly_fmult3x3_test_n_is_power_of_2(1); // 3x3 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
