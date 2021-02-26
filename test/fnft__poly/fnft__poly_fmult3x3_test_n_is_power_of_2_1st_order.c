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

 /* This code checks for 2 first-order polynomial matrices of dim 3x3
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include <stdio.h>

static INT poly_fmult3x3_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 1, n = 2;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:
    % Define the polynomial coefficients:
% p11 = [p11_1(1) p11_1(0) p11_2(1) ... p11_4(0)]
% where the subscript indicates the number of the matrix and (.) the
% coefficient. So p11_4(1) is the coefficient in front of z^1 of the upper
% left element of matrix 4. In C, we have the same, but put all pij in one
% long vector: p = [p11 p12 ... p33] and have pointers *pij indicating the
% memory address where each pij starts.
i = 0:4;
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

    format long g; result_exact = [s11 s12 s13 s21 s22 s23 s31 s32 s33].'
    */

    const UINT memsize = poly_fmult3x3_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    COMPLEX result_exact[27] = {
        -3.5031867497515 +      I*2.24765605596327,
          -4.64810742924058 +     I* 8.94837891435535,
          0.622199796943657 +     I*  9.5847779731096,
          -3.81072719638816 +     I* 1.78305829158921,
          -5.39732219478625 +     I* 9.63753237188223,
           1.30715160491115 +     I* 9.72763270802446,
          -4.08019211653745 +     I*  1.3006447981663,
          -6.09260870107913 +     I* 10.2303907916939,
           1.97904278612924 +     I* 9.77329215244861,
          -4.01094991242082 +     I* 1.11200537300603,
          -3.55813470718989 +     I*   6.187185459638,
           3.11769269143431 +     I* 8.72942108267116,
          -4.16745649563259 +     I*0.577273759586521,
          -4.34407997791979 +     I* 7.08424561504628,
           3.79110144206306 +     I* 8.68142453235659,
          -4.28232323111797 +    I*0.0367742175825609,
          -5.08662063747224 +     I* 7.91052233000927,
           4.42663076025425 +     I* 8.54668605781213,
          -4.16042686503002 -     I*0.122977438090491,
          -2.15032440876494 +     I* 2.87330915480378,
           5.33469138307117 +     I* 7.09429100532663,
          -4.15191931785452 -     I*0.680076918093759,
          -2.90279403438722 +     I* 3.89814429608133,
           5.93640347824795 +     I* 6.85973055867644,
          -4.10192716529707 -     I* 1.23038129433497,
          -3.62625990155404 +     I* 4.88403046810762,
           6.47880099218417 +     I* 6.55662995180941};
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
    printf("p = \n");
    for (UINT j=0; j<(sizeof(p)/sizeof(p[0])); j++){
        printf("%f + i%f\n", creal(p[j]), cimag(p[j]));
    }

    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult3x3(&deg, n, p, result, W_ptr);
        // Checking lengths of exact and numerical results
    UINT length_exact = sizeof(result_exact)/sizeof(result_exact[0]);
    UINT length_num = sizeof(result)/sizeof(result[0]);
    // Printing the results
    UINT j;
    printf("exact\n");
    for (j=0; j<length_exact; j++){
        printf("%f + i%f\n", creal(result_exact[j]), cimag(result_exact[j]));
    }
    printf("numerical\n");
    for (j=0; j<length_num; j++){
        printf("%f + i%f\n", creal(result[j]), cimag(result[j]));
    }
    printf("length exact = %u \n, lenght_num = %u\n", length_exact, length_num);
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

    /* Commenting out normalization, first check w/o normalization
    ret_code = poly_fmult3x3_test_n_is_power_of_2(1); // 3x3 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }
    */

    return EXIT_SUCCESS;
}
