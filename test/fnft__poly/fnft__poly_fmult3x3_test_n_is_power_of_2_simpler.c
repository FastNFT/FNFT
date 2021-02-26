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

 /* This code checks a simpler case than the original code: for now, 
   we multiply 2 scalar (degree 0) matrices
*/ 

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include <stdio.h>

static INT poly_fmult3x3_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 0, n = 2;
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
i = 0:1;
    p11 = sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
    p12 = sqrt(i+1).*(cos(i+0.1) + 1j*sin(-2*i+0.1));
    p13 = sqrt(i+1).*(cos(i+0.2) + 1j*sin(-2*i+0.2));
    p21 = sqrt(i+1).*(cos(i+0.3) + 1j*sin(-2*i+0.3));
    p22 = sqrt(i+1).*(cos(i+0.4) + 1j*sin(-2*i+0.4));
    p23 = sqrt(i+1).*(cos(i+0.5) + 1j*sin(-2*i+0.5));
    p31 = sqrt(i+1).*(cos(i+0.6) + 1j*sin(-2*i+0.6));
    p32 = sqrt(i+1).*(cos(i+0.7) + 1j*sin(-2*i+0.7));
    p33 = sqrt(i+1).*(cos(i+0.8) + 1j*sin(-2*i+0.8));
    
% So now we have the 2 scalar matrices and multiply them:
P1 = [p11(1) p12(1) p13(1);
      p21(1) p22(1) p23(1);
      p31(1) p32(1) p33(1)];
  
P2 = [p11(2) p12(2) p13(2);
      p21(2) p22(2) p23(2);
      p31(2) p32(2) p33(2)];

  result = P1*P2;
  
  format long g; result_exact = result.'
    */

    const UINT memsize = poly_fmult3x3_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    COMPLEX result_exact[9] = {
        124.004332775875 -       I*53.828612030162,
           88.0683598392946 -     I* 98.0806107862452,
            51.252436962729 -     I* 141.352620500492,
           134.373306461915 -     I* 14.7786511825712,
           113.119720060944 -     I* 67.6740064597929,
           90.7358788095321 -     I* 119.893185434521,
           132.739112879174 +     I*  25.591442560623,
           128.066432588492 -     I* 31.2222846864791,
           122.114154836521 -     I* 87.7240491857091};

    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/9 - 1;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = 10*(SQRT(i+1.0)*(COS(i) + I*SIN(-3.3*i)));
        p[i+n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.1) + I*SIN(-3.3*i+0.1)));
        p[i+2*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.2) + I*SIN(-3.3*i+0.2)));
        p[i+3*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.3) + I*SIN(-3.3*i+0.3)));
        p[i+4*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.4) + I*SIN(-3.3*i+0.4)));
        p[i+5*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.5) + I*SIN(-3.3*i+0.5)));
        p[i+6*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.6) + I*SIN(-3.3*i+0.6)));
        p[i+7*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.7) + I*SIN(-3.3*i+0.7)));
        p[i+8*n*(deg+1)] = 10*(SQRT(i+1.0)*(COS(i+0.8) + I*SIN(-3.3*i+0.8)));
    }
    printf("i= (should be 0) %d\n", i);
        for (UINT j=0; j<(sizeof(p)/sizeof(p[0])); j++){
        printf("%f + i%f\n", creal(p[j]), cimag(p[j]));
    }

    if (normalize_flag)
        W_ptr = &W;
        printf("W = %d\n", W);
        printf("W_ptr = %d\n", W_ptr);
    ret_code = poly_fmult3x3(&deg, n, p, result, W_ptr);
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
    printf("ret_code = %d\n", ret_code);
    if (ret_code != SUCCESS){
        return E_SUBROUTINE(ret_code);
        printf("compile test 3a\n");
    }
    if (deg != deg_exact){
        printf("deg test 3= %llu\n", deg);
        printf("deg_exact = %llu\n", deg_exact);
        return E_TEST_FAILED;
    }
    if (normalize_flag) {
        printf("test for W = %d", W);
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        printf("scl = %f\n", scl);
        for (i=0; i<9*(deg+1); i++)
            result[i] *= scl;
    }
    if (!(misc_rel_err(9*(deg+1), result, result_exact) <= 100*EPSILON)){
        REAL err;
        err = misc_rel_err(9*(deg+1), result, result_exact);
        printf("relative error = %d\n", err);
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
\
    return EXIT_SUCCESS;
}
