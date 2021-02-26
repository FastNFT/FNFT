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

 /* This code checks multiplication for scalar real 3x3 matrices 
*/ 

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include <stdio.h>

static INT poly_fmult3x3_test_n_is_power_of_2(INT normalize_flag)
{
    printf("test1 \n");
    UINT deg = 0, n = 2;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:

P1 = [1 2 3; 4 5 6; 7 8 9];
P2 = [10 11 12; 13 14 15; 16 17 18];

P1_reshape = reshape(P1',[1,9]);
P2_reshape = reshape(P2',[1,9]);
P = zeros(1,18);
for i = 1:9
    P(2*i-1) = P1_reshape(i);
    P(2*i) = P2_reshape(i);
end

result = P1*P2
    */

    const UINT memsize = poly_fmult3x3_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    printf("memsize = %d\n", memsize);
    COMPLEX result_exact[9] = {
        84.000000,
    90.000000,
    96.000000,
    201.000000,
    216.000000,
    231.000000,
    318.000000,
    342.000000,
    366.000000};
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/9 - 1;

    p[0] = 1;
    p[1] = 10;
    p[2] = 2;
    p[3] = 11;
    p[4] = 3;
    p[5] = 12;
    p[6] = 4;
    p[7] = 13;
    p[8] = 5;
    p[9] = 14;
    p[10] = 6;
    p[11] = 15;
    p[12] = 7;
    p[13] = 16;
    p[14] = 8;
    p[15] = 17;
    p[16] = 9;
    p[17] = 18;
    
        for (UINT j=0; j<(sizeof(p)/sizeof(p[0])); j++){
        printf("%f + i%f\n", creal(p[j]), cimag(p[j]));
    }

    if (normalize_flag)
        W_ptr = &W;
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
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
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
        printf("Compile test 3\n");
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult3x3_test_n_is_power_of_2(1); // 3x3 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        printf("Compile test 4\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS; // Replaced SUCCES FAILURE
}
