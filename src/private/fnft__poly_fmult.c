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
 * Peter J Prins (TU Delft) 2020.
 * Lianne de Vries (TU Delft) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__fft_wrapper.h"

#ifdef HAVE_PRAGMA_GCC_OPTIMIZE_OFAST
#pragma GCC optimize("Ofast")
#endif

UINT poly_fmult_numel(UINT deg, UINT n)
{
    return (deg+1)*misc_nextpowerof2(n);
}

UINT poly_fmult2x2_numel(UINT deg, UINT n)
{
    return 4*(deg+1)*misc_nextpowerof2(n);
}

// 3x3 case
UINT poly_fmult3x3_numel(UINT deg, UINT n)
{

	return 9*(deg+1)*misc_nextpowerof2(n);      // TODO: maybe change this to include fft_next_size, then we ensure "result" is big enough to set r_stride=len
}

/*
// General case, NxN is the dimension of the AKNS system
UINT poly_fmultNxN_numel(UINT deg, UINT n)
{
    return (N*N)*(deg+1)*misc_nextpowerof2(n);
}
*/

inline UINT poly_fmult_two_polys_len(const UINT deg)
{
    //return 2*(deg + 1) - 1;
    return fft_wrapper_next_fft_length(2*(deg + 1) - 1);
    /* Original code:
    return fft_wrapper_next_fft_length(2*(deg + 1) - 1);
    calls a function from the kiss_fft library and outputs the "next fast size"
    which is the next number that can be written completely in factors of 2, 3 and 5
    This is useful because the kiss_fft library calculates the fft more efficiently
    in that case. However, this caused problems when the lenght of buf0 determined
    by this function was >result_stride
    TODO: set result_stride to appropriate length if this happens
    */

}

inline INT poly_fmult_two_polys(
    const UINT deg,
    COMPLEX const * const p1,
    COMPLEX const * const p2,
    COMPLEX * const result,
    fft_wrapper_plan_t plan_fwd,
    fft_wrapper_plan_t plan_inv,
    COMPLEX * const buf0,
    COMPLEX * const buf1,
    COMPLEX * const buf2,
    const UINT mode)
{

    UINT i;
    INT ret_code = SUCCESS;

    // Zero-pad polynomials
    const UINT len = 2*(deg + 1) - 1;       // TODO: use next_fast_size
    memset(&buf0[deg+1], 0, (len - (deg+1))*sizeof(COMPLEX));

    printf("len in fmult_two_polys = %ld\n", len);
    printf("deg in fmult_two_polys = %ld\n", deg);
    printf("mode = %ld\n",mode);

    // FFT of first polynomial
    if (p1 != NULL) {
        memcpy(buf0, p1, (deg+1)*sizeof(COMPLEX));
        ret_code = fft_wrapper_execute_plan(plan_fwd, buf0, buf1);
        CHECK_RETCODE(ret_code, leave_fun);
    }
 
    // FFT of second polynomial
    if (p2 != NULL) {
        memcpy(buf0, p2, (deg+1)*sizeof(COMPLEX));
        ret_code = fft_wrapper_execute_plan(plan_fwd, buf0, buf2);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Multiply FFT's
    for (i = 0; i < len; i++)
        buf0[i] = buf1[i] * buf2[i];

    if (mode == 2) {


        // Temporarily store product of FFT's in result
        memcpy(result, buf0, len*sizeof(COMPLEX));

    } else if (mode == 3) {

        // Add product of FFT's to previously stored product in result to
        // current product
        for (i = 0; i < len; i++)
            result[i] += buf0[i];

    } else if (mode == 4) {

        // Add product of FFT's to previously stored product in result to
        // current product
        for (i = 0; i < len; i++)
            result[i] += buf0[i];

        // Temporarily store the fft of the final answer (result) in buf0 (not needed?)
        // memcpy(buf0, result, len*sizeof(COMPLEX));
    }


    if (mode != 2 && mode !=4) {
    
        // Inverse FFT of product
        ret_code = fft_wrapper_execute_plan(plan_inv, result, buf0);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    if (mode == 0 || mode == 3) {

        // Store relevant part of scaled result of inverse FFT in result
        for (i = 0; i < 2*deg + 1; i++)
            result[i] = buf0[i]/len;

    } else if (mode == 1) {

        for (i = 0; i < 2*deg + 1; i++)
            result[i] += buf0[i]/len;
    }
    
    if (mode ==3){
    }

leave_fun:
    return ret_code;
}

static inline INT poly_rescale(const UINT d, COMPLEX * const p)
{
    UINT i;
    REAL a;
    REAL scl;
    REAL cur_abs;
    REAL max_abs = 0.0;

    // Find max of absolute values of coefficients
    max_abs = 0.0;
    for (i=0; i<=d; i++) {
        cur_abs = CABS( p[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
    }

    // Return if polynomial is identical to zero
    if (max_abs == 0.0)
        return 0;

    // Otherwise, rescale
    a = FLOOR( LOG2(max_abs) );
    scl = POW( 2.0, -a );
    for (i=0; i<=d; i++)
        p[i] *= scl;

    return (INT) a;
}

INT fnft__poly_fmult(UINT * const d, UINT n, COMPLEX * const p,
    INT * const W_ptr)
{
    UINT i, j, deg, len, lenmem;
    COMPLEX *p1, *p2, *result;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    COMPLEX *buf0 = NULL, *buf1 = NULL, *buf2 = NULL;
    INT W = 0;
    INT ret_code;

    // Pad with z^deg if n is not a power of two
    const UINT n_excess = misc_nextpowerof2(n) - n;
    deg = *d;
    p1 = p + n*(deg + 1);
    for (i = 0; !(i >= n_excess); i++) { // "<" does not work because of UINT
        p1[0] = 1.0;
        for (j = 1; j<=deg; j++)
            p1[j] = 0.0;
        p1 += deg + 1;
    }
    n += n_excess;

    // Allocate memory for calls to poly_fmult_two_polys
    lenmem = poly_fmult_two_polys_len(deg * n/2) * sizeof(COMPLEX);
    buf0 = fft_wrapper_malloc(lenmem);
    buf1 = fft_wrapper_malloc(lenmem);
    buf2 = fft_wrapper_malloc(lenmem);
    if (buf0 == NULL || buf1 == NULL || buf2 == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Main loop, n is the current number of polynomials, deg is their degree
    while (n >= 2) {

        // Create FFT and IFFT config (computes twiddle factors, so reuse)
        len = poly_fmult_two_polys_len(deg);
        ret_code = fft_wrapper_create_plan(&plan_fwd, len, buf0, buf1, -1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = fft_wrapper_create_plan(&plan_inv, len, buf0, buf1, 1);
        CHECK_RETCODE(ret_code, release_mem);

        // Pointers to current pair of polynomials and their product
        p1 = p;
        p2 = p + (deg + 1);
        result = p;

        // Multiply all pairs of polynomials, normalize if desired
        for (i=0; i<n; i+=2) {
            ret_code = poly_fmult_two_polys(deg, p1, p2, result, plan_fwd,
                plan_inv, buf0, buf1, buf2, 0);
            CHECK_RETCODE(ret_code, release_mem);

            if (W_ptr != NULL)
                W += poly_rescale(2*deg, result);

            p1 += 2*deg + 2;
            p2 += 2*deg + 2;
            result += 2*deg + 1;
        }

        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);

        // Double degrees and half the number of polynomials
        deg *= 2;
        if (n%2 != 0) { // n was no power of two
            ret_code = E_INVALID_ARGUMENT(n);
            goto release_mem;
        }
        n /= 2;
    }

    // Set degree of final result, free memory and return w/o error
    *d = deg - n_excess*(*d);
    if (W_ptr != NULL)
        *W_ptr = W;
release_mem:
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    return ret_code;
}

inline INT poly_fmult_two_polys2x2(const UINT deg,
    COMPLEX const * const p1_11,
    const UINT p1_stride,
    COMPLEX const * const p2_11,
    const UINT p2_stride,
    COMPLEX * const result_11,
    const UINT result_stride,
    fft_wrapper_plan_t plan_fwd,
    fft_wrapper_plan_t plan_inv,
    COMPLEX * const buf0,
    COMPLEX * const buf1,
    COMPLEX * const buf2,
    const UINT mode_offset)
{
    INT ret_code;

    COMPLEX const * const p1_12 = p1_11 + p1_stride;
    COMPLEX const * const p1_21 = p1_12 + p1_stride;
    COMPLEX const * const p1_22 = p1_21 + p1_stride;

    COMPLEX const * const p2_12 = p2_11 + p2_stride;
    COMPLEX const * const p2_21 = p2_12 + p2_stride;
    COMPLEX const * const p2_22 = p2_21 + p2_stride;

    COMPLEX * const result_12 = result_11 + result_stride;
    COMPLEX * const result_21 = result_12 + result_stride;
    COMPLEX * const result_22 = result_21 + result_stride;

    // We compute the matrix product
    //
    //  [a b ; c d][e f ; g h]=[ae af ; ce cf]+[bg bh ; dg dh],
    //
    // where a=p1_11, b=p1_12, etc.

    UINT mode = mode_offset;

    // If mode_offset==0:
    //
    //   Store ae, af, ce and cf in result_11, result_12, result_21 and
    //   result_22, respectively
    //
    // If mode_offset==2:
    //
    //   Store FFT's of ae, af, ce and cf in result_11, result_12, result_21 and
    //   result_22, respectively. This saves us performing some inverse FFT's,
    //   but since the polynomials are zero-padded before the FFT is taken, it
    //   is only possible if the result_* arrays are larger than strictly
    //   necessary. See the doc of this fun.

    ret_code = poly_fmult_two_polys(deg, p1_11, p2_11, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_11*/, p1_21, result_21,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_21*/, p2_12, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_12*/, p1_11, result_12,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

    mode = mode_offset + 1;

    // If mode_offset==0:
    //
    //   Add bg, bh, dg and dg to result_11, result_12, result_21 and
    //   result_22, respectively
    //
    // If mode_offset==2:
    //
    //   Add result_11, result_12, result_21 and result_22 to the FFT's of ae,
    //   af, ce and cf, respectively, and perform inverse FFT's.

    ret_code = poly_fmult_two_polys(deg, p1_12, p2_21, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_21*/, p1_22, result_21,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_22*/, p2_22, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_22*/, p1_12, result_12,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}


// My code, adjusted for 3x3 case
inline INT poly_fmult_two_polys3x3(const UINT deg,
    COMPLEX const * const p1_11,
    const UINT p1_stride,
    COMPLEX const * const p2_11,
    const UINT p2_stride,
    COMPLEX * const result_11,
    const UINT result_stride,
    fft_wrapper_plan_t plan_fwd,
    fft_wrapper_plan_t plan_inv,
    COMPLEX * const buf0,
    COMPLEX * const buf1,
    COMPLEX * const buf2,
    const UINT mode_offset)
{

    INT ret_code;

    COMPLEX const * const p1_12 = p1_11 + p1_stride;
    COMPLEX const * const p1_13 = p1_12 + p1_stride;
    COMPLEX const * const p1_21 = p1_13 + p1_stride;
    COMPLEX const * const p1_22 = p1_21 + p1_stride;
    COMPLEX const * const p1_23 = p1_22 + p1_stride;
    COMPLEX const * const p1_31 = p1_23 + p1_stride;
    COMPLEX const * const p1_32 = p1_31 + p1_stride;
    COMPLEX const * const p1_33 = p1_32 + p1_stride;

    
    // Printing some stuff to see if the pointers point to the right addresses
    // Print first entry of each pij
/*    printf("First value that pij points to. Should be coefficient of 1st matrix for 6th order element\n");
        printf("%f + i%f\n", creal(*p1_11), cimag(*p1_11));
        printf("%f + i%f\n", creal(*p1_12), cimag(*p1_12));
        printf("%f + i%f\n", creal(*p1_13), cimag(*p1_13));
        printf("%f + i%f\n", creal(*p1_21), cimag(*p1_21));
        printf("%f + i%f\n", creal(*p1_22), cimag(*p1_22));
        printf("%f + i%f\n", creal(*p1_23), cimag(*p1_23));
        printf("%f + i%f\n", creal(*p1_31), cimag(*p1_31));
        printf("%f + i%f\n", creal(*p1_32), cimag(*p1_32));
        printf("%f + i%f\n", creal(*p1_33), cimag(*p1_33));*/

    // Printing coefficients of first matrix to be multiplied
/*   printf("coeffs of first matrix in fmult_two_polys3x3\n");
   UINT j;
    
    printf("p1_11 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_11[j]), cimag(p1_11[j]));

    printf("p1_12 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_12[j]), cimag(p1_12[j]));

    printf("p1_13 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_13[j]), cimag(p1_13[j]));
    
    printf("p1_21 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_21[j]), cimag(p1_21[j]));

    printf("p1_22 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_22[j]), cimag(p1_22[j]));

    printf("p1_23 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_23[j]), cimag(p1_23[j]));

    printf("p1_31 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_31[j]), cimag(p1_31[j]));

    printf("p1_32 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_32[j]), cimag(p1_32[j]));

    printf("p1_33 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p1_33[j]), cimag(p1_33[j]));*/
    
    

    
    COMPLEX const * const p2_12 = p2_11 + p2_stride;
    COMPLEX const * const p2_13 = p2_12 + p2_stride;
    COMPLEX const * const p2_21 = p2_13 + p2_stride;
    COMPLEX const * const p2_22 = p2_21 + p2_stride;
    COMPLEX const * const p2_23 = p2_22 + p2_stride;
    COMPLEX const * const p2_31 = p2_23 + p2_stride;
    COMPLEX const * const p2_32 = p2_31 + p2_stride;
    COMPLEX const * const p2_33 = p2_32 + p2_stride;

/*   printf("coeffs of second matrix in fmult_two_polys3x3\n");
    
    printf("p2_11 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_11[j]), cimag(p2_11[j]));

    printf("p2_12 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_12[j]), cimag(p2_12[j]));

    printf("p2_13 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_13[j]), cimag(p2_13[j]));
    
    printf("p2_21 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_21[j]), cimag(p2_21[j]));

    printf("p2_22 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_22[j]), cimag(p2_22[j]));

    printf("p2_23 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_23[j]), cimag(p2_23[j]));

    printf("p2_31 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_31[j]), cimag(p2_31[j]));

    printf("p2_32 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_32[j]), cimag(p2_32[j]));

    printf("p2_33 = \n");
    for (j = 0; j<7; j++)
        printf("%f + i%f\n", creal(p2_33[j]), cimag(p2_33[j]));*/

    COMPLEX * const result_12 = result_11 + result_stride;
    COMPLEX * const result_13 = result_12 + result_stride;
    COMPLEX * const result_21 = result_13 + result_stride;
    COMPLEX * const result_22 = result_21 + result_stride;
    COMPLEX * const result_23 = result_22 + result_stride;
    COMPLEX * const result_31 = result_23 + result_stride;
    COMPLEX * const result_32 = result_31 + result_stride;
    COMPLEX * const result_33 = result_32 + result_stride;
    
    printf("result_stride=%ld\n",result_stride);

    // Printing results
/*    printf("result BEFORE multipication in fmult_two_polys3x3\n");
    UINT j;
    printf("result_11 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_11[j]), cimag(result_11[j]));

    printf("result_12 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_12[j]), cimag(result_12[j]));

    printf("result_13 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_13[j]), cimag(result_13[j]));
    
    printf("result_21 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_21[j]), cimag(result_21[j]));

    printf("result_22 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_22[j]), cimag(result_22[j]));

    printf("result_23 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_23[j]), cimag(result_23[j]));

    printf("result_31 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_31[j]), cimag(result_31[j]));

    printf("result_32 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_32[j]), cimag(result_32[j]));

    printf("result_33 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_33[j]), cimag(result_33[j]));

    printf("deg in polys3x3 = %d\n",deg);*/

    // Second time, these memory addresses already contain values, but they shouldn't (I think)
    // However, this should not matter? They are overwritten in first round of multiplication when "mode"=2?

    // We compute the matrix product
    //
    //  [a b ; c d][e f ; g h]=[ae af ; ce cf]+[bg bh ; dg dh],
    //
    // where a=p1_11, b=p1_12, etc.

    UINT mode = mode_offset;

    // If mode_offset==0:
    //
    //   Store ae, af, ce and cf in result_11, result_12, result_21 and
    //   result_22, respectively
    //
    // If mode_offset==2:
    //
    //   Store FFT's of ae, af, ce and cf in result_11, result_12, result_21 and
    //   result_22, respectively. This saves us performing some inverse FFT's,
    //   but since the polynomials are zero-padded before the FFT is taken, it
    //   is only possible if the result_* arrays are larger than strictly
    //   necessary. See the doc of this fun.

    ret_code = poly_fmult_two_polys(deg, p1_11, p2_11, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_11*/, p2_12, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_11*/, p2_13, result_13,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_13*/, p1_21, result_23,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_21*/, p2_12, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_21*/, p2_11, result_21,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_11*/, p1_31, result_31,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_31*/, p2_12, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_31*/, p2_13, result_33,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

/*    printf("result_11 after calculating first term of multiplication= \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_11[j]), cimag(result_11[j]));*/


    mode = mode_offset + 2;     // We are now switching to the mode:
                                // not calculating 1st term, but also not calculating last term,
                                // so add result to previous result but don't take inverse fft yet

    // If mode_offset==0:
    //
    //   Add bg, bh, dg and dg to result_11, result_12, result_21 and
    //   result_22, respectively
    //
    // If mode_offset==2:
    //
    //   Add result_11, result_12, result_21 and result_22 to the FFT's of ae,
    //   af, ce and cf, respectively, and perform inverse FFT's.

    // add 2nd term
    ret_code = poly_fmult_two_polys(deg, p1_12, p2_21, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_12*/, p2_22, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_12*/, p2_23, result_13,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_23*/, p1_22, result_23,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_22*/, p2_22, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_22*/, p2_21, result_21,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_21*/, p1_32, result_31,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_32*/, p2_22, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_32*/, p2_23, result_33,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

/*    printf("result_11 after calculating 2nd term of multiplication= \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_11[j]), cimag(result_11[j]));*/


    // add 3rd term
    mode = mode - 1;        // We are now switching to the mode:
                            // caluclating last term. Add to previous result AND take inverse fft
                            // (mode 1 or 3)

    // For 3x3, it is not impossible to avoid calculating some values twice,
    // unless we use an extra buffer. I chose not to, so in the last two calls
    // of poly_fmult_two_polys we cannot use the NULL pointer and reuse values
        ret_code = poly_fmult_two_polys(deg, p1_13, p2_31, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    
        ret_code = poly_fmult_two_polys(deg, NULL /*p1_13*/, p2_32, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p1_13*/, p2_33, result_13,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p2_33*/, p1_23, result_23,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p1_23*/, p2_32, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p1_23*/, p2_31, result_21,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p2_31*/, p1_33, result_31,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p1_33*/, p2_32, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_fmult_two_polys(deg, NULL /*p1_33*/, p2_33, result_33,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);


    // Printing results
/*    UINT j;
    printf("result AFTER multiplication in fmult_two_polys3x3\n");
    printf("result_11 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_11[j]), cimag(result_11[j]));

    printf("result_12 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_12[j]), cimag(result_12[j]));

    printf("result_13 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_13[j]), cimag(result_13[j]));
    
    printf("result_21 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_21[j]), cimag(result_21[j]));

    printf("result_22 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_22[j]), cimag(result_22[j]));

    printf("result_23 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_23[j]), cimag(result_23[j]));

    printf("result_31 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_31[j]), cimag(result_31[j]));

    printf("result_32 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_32[j]), cimag(result_32[j]));

    printf("result_33 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(result_33[j]), cimag(result_33[j]));*/


leave_fun:
    return ret_code;
}



static inline INT poly_rescale2x2(const UINT d,
    COMPLEX * const p11,
    COMPLEX * const p12,
    COMPLEX * const p21,
    COMPLEX * const p22)
{
    UINT i;
    REAL a;
    REAL scl;
    REAL cur_abs;
    REAL max_abs = 0.0;

    // Find max of absolute values of coefficients
    max_abs = 0.0;
    for (i=0; i<=d; i++) {
        cur_abs = CABS( p11[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p12[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p21[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p22[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
    }

    // Return if polynomials are all identical to zero
    if (max_abs == 0.0)
        return 0;

    // Otherwise, rescale
    a = FLOOR( LOG2(max_abs) );
    scl = POW( 2.0, -a );
    for (i=0; i<=d; i++) {
        p11[i] *= scl;
        p12[i] *= scl;
        p21[i] *= scl;
        p22[i] *= scl;
    }

    return (INT) a;
}

// rescale function for 3x3.
static inline INT poly_rescale3x3(const UINT d,
    COMPLEX * const p11,
    COMPLEX * const p12,
    COMPLEX * const p13,
    COMPLEX * const p21,
    COMPLEX * const p22,
    COMPLEX * const p23,
    COMPLEX * const p31,
    COMPLEX * const p32,
    COMPLEX * const p33)
{
    UINT i;
    REAL a;
    REAL scl;
    REAL cur_abs;
    REAL max_abs = 0.0;

    // Find max of absolute values of coefficients
    max_abs = 0.0;
    for (i=0; i<=d; i++) {
        cur_abs = CABS( p11[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p12[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p13[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p21[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p22[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p23[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p31[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p32[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p33[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
    }

    // Return if polynomials are all identical to zero
    if (max_abs == 0.0)
        return 0;

    // Otherwise, rescale
    a = FLOOR( LOG2(max_abs) );
    scl = POW( 2.0, -a );
    for (i=0; i<=d; i++) {
        p11[i] *= scl;
        p12[i] *= scl;
        p13[i] *= scl;
        p21[i] *= scl;
        p22[i] *= scl;
        p23[i] *= scl;
        p31[i] *= scl;
        p32[i] *= scl;
        p33[i] *= scl;
    }

    return (INT) a;
}

/*
* length of p = m*m*n*(deg+1)
* length of result = m*m*(n/2)*(2*deg+1)
* WARNING: p is overwritten
*/
INT fnft__poly_fmult2x2(UINT * const d, UINT n, COMPLEX * const p,
    COMPLEX * const result, INT * const W_ptr)
{
    UINT i, j, deg, lenmem, len;
    UINT o1, o2, or; // pointer offsets
    COMPLEX *p11, *p12, *p21, *p22;
    COMPLEX *p11_pad, *p12_pad, *p21_pad, *p22_pad;
    COMPLEX *r11 = NULL, *r12 = NULL, *r21 = NULL, *r22 = NULL;
    COMPLEX *r12_pad, *r21_pad, *r22_pad;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    COMPLEX *buf0 = NULL, *buf1 = NULL, *buf2 = NULL;
    INT W = 0;
    INT ret_code;

    // Setup pointers to the individual polynomials in p
    deg = *d;
    p11 = p;
    p12 = p11 + n*(deg+1);
    p21 = p12 + n*(deg+1);
    p22 = p21 + n*(deg+1);

    // Pad if n is not a power of two
    const UINT n_excess = misc_nextpowerof2(n) - n;
    if (n_excess > 0) {

        // Pointers to beginning of polynomials after padding
        p11_pad = p;
        p12_pad = p11_pad + (n+n_excess)*(deg+1);
        p21_pad = p12_pad + (n+n_excess)*(deg+1);
        p22_pad = p21_pad + (n+n_excess)*(deg+1);

        memmove(p22_pad, p22, n*(deg+1)*sizeof(COMPLEX));
        memmove(p21_pad, p21, n*(deg+1)*sizeof(COMPLEX));
        memmove(p12_pad, p12, n*(deg+1)*sizeof(COMPLEX));

        // Reuse orig points as buffers to new poly's that are padded
        p11 = p11_pad + n*(deg + 1);
        p12 = p12_pad + n*(deg + 1);
        p21 = p21_pad + n*(deg + 1);
        p22 = p22_pad + n*(deg + 1);
        for (i = 0; !(i >= n_excess); i++) {
            // We pad with z^deg*[1 0; 0 1]
            p11[0] = 1.0;
            p12[0] = 0.0;
            p21[0] = 0.0;
            p22[0] = 1.0;
            for (j = 1; j<=deg; j++) {
                p11[j] = 0.0;
                p12[j] = 0.0;
                p21[j] = 0.0;
                p22[j] = 0.0;
            }
            p11 += deg + 1;
            p12 += deg + 1;
            p21 += deg + 1;
            p22 += deg + 1;
        }

        p11 = p11_pad;
        p12 = p12_pad;
        p21 = p21_pad;
        p22 = p22_pad;
        n += n_excess;
    }
    
    // Allocate memory for calls to poly_fmult_two_polys2x2
    lenmem = poly_fmult_two_polys_len(deg * n/2) * sizeof(COMPLEX);
    buf0 = fft_wrapper_malloc(lenmem);
    buf1 = fft_wrapper_malloc(lenmem);
    buf2 = fft_wrapper_malloc(lenmem);
    if (buf0 == NULL || buf1 == NULL || buf2 == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    const UINT p_stride = n*(deg + 1);

    // Main loop, n is the current number of polynomials, deg is their degree
    while (n >= 2) {

        // Create FFT and IFFT config (computes twiddle factors, so reuse)
        len = poly_fmult_two_polys_len(deg);
        ret_code = fft_wrapper_create_plan(&plan_fwd, len, buf0, buf1, -1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = fft_wrapper_create_plan(&plan_inv, len, buf0, buf2, 1);
        CHECK_RETCODE(ret_code, release_mem);

        // Offsets for the current pair of polynomials and their product
        o1 = 0;
        o2 = deg + 1;
        or = 0;

        // Setup pointers to the individual polynomials in result
        const UINT r_stride = (n/2)*(2*deg+1);
        r11 = result;
        r12 = r11 + r_stride;
        r21 = r12 + r_stride;
        r22 = r21 + r_stride;

        // Multiply all pairs of polynomials, normalize if desired
        for (i=0; i<n; i+=2) {

            const UINT mode_offset  = r_stride - or < len ? 0 : 2;
            ret_code = poly_fmult_two_polys2x2(deg, p+o1, p_stride, p+o2,
                                               p_stride, result+or, r_stride,
                                               plan_fwd, plan_inv, buf0, buf1,
                                               buf2, mode_offset);
            CHECK_RETCODE(ret_code, release_mem);

            // Normalize if desired
            if (W_ptr != NULL)
                W += poly_rescale2x2(2*deg, r11+or, r12+or, r21+or, r22+or);

            // Move pointers to next pair
            o1 += 2*deg + 2;
            o2 += 2*deg + 2;
            or += 2*deg + 1;
        }

        // Update degrees and number of polynomials
        deg *= 2;
        if (n%2 != 0) {
            ret_code = E_INVALID_ARGUMENT(n); // n was no power of two
            goto release_mem;
        }
        n /= 2;

        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);

        // Prepare for the next iteration
        if (n>1) {
            memcpy(p11, r11, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p12, r12, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p21, r21, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p22, r22, n*(deg+1)*sizeof(COMPLEX));
        }
    }

    // If padding was applied, reduce degree of the result
    if (n_excess > 0) {
        r12_pad = r12;
        r21_pad = r21;
        r22_pad = r22;
        deg -= n_excess*(*d);
        r12 = r11 + (deg+1);
        r21 = r12 + (deg+1);
        r22 = r21 + (deg+1);
        memmove(r12, r12_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r21, r21_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r22, r22_pad, (deg+1)*sizeof(COMPLEX));
    }

    // Set degree of final result, free memory and return w/o error
    *d = deg;
    if (W_ptr != NULL)
        *W_ptr = W;
release_mem:
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    return ret_code;
}


// My code, for 3x3 case %
INT fnft__poly_fmult3x3(UINT * const d, UINT n, COMPLEX * const p,
    COMPLEX * const result, INT * const W_ptr)
{
    UINT i, j, deg, lenmem, len;
    UINT o1, o2, or; // pointer offsets
    COMPLEX *p11, *p12, *p13, *p21, *p22, *p23, *p31, *p32, *p33;
    COMPLEX *p11_pad, *p12_pad, *p13_pad, *p21_pad, *p22_pad,
            *p23_pad, *p31_pad, *p32_pad, *p33_pad;
    COMPLEX *r11 = NULL, *r12 = NULL, *r13 = NULL, *r21 = NULL,
        *r22 = NULL, *r23 = NULL, *r31 = NULL, *r32 = NULL, *r33 = NULL;
    COMPLEX *r12_pad, *r13_pad, *r21_pad, *r22_pad, *r23_pad,
            *r31_pad, *r32_pad, *r33_pad;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    COMPLEX *buf0 = NULL, *buf1 = NULL, *buf2 = NULL;
    INT W = 0;
    INT ret_code;

    // Setup pointers to the individual polynomials in p
    deg = *d;
    printf("deg=%ld\n",deg);
    p11 = p;
    p12 = p11 + n*(deg+1);
    p13 = p12 + n*(deg+1);
    p21 = p13 + n*(deg+1);
    p22 = p21 + n*(deg+1);
    p23 = p22 + n*(deg+1);
    p31 = p23 + n*(deg+1);
    p32 = p31 + n*(deg+1);
    p33 = p32 + n*(deg+1);

    // Pad if n is not a power of two
    const UINT n_excess = misc_nextpowerof2(n) - n;
    if (n_excess > 0) {

        // Pointers to beginning of polynomials after padding
        p11_pad = p;
        p12_pad = p11_pad + (n+n_excess)*(deg+1);
        p13_pad = p12_pad + (n+n_excess)*(deg+1);
        p21_pad = p13_pad + (n+n_excess)*(deg+1);
        p22_pad = p21_pad + (n+n_excess)*(deg+1);
        p23_pad = p22_pad + (n+n_excess)*(deg+1);
        p31_pad = p23_pad + (n+n_excess)*(deg+1);
        p32_pad = p31_pad + (n+n_excess)*(deg+1);
        p33_pad = p32_pad + (n+n_excess)*(deg+1);

        memmove(p33_pad, p33, n*(deg+1)*sizeof(COMPLEX));
        memmove(p32_pad, p32, n*(deg+1)*sizeof(COMPLEX));
        memmove(p31_pad, p31, n*(deg+1)*sizeof(COMPLEX));
        memmove(p23_pad, p23, n*(deg+1)*sizeof(COMPLEX));
        memmove(p22_pad, p22, n*(deg+1)*sizeof(COMPLEX));
        memmove(p21_pad, p21, n*(deg+1)*sizeof(COMPLEX));
        memmove(p13_pad, p13, n*(deg+1)*sizeof(COMPLEX));
        memmove(p12_pad, p12, n*(deg+1)*sizeof(COMPLEX));

        // Reuse orig points as buffers to new poly's that are padded
        p11 = p11_pad + n*(deg + 1);
        p12 = p12_pad + n*(deg + 1);
        p13 = p13_pad + n*(deg + 1);
        p21 = p21_pad + n*(deg + 1);
        p22 = p22_pad + n*(deg + 1);
        p23 = p23_pad + n*(deg + 1);
        p31 = p31_pad + n*(deg + 1);
        p32 = p32_pad + n*(deg + 1);
        p33 = p33_pad + n*(deg + 1);

        for (i = 0; !(i >= n_excess); i++) {
            // We pad with z^deg*[1 0; 0 1]
            p11[0] = 1.0;
            p12[0] = 0.0;
            p13[0] = 0.0;
            p21[0] = 0.0;
            p22[0] = 1.0;
            p23[0] = 0.0;
            p31[0] = 0.0;
            p32[0] = 0.0;
            p33[0] = 1.0;
            for (j = 1; j<=deg; j++) {
                p11[j] = 0.0;
                p12[j] = 0.0;
                p13[j] = 0.0;
                p21[j] = 0.0;
                p22[j] = 0.0;
                p23[j] = 0.0;
                p31[j] = 0.0;
                p32[j] = 0.0;
                p33[j] = 0.0;
            }
            p11 += deg + 1;
            p12 += deg + 1;
            p13 += deg + 1;
            p21 += deg + 1;
            p22 += deg + 1;
            p23 += deg + 1;
            p31 += deg + 1;
            p32 += deg + 1;
            p33 += deg + 1;
        }

        p11 = p11_pad;
        p12 = p12_pad;
        p13 = p13_pad;
        p21 = p21_pad;
        p22 = p22_pad;
        p23 = p23_pad;
        p31 = p31_pad;
        p32 = p32_pad;
        p33 = p33_pad;
        n += n_excess;
    }
    
    // Allocate memory for calls to poly_fmult_two_polys3x3
    lenmem = poly_fmult_two_polys_len(deg * n/2) * sizeof(COMPLEX); //TODO: change how buf0 is allocated here?
    printf("len buf0= %ld\n",poly_fmult_two_polys_len(deg * n/2));
    buf0 = fft_wrapper_malloc(lenmem);
    buf1 = fft_wrapper_malloc(lenmem);
    buf2 = fft_wrapper_malloc(lenmem);
    if (buf0 == NULL || buf1 == NULL || buf2 == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    const UINT p_stride = n*(deg + 1);

    // Main loop, n is the current number of polynomials, deg is their degree
    while (n >= 2) {

        // Create FFT and IFFT config (computes twiddle factors, so reuse)
        len = 2*(deg + 1) - 1;    // TODO: original: poly_fmult_two_polys_len(deg);
        printf("len= %ld\n",poly_fmult_two_polys_len(deg));
        ret_code = fft_wrapper_create_plan(&plan_fwd, len, buf0, buf1, -1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = fft_wrapper_create_plan(&plan_inv, len, buf0, buf2, 1);
        CHECK_RETCODE(ret_code, release_mem);

        // Offsets for the current pair of polynomials and their product
        o1 = 0;
        o2 = deg + 1;
        or = 0;

        // Setup pointers to the individual polynomials in result
        const UINT r_stride = (n/2)*(2*deg+1);  // TODO: adjust using next_fast_size
        printf("r_stride = %d\n",r_stride);
        printf("len = %d\n",len);
        r11 = result;
        r12 = r11 + r_stride;
        r13 = r12 + r_stride;
        r21 = r13 + r_stride;
        r22 = r21 + r_stride;
        r23 = r22 + r_stride;
        r31 = r23 + r_stride;
        r32 = r31 + r_stride;
        r33 = r32 + r_stride;


        // Printing results
/*    printf("result before starting treewise multiplication\n");
    printf("r11 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r11[j]), cimag(r11[j]));

    printf("r12 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r12[j]), cimag(r12[j]));

    printf("r13 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r13[j]), cimag(r13[j]));
    
    printf("r21 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r21[j]), cimag(r21[j]));

    printf("r22 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r22[j]), cimag(r22[j]));

    printf("r23 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r23[j]), cimag(r23[j]));

    printf("r31 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r31[j]), cimag(r31[j]));

    printf("r32 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r32[j]), cimag(r32[j]));

    printf("r33 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r33[j]), cimag(r33[j]));*/


        // Multiply all pairs of polynomials, normalize if desired
        for (i=0; i<n; i+=2) {

            const UINT mode_offset  = r_stride - or < len ? 0 : 2;
            ret_code = poly_fmult_two_polys3x3(deg, p+o1, p_stride, p+o2,
                                               p_stride, result+or, r_stride,
                                               plan_fwd, plan_inv, buf0, buf1,
                                               buf2, mode_offset);
            CHECK_RETCODE(ret_code, release_mem);

            // Normalize if desired
            if (W_ptr != NULL)
                W += poly_rescale3x3(2*deg, r11+or, r12+or, r13+or, r21+or, r22+or,
                                     r23+or, r31+or, r32+or, r33+or);

            // Move pointers to next pair
            o1 += 2*deg + 2;
            o2 += 2*deg + 2;
            or += 2*deg + 1;

            printf("In loop %ld for n=%ld\n",i,n);
        }

        // Update degrees and number of polynomials
        deg *= 2;
        if (n%2 != 0) {
            ret_code = E_INVALID_ARGUMENT(n); // n was no power of two
            goto release_mem;
        }
        n /= 2;

        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);

       
        // Prepare for the next iteration
        if (n>1) {
            memcpy(p11, r11, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p12, r12, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p13, r13, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p21, r21, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p22, r22, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p23, r23, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p31, r31, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p32, r32, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p33, r33, n*(deg+1)*sizeof(COMPLEX));
        }
    }

    printf("result after treewise multiplication\n");
    printf("r11 = \n");
    for (UINT j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r11[j]), cimag(r11[j]));

    printf("r12 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r12[j]), cimag(r12[j]));

    printf("r13 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r13[j]), cimag(r13[j]));
    
    printf("r21 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r21[j]), cimag(r21[j]));

    printf("r22 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r22[j]), cimag(r22[j]));

    printf("r23 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r23[j]), cimag(r23[j]));

    printf("r31 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r31[j]), cimag(r31[j]));

    printf("r32 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r32[j]), cimag(r32[j]));

    printf("r33 = \n");
    for (j = 0; j<26; j++)
        printf("%f + i%f\n", creal(r33[j]), cimag(r33[j]));



    misc_print_buf(15*9,r11,"r");
    // If padding was applied, reduce degree of the result
    if (n_excess > 0) {
        r12_pad = r12;
        r13_pad = r13;
        r21_pad = r21;
        r22_pad = r22;
        r23_pad = r23;
        r31_pad = r31;
        r32_pad = r32;
        r33_pad = r33;
        deg -= n_excess*(*d);
        r12 = r11 + (deg+1);
        r13 = r12 + (deg+1);
        r21 = r13 + (deg+1);
        r22 = r21 + (deg+1);
        r23 = r22 + (deg+1);
        r31 = r23 + (deg+1);
        r32 = r31 + (deg+1);
        r33 = r32 + (deg+1);
        memmove(r12, r12_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r13, r13_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r21, r21_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r22, r22_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r23, r23_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r31, r31_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r32, r32_pad, (deg+1)*sizeof(COMPLEX));
        memmove(r33, r33_pad, (deg+1)*sizeof(COMPLEX));
    }

    // Set degree of final result, free memory and return w/o error
    *d = deg;
    printf("new deg = %ld\n",deg);
    if (W_ptr != NULL)
        *W_ptr = W;
release_mem:
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    return ret_code;
}


