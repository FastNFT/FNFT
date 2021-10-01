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
* Lianne de Vries (TU Delft student) 2021.
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
        return 9*(deg+1)*misc_nextpowerof2(n);
}

inline UINT poly_fmult_two_polys_len(const UINT deg)
{
    return 2*(deg + 1) - 1;
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
    UINT len = 0;

    // Zero-pad polynomials
        len = 2*(deg + 1) - 1;                   // This line works for the 3x3 case as well as 2x2 case

    memset(&buf0[deg+1], 0, (len - (deg+1))*sizeof(COMPLEX));

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

    } else if (mode == 3 || mode == 5) {

        // Add product of FFT's to previously stored product in result to
        // current product
        for (i = 0; i < len; i++)
            buf0[i] += result[i];
    } else if (mode == 4){
        // Add product of FFT's to previously stored product in result to
        // current product, but store in result so we can add the last term to result later
        for (i = 0; i < len; i++)
            result[i] += buf0[i];
    }

    if (mode != 2 && mode != 4 && mode != 5) {

        // Inverse FFT of product
        ret_code = fft_wrapper_execute_plan(plan_inv, buf0, buf1);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    if (mode == 5) {

        // Inverse FFT of product, but store in buf0 so we can re-use buf1
        memcpy(result, buf0, len*sizeof(COMPLEX));
        ret_code = fft_wrapper_execute_plan(plan_inv, result, buf0);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    if (mode == 0 || mode == 3) {

        // Store relevant part of scaled result of inverse FFT in result
        for (i = 0; i < 2*deg + 1; i++)
            result[i] = buf1[i]/len;

    } else if (mode == 1) {

        for (i = 0; i < 2*deg + 1; i++)
            result[i] += buf1[i]/len;
    } else if (mode == 5) {
        // Store relevant part of scaled result of inverse FFT in result
        for (i = 0; i < 2*deg + 1; i++)
            result[i] = buf0[i]/len;
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


// 3x3 case
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
    
    COMPLEX const * const p2_12 = p2_11 + p2_stride;
    COMPLEX const * const p2_13 = p2_12 + p2_stride;
    COMPLEX const * const p2_21 = p2_13 + p2_stride;
    COMPLEX const * const p2_22 = p2_21 + p2_stride;
    COMPLEX const * const p2_23 = p2_22 + p2_stride;
    COMPLEX const * const p2_31 = p2_23 + p2_stride;
    COMPLEX const * const p2_32 = p2_31 + p2_stride;
    COMPLEX const * const p2_33 = p2_32 + p2_stride;

    COMPLEX * const result_12 = result_11 + result_stride;
    COMPLEX * const result_13 = result_12 + result_stride;
    COMPLEX * const result_21 = result_13 + result_stride;
    COMPLEX * const result_22 = result_21 + result_stride;
    COMPLEX * const result_23 = result_22 + result_stride;
    COMPLEX * const result_31 = result_23 + result_stride;
    COMPLEX * const result_32 = result_31 + result_stride;
    COMPLEX * const result_33 = result_32 + result_stride;

    // We compute the matrix product
    //
    //  [a b c; d e f; g h i][j k l; m n o; p q r] = [aj ak al; dj dk dl; gj gk gl]+
    //      [bm bn bo; em en eo; hm hn ho] + [cp cq cr; fp fq fr; ip iq ir],
    // where a=p1_11, b=p1_12, etc.

    UINT mode = mode_offset;

    // If mode_offset==0:
    //
    //   Store values of [aj ak al; dj dk dl; gj gk gl] in result_11, result_12, result_13... result_33
    //
    // If mode_offset==2:
    //
    //   Store FFT's of [aj ak al; dj dk dl; gj gk gl] in result_11, result_12, result_13... result_33

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

    mode = mode_offset + 2;     // We are now switching to the mode:
                                // not calculating 1st term, but also not calculating last term,
                                // so add result to previous result but don't take inverse fft yet

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

    // add 3rd term
    mode = mode + 1;        // We are now switching to the mode:
                            // calculating last term. Add to previous result AND take inverse fft
                            // (mode 1 or 3)

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


// 3x3 multiplication
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
            // We pad with z^deg*[1 0 0; 0 1 0; 0 0 1]
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
        UINT r_stride_temp = (n/2)*(2*deg+1);
        if (r_stride_temp < len)
            r_stride_temp = len;
        const UINT r_stride = r_stride_temp;

        r11 = result;
        r12 = r11 + r_stride;
        r13 = r12 + r_stride;
        r21 = r13 + r_stride;
        r22 = r21 + r_stride;
        r23 = r22 + r_stride;
        r31 = r23 + r_stride;
        r32 = r31 + r_stride;
        r33 = r32 + r_stride;

        // Multiply all pairs of polynomials, normalize if desired
        for (i=0; i<n; i+=2) {

            const UINT mode_offset  = r_stride - or < len ? 0 : 2;
            ret_code = poly_fmult_two_polys3x3(deg, p+o1, p_stride, p+o2,
                                               p_stride, result+or, r_stride,
                                               plan_fwd, plan_inv, buf0, buf1,
                                               buf2, mode_offset);
            CHECK_RETCODE(ret_code, release_mem);

            // Normalize if desired
            if (W_ptr != NULL){
                W += poly_rescale3x3(2*deg, r11+or, r12+or, r13+or, r21+or, r22+or,
                                     r23+or, r31+or, r32+or, r33+or);
            }

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
            memcpy(p13, r13, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p21, r21, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p22, r22, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p23, r23, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p31, r31, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p32, r32, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p33, r33, n*(deg+1)*sizeof(COMPLEX));
        }
    }

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


