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

inline INT poly_fmult_two_polys_len(const UINT deg)
{
    return fft_wrapper_next_fft_length(2*(deg + 1) - 1);
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
    const INT mode)
{
    UINT i;
    INT ret_code = SUCCESS;

    // Zero-pad polynomials
    const UINT len = poly_fmult_two_polys_len(deg);
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

    } else if (mode == 3) {

        // Add product of FFT's to previously stored product in result to
        // current product
        for (i = 0; i < len; i++)
            buf0[i] += result[i];
    }

    if (mode != 2) {

        // Inverse FFT of product
        ret_code = fft_wrapper_execute_plan(plan_inv, buf0, buf1);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    if (mode == 0 || mode == 3) {

        // Store relevant part of scaled result of inverse FFT in result
        for (i = 0; i < 2*deg + 1; i++)
            result[i] = buf1[i]/len;

    } else if (mode == 1) {

        for (i = 0; i < 2*deg + 1; i++)
            result[i] += buf1[i]/len;
    }

leave_fun:
    return ret_code;
}

static inline INT poly_rescale(const UINT d, COMPLEX * const p)
{
    UINT i;
    INT a;
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

    return a;
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

    // Allocate memory for for calls to poly_fmult2
    lenmem = poly_fmult_two_polys_len(deg * n) * sizeof(COMPLEX);
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

static inline INT poly_rescale2x2(const UINT d,
    COMPLEX * const p11,
    COMPLEX * const p12,
    COMPLEX * const p21,
    COMPLEX * const p22)
{
    UINT i;
    INT a;
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

    return a;
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

    // Allocate memory for for calls to poly_fmult2
    lenmem = poly_fmult_two_polys_len(deg * n) * sizeof(COMPLEX);
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
