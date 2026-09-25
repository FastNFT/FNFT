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
 * Sander Wahls (TU Delft) 2017-2018, 2021, 2023.
 * Peter J Prins (TU Delft) 2020.
 * Lianne de Vries (TU Delft student) 2021.
 * Igor Chekhovskoy (NSU, FRC ICT) 2026.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdbool.h>
#include <limits.h>
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

UINT poly_fmult3x3_numel(UINT deg, UINT n)
{
     return 9*(deg+1)*misc_nextpowerof2(n);
}

inline UINT poly_fmult_two_polys_len(const UINT deg)
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
    const UINT mode)
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

        // Sum product of FFT's with previously stored product in result and
        // save the result in buf0
        for (i = 0; i < len; i++)
            buf0[i] += result[i];

    } else if (mode == 4) {

        // Add product of FFT's to previously stored product in result
        for (i = 0; i < len; i++)
            result[i] += buf0[i];
    }

    if (mode != 2 && mode != 4) {

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
    return misc_normalize_vector(d+1, p);
}

INT fnft__poly_power_to_chebyshev(const UINT degree,
    COMPLEX const * const power, COMPLEX * const chebyshev)
{
    const UINT complex_max = ((UINT)-1)/sizeof(COMPLEX);
    COMPLEX *copy = NULL;
    COMPLEX const *source = power;
    REAL scale;
    UINT i, j;

    if (power == NULL || chebyshev == NULL)
        return E_INVALID_ARGUMENT(power);
    if (degree >= complex_max)
        return E_INVALID_ARGUMENT(degree);
    if (power == chebyshev) {
        copy = malloc((degree + 1)*sizeof(COMPLEX));
        if (copy == NULL)
            return E_NOMEM;
        memcpy(copy, power, (degree + 1)*sizeof(COMPLEX));
        source = copy;
    }
    memset(chebyshev, 0, (degree + 1)*sizeof(COMPLEX));
    chebyshev[0] = source[0];
    for (i = 1; i <= degree; i++) {
        REAL coefficient = 1.0;

        scale = ldexp(1.0, 1 - (INT)i);
        for (j = 0; j < (i + 1)/2; j++) {
            chebyshev[i - 2*j] += source[i]*scale*coefficient;
            coefficient *= (REAL)(i - j)/(REAL)(j + 1);
        }
        if ((i & 1U) == 0)
            chebyshev[0] += source[i]*0.5*scale*coefficient;
    }
    free(copy);
    return SUCCESS;
}

COMPLEX fnft__poly_eval_chebyshev(const UINT degree,
    COMPLEX const * const chebyshev, const COMPLEX x)
{
    COMPLEX b1 = 0.0, b2 = 0.0;
    UINT i;

    if (chebyshev == NULL)
        return NAN + I*NAN;
    for (i = degree; i > 0; i--) {
        const COMPLEX b0 = 2.0*x*b1 - b2 + chebyshev[i];
        b2 = b1;
        b1 = b0;
    }
    return x*b1 - b2 + chebyshev[0];
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
    const UINT sufficent_space_flag)
{
    UINT mode;
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

    if (sufficent_space_flag == 0)
        mode = 0;
    else
        mode = 2;

    // If sufficent_space_flag==0:
    //
    //   Store values of [aj ak al; dj dk dl; gj gk gl] in result_11, result_12, result_13... result_33
    //
    // If sufficent_space_flag!=0:
    //
    //   Store FFT's of [aj ak al; dj dk dl; gj gk gl] in result_11, result_12, result_13... result_33

    ret_code = poly_fmult_two_polys(deg, p2_11, p1_11, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_12, NULL /*p1_11*/, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_11*/, p2_13, result_13,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_13*/, p1_21, result_23,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_12, NULL /*p1_21*/, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_21*/, p2_11, result_21,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_11*/, p1_31, result_31,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_12, NULL /*p1_31*/, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_31*/, p2_13, result_33,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

    if (sufficent_space_flag == 0)
        mode = 1;
    else
        mode = 4;

    // add 2nd term
    ret_code = poly_fmult_two_polys(deg, p2_21, p1_12, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_22, NULL /*p1_12*/, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_12*/, p2_23, result_13,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_23*/, p1_22, result_23,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_22, NULL /*p1_22*/, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_22*/, p2_21, result_21,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_21*/, p1_32, result_31,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_22, NULL /*p1_32*/, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_32*/, p2_23, result_33,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);

    // add 3rd term
    if (sufficent_space_flag)
        mode = 3;

    ret_code = poly_fmult_two_polys(deg, p2_31, p1_13, result_11,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_32, NULL /*p1_13*/, result_12,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_13*/, p2_33, result_13,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_33*/, p1_23, result_23,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_32, NULL /*p1_23*/, result_22,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_23*/, p2_31, result_21,
                                    plan_fwd, plan_inv, buf0, buf2, buf1, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p2_31*/, p1_33, result_31,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, p2_32, NULL /*p1_33*/, result_32,
                                    plan_fwd, plan_inv, buf0, buf1, buf2, mode);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = poly_fmult_two_polys(deg, NULL /*p1_33*/, p2_33, result_33,
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

static void chebyshev_to_laurent(const UINT degree,
    COMPLEX const * const chebyshev, COMPLEX * const laurent)
{
    UINT i;

    memset(laurent, 0, (2*degree + 1)*sizeof(COMPLEX));
    laurent[degree] = chebyshev[0];
    for (i = 1; i <= degree; i++) {
        laurent[degree - i] = 0.5*chebyshev[i];
        laurent[degree + i] = 0.5*chebyshev[i];
    }
}

static void laurent_to_chebyshev(const UINT degree,
    COMPLEX const * const laurent, COMPLEX * const chebyshev)
{
    UINT i;

    chebyshev[0] = laurent[degree];
    for (i = 1; i <= degree; i++)
        chebyshev[i] = laurent[degree - i] + laurent[degree + i];
}

static INT multiply_two_chebyshev(const UINT degree,
    COMPLEX const * const first, COMPLEX const * const second,
    COMPLEX * const product, fft_wrapper_plan_t plan_fwd,
    fft_wrapper_plan_t plan_inv, COMPLEX * const buffer0,
    COMPLEX * const buffer1, COMPLEX * const buffer2,
    COMPLEX * const laurent_first, COMPLEX * const laurent_second,
    COMPLEX * const laurent_product)
{
    INT ret_code;

    chebyshev_to_laurent(degree, first, laurent_first);
    chebyshev_to_laurent(degree, second, laurent_second);
    ret_code = poly_fmult_two_polys(2*degree, laurent_first,
            laurent_second, laurent_product, plan_fwd, plan_inv, buffer0,
            buffer1, buffer2, 0);
    if (ret_code == SUCCESS)
        laurent_to_chebyshev(2*degree, laurent_product, product);
    return ret_code;
}

INT fnft__poly_fmult_chebyshev(UINT * const d, UINT n,
    COMPLEX * const p, INT * const W_ptr)
{
    const UINT uint_max = (UINT)-1;
    const UINT complex_max = uint_max/sizeof(COMPLEX);
    const UINT original_degree = d == NULL ? 0 : *d;
    UINT degree, excess, i, j, length, max_degree, padded_n;
    COMPLEX *first, *second, *product;
    COMPLEX *buffer0 = NULL, *buffer1 = NULL, *buffer2 = NULL;
    COMPLEX *laurent_first = NULL, *laurent_second = NULL;
    COMPLEX *laurent_product = NULL;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    INT W = 0, ret_code = SUCCESS;

    if (d == NULL || p == NULL || n == 0)
        return E_INVALID_ARGUMENT(p);
    if (original_degree == uint_max)
        return E_INVALID_ARGUMENT(*d);
    if (n == 1) {
        if (W_ptr != NULL)
            *W_ptr = 0;
        return SUCCESS;
    }
    if (n > uint_max/2 + 1U)
        return E_INVALID_ARGUMENT(n);
    padded_n = misc_nextpowerof2(n);
    if (padded_n == 0
            || original_degree + 1 > complex_max/padded_n)
        return E_INVALID_ARGUMENT(n);
    excess = padded_n - n;
    degree = original_degree;
    first = p + n*(degree + 1);
    for (i = 0; i < excess; i++) {
        first[0] = 1.0;
        for (j = 1; j <= degree; j++)
            first[j] = 0.0;
        first += degree + 1;
    }
    n = padded_n;
    max_degree = degree*n/2;
    if (max_degree > (complex_max - 1)/4)
        return E_INVALID_ARGUMENT(n);
    length = poly_fmult_two_polys_len(2*max_degree);
    if (length == 0 || length > complex_max)
        return E_INVALID_ARGUMENT(n);
    buffer0 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    buffer1 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    buffer2 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    laurent_first = malloc((2*max_degree + 1)*sizeof(COMPLEX));
    laurent_second = malloc((2*max_degree + 1)*sizeof(COMPLEX));
    laurent_product = malloc((4*max_degree + 1)*sizeof(COMPLEX));
    if (buffer0 == NULL || buffer1 == NULL || buffer2 == NULL
            || laurent_first == NULL || laurent_second == NULL
            || laurent_product == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    while (n >= 2) {
        length = poly_fmult_two_polys_len(2*degree);
        ret_code = fft_wrapper_create_plan(&plan_fwd, length, buffer0,
                buffer1, -1);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = fft_wrapper_create_plan(&plan_inv, length, buffer0,
                buffer1, 1);
        CHECK_RETCODE(ret_code, leave_fun);
        first = p;
        second = p + degree + 1;
        product = p;
        for (i = 0; i < n; i += 2) {
            ret_code = multiply_two_chebyshev(degree, first, second, product,
                    plan_fwd, plan_inv, buffer0, buffer1, buffer2,
                    laurent_first, laurent_second, laurent_product);
            CHECK_RETCODE(ret_code, leave_fun);
            if (W_ptr != NULL)
                W += poly_rescale(2*degree, product);
            first += 2*degree + 2;
            second += 2*degree + 2;
            product += 2*degree + 1;
        }
        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);
        degree *= 2;
        n /= 2;
    }
    *d = degree - excess*original_degree;
    if (W_ptr != NULL)
        *W_ptr = W;

leave_fun:
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buffer0);
    fft_wrapper_free(buffer1);
    fft_wrapper_free(buffer2);
    free(laurent_first);
    free(laurent_second);
    free(laurent_product);
    return ret_code;
}

static INT multiply_two_chebyshev2x2(const UINT degree,
    COMPLEX const * const first, const UINT first_stride,
    COMPLEX const * const second, const UINT second_stride,
    COMPLEX * const product, const UINT product_stride,
    fft_wrapper_plan_t plan_fwd, fft_wrapper_plan_t plan_inv,
    COMPLEX * const buffer0, COMPLEX * const buffer1,
    COMPLEX * const buffer2, COMPLEX * const laurent_first,
    COMPLEX * const laurent_second, COMPLEX * const laurent_product)
{
    const UINT input_stride = 2*degree + 1;
    const UINT output_stride = 4*degree + 1;
    UINT component;
    INT ret_code;

    for (component = 0; component < 4; component++) {
        chebyshev_to_laurent(degree, first + component*first_stride,
                laurent_first + component*input_stride);
        chebyshev_to_laurent(degree, second + component*second_stride,
                laurent_second + component*input_stride);
    }
    ret_code = poly_fmult_two_polys2x2(2*degree, laurent_first,
            input_stride, laurent_second, input_stride, laurent_product,
            output_stride, plan_fwd, plan_inv, buffer0, buffer1, buffer2, 0);
    if (ret_code != SUCCESS)
        return ret_code;
    for (component = 0; component < 4; component++)
        laurent_to_chebyshev(2*degree,
                laurent_product + component*output_stride,
                product + component*product_stride);
    return SUCCESS;
}

INT fnft__poly_fmult2x2_chebyshev(UINT * const d, UINT n,
    COMPLEX * const p, COMPLEX * const result, INT * const W_ptr)
{
    const UINT uint_max = (UINT)-1;
    const UINT complex_max = uint_max/sizeof(COMPLEX);
    const UINT original_degree = d == NULL ? 0 : *d;
    UINT degree, excess, padded_n, max_degree, length, i, j;
    UINT first_offset, second_offset, product_offset;
    UINT input_stride, output_stride;
    COMPLEX *p11, *p12, *p21, *p22;
    COMPLEX *buffer0 = NULL, *buffer1 = NULL, *buffer2 = NULL;
    COMPLEX *laurent_first = NULL, *laurent_second = NULL;
    COMPLEX *laurent_product = NULL;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    INT W = 0, ret_code = SUCCESS;

    if (d == NULL || p == NULL || result == NULL || n == 0)
        return E_INVALID_ARGUMENT(p);
    if (original_degree == uint_max)
        return E_INVALID_ARGUMENT(*d);
    if (n == 1) {
        if (original_degree + 1 > complex_max/4)
            return E_INVALID_ARGUMENT(*d);
        memcpy(result, p, 4*(original_degree + 1)*sizeof(COMPLEX));
        if (W_ptr != NULL)
            *W_ptr = 0;
        return SUCCESS;
    }
    if (n > uint_max/2 + 1U)
        return E_INVALID_ARGUMENT(n);
    padded_n = misc_nextpowerof2(n);
    if (padded_n == 0 || padded_n > complex_max/4
            || original_degree + 1 > complex_max/(4*padded_n))
        return E_INVALID_ARGUMENT(n);
    excess = padded_n - n;
    degree = original_degree;
    p11 = p;
    p12 = p11 + n*(degree + 1);
    p21 = p12 + n*(degree + 1);
    p22 = p21 + n*(degree + 1);
    if (excess > 0) {
        COMPLEX * const p12_padded = p + padded_n*(degree + 1);
        COMPLEX * const p21_padded = p12_padded + padded_n*(degree + 1);
        COMPLEX * const p22_padded = p21_padded + padded_n*(degree + 1);

        memmove(p22_padded, p22, n*(degree + 1)*sizeof(COMPLEX));
        memmove(p21_padded, p21, n*(degree + 1)*sizeof(COMPLEX));
        memmove(p12_padded, p12, n*(degree + 1)*sizeof(COMPLEX));
        p12 = p12_padded;
        p21 = p21_padded;
        p22 = p22_padded;
        for (i = n; i < padded_n; i++) {
            p11[i*(degree + 1)] = 1.0;
            p12[i*(degree + 1)] = 0.0;
            p21[i*(degree + 1)] = 0.0;
            p22[i*(degree + 1)] = 1.0;
            for (j = 1; j <= degree; j++) {
                p11[i*(degree + 1) + j] = 0.0;
                p12[i*(degree + 1) + j] = 0.0;
                p21[i*(degree + 1) + j] = 0.0;
                p22[i*(degree + 1) + j] = 0.0;
            }
        }
    }
    n = padded_n;
    max_degree = degree*n/2;
    if (max_degree > (complex_max/4 - 1)/4) {
        ret_code = E_INVALID_ARGUMENT(n);
        goto leave_fun;
    }
    length = poly_fmult_two_polys_len(2*max_degree);
    if (length == 0 || length > complex_max) {
        ret_code = E_INVALID_ARGUMENT(n);
        goto leave_fun;
    }
    buffer0 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    buffer1 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    buffer2 = fft_wrapper_malloc(length*sizeof(COMPLEX));
    laurent_first = malloc(4*(2*max_degree + 1)*sizeof(COMPLEX));
    laurent_second = malloc(4*(2*max_degree + 1)*sizeof(COMPLEX));
    laurent_product = malloc(4*(4*max_degree + 1)*sizeof(COMPLEX));
    if (buffer0 == NULL || buffer1 == NULL || buffer2 == NULL
            || laurent_first == NULL || laurent_second == NULL
            || laurent_product == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    while (n >= 2) {
        length = poly_fmult_two_polys_len(2*degree);
        ret_code = fft_wrapper_create_plan(&plan_fwd, length, buffer0,
                buffer1, -1);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = fft_wrapper_create_plan(&plan_inv, length, buffer0,
                buffer1, 1);
        CHECK_RETCODE(ret_code, leave_fun);
        input_stride = n*(degree + 1);
        output_stride = (n/2)*(2*degree + 1);
        first_offset = 0;
        second_offset = degree + 1;
        product_offset = 0;
        for (i = 0; i < n; i += 2) {
            ret_code = multiply_two_chebyshev2x2(degree,
                    p + first_offset, input_stride, p + second_offset,
                    input_stride, result + product_offset, output_stride,
                    plan_fwd, plan_inv, buffer0, buffer1, buffer2,
                    laurent_first, laurent_second, laurent_product);
            CHECK_RETCODE(ret_code, leave_fun);
            if (W_ptr != NULL)
                W += poly_rescale2x2(2*degree, result + product_offset,
                        result + output_stride + product_offset,
                        result + 2*output_stride + product_offset,
                        result + 3*output_stride + product_offset);
            first_offset += 2*degree + 2;
            second_offset += 2*degree + 2;
            product_offset += 2*degree + 1;
        }
        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);
        degree *= 2;
        n /= 2;
        if (n > 1) {
            const UINT count = n*(degree + 1);
            memcpy(p, result, count*sizeof(COMPLEX));
            memcpy(p + count, result + count, count*sizeof(COMPLEX));
            memcpy(p + 2*count, result + 2*count, count*sizeof(COMPLEX));
            memcpy(p + 3*count, result + 3*count, count*sizeof(COMPLEX));
        }
    }
    if (excess > 0) {
        const UINT true_degree = degree - excess*original_degree;
        const UINT padded_stride = degree + 1;
        const UINT true_stride = true_degree + 1;

        memmove(result + true_stride, result + padded_stride,
                true_stride*sizeof(COMPLEX));
        memmove(result + 2*true_stride, result + 2*padded_stride,
                true_stride*sizeof(COMPLEX));
        memmove(result + 3*true_stride, result + 3*padded_stride,
                true_stride*sizeof(COMPLEX));
        degree = true_degree;
    }
    *d = degree;
    if (W_ptr != NULL)
        *W_ptr = W;

leave_fun:
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buffer0);
    fft_wrapper_free(buffer1);
    fft_wrapper_free(buffer2);
    free(laurent_first);
    free(laurent_second);
    free(laurent_product);
    return ret_code;
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
        const UINT r_stride = (n/2)*(2*deg+1);

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

            const UINT sufficent_space_flag = r_stride - or < len ? 0 : 1;
            ret_code = poly_fmult_two_polys3x3(deg, p+o1, p_stride, p+o2,
                                               p_stride, result+or, r_stride,
                                               plan_fwd, plan_inv, buf0, buf1,
                                               buf2, sufficent_space_flag);
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
