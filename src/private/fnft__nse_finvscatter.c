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
 * Sander Wahls (TU Delft) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__nse_finvscatter.h"
#include "fnft__nse_discretization.h"
#include "fnft__misc.h"
#include "fnft__poly_fmult.h"
#include "fnft__fft_wrapper.h"
#include <stdio.h>

static inline INT create_fft_plans(const UINT deg,
    fft_wrapper_plan_t * const plan_fwd_ptr,
    fft_wrapper_plan_t * const plan_inv_ptr,
    COMPLEX * const buf0, COMPLEX * const buf1)
{
    INT ret_code = SUCCESS;

    const UINT len = poly_fmult_two_polys_len(deg);
    ret_code = fft_wrapper_create_plan(plan_fwd_ptr, len, buf0, buf1, -1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_create_plan(plan_inv_ptr, len, buf0, buf1, 1);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}

static INT nse_finvscatter_recurse(
    const UINT deg,
    COMPLEX const * const T,
    const UINT T_stride,
    COMPLEX * const Ti,
    const UINT Ti_stride,
    COMPLEX * const q,
    const REAL eps_t,
    const INT kappa,
    const nse_discretization_t discretization,
    COMPLEX * const buf0,
    COMPLEX * const buf1,
    COMPLEX * const buf2)
{
    INT ret_code = SUCCESS;
    COMPLEX * T1 = NULL;
    COMPLEX * T1i = NULL;
    COMPLEX * T2i = NULL;
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
 
    if (deg > 1) { // recursive case

        // The transfer matrix for all samples is
        //
        //   T(z) = T2(z)*T1(z), where
        //
        // T1(z) is built from the samples
        //
        //   q[0],...,q[D/2-1]
        //
        // and T2(z) is build from
        //
        //   q[D/2],...,q[D-1],
        //
        // respectively. The inverse of T(z), up to some power of z, is
        //
        //   Ti(z) = T1i(z)*T2i(z),
        //
        // where Tni(z) is the inverse of Tn(z) [also up to some power of z].

        // Allocate memory, prepare pointers
        T1 = malloc(4*(2*deg + 1) * sizeof(COMPLEX));
        T1i = malloc(4*(deg/2 + 1) * sizeof(COMPLEX));
        T2i = calloc(4*(deg + 1), sizeof(COMPLEX));
        if (T1 == NULL || T1i == NULL || T2i == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }

        // Step 1: Determine T2i(z) and the samples q[D/2],...,q[D-1] from the
        // lower part of T(z).
        
        ret_code = nse_finvscatter_recurse(deg/2, T+deg/2, T_stride, T2i+deg/2,
            deg+1, q+(deg/2), eps_t, kappa, discretization, buf0, buf1, buf2);
        CHECK_RETCODE(ret_code, leave_fun);

        // Step 2: Determine T1(z) = T2i(z)*T(z).

        create_fft_plans(deg, &plan_fwd, &plan_inv, buf0, buf1);
        ret_code = poly_fmult_two_polys2x2(deg, T2i, deg+1, T, T_stride, T1,
            2*deg+1, plan_fwd, plan_inv, buf0, buf1, buf2);
        CHECK_RETCODE(ret_code, leave_fun);
        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_destroy_plan(&plan_inv);

        // Step 3: Determine T1i(z) and q[0],...,q[D/2-1] from T1(z).

        ret_code = nse_finvscatter_recurse(deg/2, T1+deg, 2*deg+1, T1i,
            deg/2+1, q, eps_t, kappa, discretization, buf0, buf1, buf2);
        CHECK_RETCODE(ret_code, leave_fun);
 
        // Step 4: Determine Ti(z) = T1i(z)*T2i(z) if requested
        if (Ti != NULL) {
            create_fft_plans(deg/2, &plan_fwd, &plan_inv, buf0, buf1);
            ret_code = poly_fmult_two_polys2x2(deg/2,
                T1i, deg/2+1,
                T2i+deg/2, deg+1,
                Ti, Ti_stride,
                plan_fwd, plan_inv, buf0, buf1, buf2);
            CHECK_RETCODE(ret_code, leave_fun);
            fft_wrapper_destroy_plan(&plan_fwd);
            fft_wrapper_destroy_plan(&plan_inv);
        }

    } else if (deg == 1) { // base case

        COMPLEX const * const T_21 = T + 2*T_stride;

        const COMPLEX Q = -kappa*CONJ(T_21[1] / T[1]);
        *q = Q/eps_t;

        const REAL scl_den = SQRT( 1.0 + kappa*CABS(Q)*CABS(Q) );
        if (scl_den == 0.0)
            return E_DIV_BY_ZERO;
        const REAL scl = 1.0 / scl_den;

        COMPLEX * const Ti_12 = Ti + Ti_stride;
        COMPLEX * const Ti_21 = Ti_12 + Ti_stride;
        COMPLEX * const Ti_22 = Ti_21 + Ti_stride;

        Ti[0] = scl;
        Ti[1] = 0.0;
        Ti_12[0] = -scl*Q;
        Ti_12[1] = 0.0;
        Ti_21[0] = 0.0;
        Ti_21[1] = scl*kappa*CONJ(Q);
        Ti_22[0] = 0.0;
        Ti_22[1] = scl;

    } else { // deg == 0

        ret_code = E_ASSERTION_FAILED;

    }

leave_fun:
    free(T1);
    free(T2i);
    return ret_code;
}

INT nse_finvscatter(
    const UINT deg,
    COMPLEX const * const transfer_matrix,
    COMPLEX * const q,
    const REAL eps_t,
    const INT kappa,
    const nse_discretization_t discretization)
{
    if (deg == 0)
        return E_INVALID_ARGUMENT(de);
    if (transfer_matrix == NULL)
        return E_INVALID_ARGUMENT(transfer_matrix);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (!(eps_t > 0.0))
        return E_INVALID_ARGUMENT(eps_t);
    if (kappa != -1 && kappa != 1)
        return E_INVALID_ARGUMENT(kappa);
    if (discretization != fnft_nse_discretization_2SPLIT2_MODAL)
        return E_INVALID_ARGUMENT(discretization);

    const UINT discretization_degree = nse_discretization_degree(
        discretization);
    if (discretization_degree == 0)
        return E_INVALID_ARGUMENT(discretization);
    const UINT D = deg / discretization_degree;
    if (D < 2 || (D&(D-1)) != 0) {
        return E_OTHER("Number of samples D used to build the transfer matrix was no positive power of two.");
    }

    INT ret_code = SUCCESS;

    const UINT max_len = poly_fmult_two_polys_len(deg);
    COMPLEX * const buf0 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));
    COMPLEX * const buf1 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));
    COMPLEX * const buf2 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));
    if (buf0 == NULL || buf1 == NULL || buf2 == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    ret_code = nse_finvscatter_recurse(deg, transfer_matrix, deg+1,
        NULL, 0, q, eps_t, kappa, discretization,
        buf0, buf1, buf2);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    return ret_code;
}
