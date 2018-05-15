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

// poly_fmult_p1xp2 for lazy people, takes care of all buffers (ineffecient
// since they are not reused => TODO)
static INT poly_fmult_p1xp2_X(
    const UINT deg, 
    COMPLEX const * const p1, 
    COMPLEX const * const p2,
    COMPLEX * const result,
    const INT add_flag)
{
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    COMPLEX * buf0 = NULL;
    COMPLEX * buf1 = NULL;
    COMPLEX * buf2 = NULL;
    INT ret_code = SUCCESS;

    const UINT len = poly_fmult_p1xp2_len(deg);

    buf0 = fft_wrapper_malloc(len*sizeof(COMPLEX));
    buf1 = fft_wrapper_malloc(len*sizeof(COMPLEX));
    buf2 = fft_wrapper_malloc(len*sizeof(COMPLEX));
    if (buf0 == NULL || buf1 == NULL || buf2 == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    ret_code = fft_wrapper_create_plan(&plan_fwd, len, buf0, buf1, -1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_create_plan(&plan_inv, len, buf0, buf1, 1);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = poly_fmult_p1xp2(deg, p1, p2, result, plan_fwd, plan_inv, buf0,
        buf1, buf2, add_flag);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:    
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    return ret_code;
}

static INT nse_finvscatter_recurse(
    const UINT deg,
    COMPLEX const * const T_11,
    COMPLEX const * const T_12,
    COMPLEX const * const T_21,
    COMPLEX const * const T_22,
    COMPLEX * const Ti_11,
    COMPLEX * const Ti_12,
    COMPLEX * const Ti_21,
    COMPLEX * const Ti_22,
    COMPLEX * const q,
    const REAL eps_t,
    const INT kappa,
    const nse_discretization_t discretization)
{
    INT ret_code = SUCCESS;
    COMPLEX * T1 = NULL;
    COMPLEX * T1i = NULL;
    COMPLEX * T2i = NULL;

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
        T1i = calloc(4*(deg/2 + 1), sizeof(COMPLEX));
        T2i = calloc(4*(deg + 1), sizeof(COMPLEX));
        if (T1 == NULL || T1i == NULL || T2i == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }

        COMPLEX * T1_11 = T1;
        COMPLEX * T1_12 = T1_11 + (2*deg + 1);
        COMPLEX * T1_21 = T1_12 + (2*deg + 1);
        COMPLEX * T1_22 = T1_21 + (2*deg + 1);

        COMPLEX * T2i_11 = T2i + deg/2;
        COMPLEX * T2i_12 = T2i_11 + (deg + 1);
        COMPLEX * T2i_21 = T2i_12 + (deg + 1);
        COMPLEX * T2i_22 = T2i_21 + (deg + 1);

        // Step 1: Determine T2i(z) and the samples q[D/2],...,q[D-1] from the
        // upper part of T(z).
        
        ret_code = nse_finvscatter_recurse(deg/2, T_11+deg/2, T_12+deg/2,
            T_21+deg/2, T_22+deg/2, T2i_11, T2i_12, T2i_21, T2i_22,
            q + (deg/2), eps_t, kappa, discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        T2i_11 = T2i;
        T2i_12 = T2i_11 + (deg + 1);
        T2i_21 = T2i_12 + (deg + 1);
        T2i_22 = T2i_21 + (deg + 1);

        // Step 2: Determine T1(z) = T2i(z)*T(z). Note that since the degree
        // of T2i(z) is deg/2 
        
        ret_code = poly_fmult_p1xp2_X(deg, T2i_11, T_11, T1_11, 0);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_12, T_21, T1_11, 1);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_11, T_12, T1_12, 0);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_12, T_22, T1_12, 1);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_21, T_11, T1_21, 0);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_22, T_21, T1_21, 1);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_21, T_12, T1_22, 0);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = poly_fmult_p1xp2_X(deg, T2i_22, T_22, T1_22, 1);
        CHECK_RETCODE(ret_code, leave_fun);

        T1_11 = T1 + deg;
        T1_12 = T1_11 + (2*deg + 1);
        T1_21 = T1_12 + (2*deg + 1);
        T1_22 = T1_21 + (2*deg + 1);

        T2i_11 = T2i + deg/2;
        T2i_12 = T2i_11 + (deg + 1);
        T2i_21 = T2i_12 + (deg + 1);
        T2i_22 = T2i_21 + (deg + 1);

        // Step 3: Determine T1i(z) and q[0],...,q[D/2-1] from T1(z)
 
        COMPLEX * T1i_11 = T1i;
        COMPLEX * T1i_12 = T1i_11 + (deg/2 + 1);
        COMPLEX * T1i_21 = T1i_12 + (deg/2 + 1);
        COMPLEX * T1i_22 = T1i_21 + (deg/2 + 1);
       
        ret_code = nse_finvscatter_recurse(deg/2, T1_11, T1_12,
            T1_21, T1_22, T1i_11, T1i_12, T1i_21, T1i_22,
            q, eps_t, kappa, discretization);
        CHECK_RETCODE(ret_code, leave_fun);
 
        // Step 4: Determine Ti(z) = T1i(z)*T2i(z) if requested
        if (Ti_11 != NULL) {
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_11, T2i_11, Ti_11, 0);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_12, T2i_21, Ti_11, 1);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_11, T2i_12, Ti_12, 0);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_12, T2i_22, Ti_12, 1);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_21, T2i_11, Ti_21, 0);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_22, T2i_21, Ti_21, 1);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_21, T2i_12, Ti_22, 0);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = poly_fmult_p1xp2_X(deg/2, T1i_22, T2i_22, Ti_22, 1);
            CHECK_RETCODE(ret_code, leave_fun);
        }

    } else if (deg == 1) { // base case

        const COMPLEX Q = -kappa*CONJ(T_21[1] / T_11[1]);
        *q = Q/eps_t;

        const REAL scl_den = SQRT( 1.0 + kappa*CABS(Q)*CABS(Q) );
        if (scl_den == 0.0)
            return E_DIV_BY_ZERO;
        const REAL scl = 1.0 / scl_den;

        Ti_11[0] = scl;
        Ti_11[1] = 0.0;
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

    COMPLEX const * const T_11 = transfer_matrix;
    COMPLEX const * const T_12 = T_11 + (deg+1);
    COMPLEX const * const T_21 = T_12 + (deg+1);
    COMPLEX const * const T_22 = T_21 + (deg+1);

    ret_code = nse_finvscatter_recurse(deg, T_11, T_12, T_21, T_22,
        NULL, NULL, NULL, NULL, q, eps_t, kappa, discretization);
    if (ret_code != SUCCESS)
        ret_code = E_SUBROUTINE(ret_code);

    return ret_code; 
}
