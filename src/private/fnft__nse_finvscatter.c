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

struct fnft__nse_finvscatter_stack_elem {
    UINT deg;
    COMPLEX * T;
    UINT T_stride;
    COMPLEX * Ti;
    UINT Ti_stride;
    COMPLEX * q;
    COMPLEX * T1;
    COMPLEX * T1i;
    COMPLEX * T2i;
    fft_wrapper_plan_t plan_fwd_1;
    fft_wrapper_plan_t plan_inv_1;
    fft_wrapper_plan_t plan_fwd_2;
    fft_wrapper_plan_t plan_inv_2;
    INT ret_code;
    UINT start_pos;
};

// This is an iterative implementation of a recursive algorithm. We use our own
// custom stack to simulate what the routine would do when implemented using
// straight-forward recursion. This lets us avoid stack overflows that can
// happen for (very?) large degrees if recursion is implemented naively.
static INT nse_finvscatter_recurse(
    struct fnft__nse_finvscatter_stack_elem * const s,
    const REAL eps_t,
    const INT kappa,
    const nse_discretization_t discretization,
    COMPLEX * const buf0,
    COMPLEX * const buf1,
    COMPLEX * const buf2)
{
    INT ret_code;
    INT i = 0;

    while (i >= 0) {

        if (s[i].start_pos == 1)
            goto start_pos_1;
        else if (s[i].start_pos == 2)
            goto start_pos_2;
    
        if (s[i].deg > 1) { // recursive case
    
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
            // where Tni(z) is the inverse of Tn(z), also up to a power of z.
    
            // Allocate memory for level i if not done before
            if (s[i].T1 == NULL) {
                s[i].T1 = malloc(4*(2*s[i].deg + 1) * sizeof(COMPLEX));
                s[i].T1i = malloc(4*(s[i].deg/2 + 1) * sizeof(COMPLEX));
                s[i].T2i = calloc(4*(s[i].deg + 1), sizeof(COMPLEX));
                if (s[i].T1 == NULL || s[i].T1i == NULL || s[i].T2i == NULL) {
                    ret_code = E_NOMEM;
                    goto leave_fun;
                }
            }
    
            // Step 1: Determine T2i(z) and the samples q[D/2],...,q[D-1] from
            // the lower part of T(z) with a recursive call.
            
            s[i+1].deg = s[i].deg/2;
            s[i+1].T = s[i].T + s[i].deg/2;
            s[i+1].T_stride = s[i].T_stride;
            s[i+1].Ti = s[i].T2i + s[i].deg/2;
            s[i+1].Ti_stride = s[i].deg+1;
            s[i+1].q = s[i].q + s[i].deg/2;
            s[i+1].start_pos = 0;
            s[i].start_pos = 1;
            i++;
            continue;

start_pos_1:

            // Step 2: Determine T1(z) = T2i(z)*T(z).
   
            // Allocate 1st pair of FFT plans for level i if not done before
            if (s[i].plan_fwd_1 == fft_wrapper_safe_plan_init()) {
                create_fft_plans(s[i].deg, &s[i].plan_fwd_1, &s[i].plan_inv_1,
                    buf0, buf1);
            }

            ret_code = poly_fmult_two_polys2x2(s[i].deg,
                s[i].T2i, s[i].deg+1,
                s[i].T, s[i].T_stride,
                s[i].T1, 2*s[i].deg+1,
                s[i].plan_fwd_1, s[i].plan_inv_1, buf0, buf1, buf2);
            CHECK_RETCODE(ret_code, leave_fun);
    
            // Step 3: Determine T1i(z) and q[0],...,q[D/2-1] from T1(z) with
            // a recursive call.
    
            s[i+1].deg = s[i].deg/2;
            s[i+1].T = s[i].T1 + s[i].deg;
            s[i+1].T_stride = 2*s[i].deg + 1;
            s[i+1].Ti = s[i].T1i;
            s[i+1].Ti_stride = s[i].deg/2+1;
            s[i+1].q = s[i].q;
            s[i+1].start_pos = 0;
            s[i].start_pos = 2;
            i++;
            continue;

start_pos_2:

            // Step 4: Determine Ti(z) = T1i(z)*T2i(z) if requested
            
            if (s[i].Ti != NULL) {

                // Allocate 2nd pair of FFT plans for level i if not done
                // before
                if (s[i].plan_fwd_2 == fft_wrapper_safe_plan_init()) {
                    create_fft_plans(s[i].deg/2, &s[i].plan_fwd_2,
                        &s[i].plan_inv_2, buf0, buf1);
                }

                ret_code = poly_fmult_two_polys2x2(s[i].deg/2,
                    s[i].T1i, s[i].deg/2+1,
                    s[i].T2i+s[i].deg/2, s[i].deg+1,
                    s[i].Ti, s[i].Ti_stride,
                    s[i].plan_fwd_2, s[i].plan_inv_2, buf0, buf1, buf2);
                CHECK_RETCODE(ret_code, leave_fun);
            }
    
        } else if (s[i].deg == 1) { // base case
    
            COMPLEX const * const T_12 = s[i].T + s[i].T_stride;
            COMPLEX const * const T_22 = T_12 + 2*s[i].T_stride;

            const COMPLEX Q = (T_12[1] / T_22[1]);
            const REAL scl_den = 1.0 + kappa*CABS(Q)*CABS(Q);
            if (scl_den <= 0.0) {
                ret_code = E_OTHER("The step size eps_t is not small enough. Use more samples.");
                goto leave_fun;
            }
            const REAL scl = 1.0 / SQRT(scl_den);

            COMPLEX * const Ti_11 = s[i].Ti;
            COMPLEX * const Ti_12 = Ti_11 + s[i].Ti_stride;
            COMPLEX * const Ti_21 = Ti_12 + s[i].Ti_stride;
            COMPLEX * const Ti_22 = Ti_21 + s[i].Ti_stride;

            s[i].Ti[0] = 0.0;
            s[i].Ti[1] = scl;
            Ti_12[0] = 0.0;
            Ti_12[1] = -scl*Q;
            Ti_21[0] = scl*kappa*CONJ(Q);
            Ti_21[1] = 0.0;
            Ti_22[0] = scl;
            Ti_22[1] = 0.0;

            if (discretization == fnft_nse_discretization_2SPLIT2_MODAL) {

                *s[i].q = Q/eps_t;

            } else if (discretization == fnft_nse_discretization_2SPLIT2A) {

                const REAL Qabs = CABS(Q);
                if (Qabs >= FNFT_PI) {
                    ret_code = E_OTHER("The step size eps_t is not small enough. Use more samples.");
                    goto leave_fun;
                }
                *s[i].q = atan(Qabs)*CEXP(I*CARG(Q))/eps_t;

            } else { // discretization is unknown or not supported

                ret_code = E_INVALID_ARGUMENT(discretization);
                goto leave_fun;

            }
    
        } else { // deg == 0
    
            ret_code = E_ASSERTION_FAILED;
            goto leave_fun;
    
        }
    
        i--;
    }

leave_fun:
    return ret_code;
}

INT nse_finvscatter(
    const UINT deg,
    COMPLEX * const transfer_matrix,
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

    const UINT discretization_degree = nse_discretization_degree(
        discretization);
    if (discretization_degree == 0) // unknown discretization
        return E_INVALID_ARGUMENT(discretization);
    const UINT D = deg / discretization_degree;
    if (D < 2 || (D&(D-1)) != 0) {
        return E_OTHER("Number of samples D used to build the transfer matrix was no positive power of two.");
    }

    INT ret_code = SUCCESS;
    UINT i;

    const UINT max_len = poly_fmult_two_polys_len(deg);
    COMPLEX * const buf0 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));
    COMPLEX * const buf1 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));
    COMPLEX * const buf2 = fft_wrapper_malloc(max_len*sizeof(COMPLEX));

    const UINT stack_size = LOG2(D) + 1;
    struct fnft__nse_finvscatter_stack_elem * const s = malloc(
        stack_size * sizeof(struct fnft__nse_finvscatter_stack_elem));

    if (buf0 == NULL || buf1 == NULL || buf2 == NULL || s == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    s[0].deg = deg;
    s[0].T = transfer_matrix;
    s[0].T_stride = deg+1;
    s[0].Ti = NULL;
    s[0].Ti_stride = 0;
    s[0].q = q;
    s[0].start_pos = 0;
    s[0].ret_code = SUCCESS;

    for (i=0; i<stack_size; i++) {
        s[i].T1 = NULL;
        s[i].T1i = NULL;
        s[i].T2i = NULL;
        s[i].plan_fwd_1 = fft_wrapper_safe_plan_init();
        s[i].plan_inv_1 = fft_wrapper_safe_plan_init();
        s[i].plan_fwd_2 = fft_wrapper_safe_plan_init();
        s[i].plan_inv_2 = fft_wrapper_safe_plan_init();
    }

    ret_code = nse_finvscatter_recurse(s, eps_t, kappa, discretization, buf0,
        buf1, buf2);
    if (ret_code != SUCCESS)
        ret_code = E_SUBROUTINE(ret_code);

    for (i=0; i<stack_size; i++) {
        free(s[i].T1);
        free(s[i].T1i);
        free(s[i].T2i);
        fft_wrapper_destroy_plan(&s[i].plan_fwd_1);
        fft_wrapper_destroy_plan(&s[i].plan_inv_1);
        fft_wrapper_destroy_plan(&s[i].plan_fwd_2);
        fft_wrapper_destroy_plan(&s[i].plan_inv_2);
    }

leave_fun:
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    fft_wrapper_free(buf2);
    free(s);
    return ret_code;
}
