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
 * Igor Chekhovskoy (NSU, FRC ICT) 2026.
 */

#include "fnft__poly_eval_nfft.h"
#include "fnft_config.h"
#include "fnft__errwarn.h"

#include <limits.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_NFFT3
#include <nfft3.h>

enum { NFFT_CUTOFF = 12 };

static FNFT_REAL fnft__poly_eval_nfft_node(FNFT_REAL angle)
{
    const FNFT_REAL pi = FNFT_PI;
    FNFT_REAL x;

    /* The caller supplies principal angles; also accept other finite angles
     * without overflowing the phase multiplication or losing full periods. */
    if (angle < -pi || angle > pi)
        angle = atan2(sin(angle), cos(angle));
    x = -angle/(2.0*pi);
    if (x >= 0.5)
        x -= 1.0;
    return x;
}

static FNFT_INT fnft__poly_eval_nfft_prepare_plan(
        const FNFT_UINT N,
        const FNFT_UINT degree,
        const FNFT_UINT npoly,
        FNFT_COMPLEX const * const * const polynomials,
        const FNFT_UINT M,
        FNFT_REAL const * const angles,
        FNFT_COMPLEX * const * const values,
        nfft_plan * const plan)
{
    FNFT_UINT i, j;
    int bandwidth, fft_length = 32;
    const char *error_message;

    if (N == 0 || N > INT_MAX || (N & 1U) != 0)
        return FNFT__E_INVALID_ARGUMENT(N);
    if (M == 0 || M > INT_MAX/(2*NFFT_CUTOFF + 2)
            || M > (size_t)-1/((2*NFFT_CUTOFF + 2)*sizeof(FNFT_REAL)))
        return FNFT__E_INVALID_ARGUMENT(M);
    if (npoly == 0 || polynomials == NULL || values == NULL)
        return FNFT__E_INVALID_ARGUMENT(polynomials);
    if (angles == NULL || plan == NULL)
        return FNFT__E_INVALID_ARGUMENT(angles);

    /* Check the oversampled FFT size before NFFT3 can allocate any memory. */
    while ((FNFT_UINT)fft_length < N) {
        if (fft_length > INT_MAX/2)
            return FNFT__E_INVALID_ARGUMENT(N);
        fft_length *= 2;
    }
    if (fft_length > INT_MAX/2
            || (size_t)fft_length > (size_t)-1/(2*sizeof(FNFT_COMPLEX)))
        return FNFT__E_INVALID_ARGUMENT(N);
    fft_length *= 2;
    bandwidth = (int)N;

    /* nfft_check does not reject NaN nodes. Validate before initialization
     * and, in particular, before node-dependent window precomputation. */
    for (i = 0; i < M; i++) {
        if (!isfinite(angles[i]))
            return FNFT__E_INVALID_ARGUMENT(angles);
    }
    for (i = 0; i < npoly; i++) {
        if (polynomials[i] == NULL || values[i] == NULL)
            return FNFT__E_INVALID_ARGUMENT(polynomials);
        for (j = 0; j <= degree; j++) {
            if (!isfinite(FNFT_CREAL(polynomials[i][j]))
                    || !isfinite(FNFT_CIMAG(polynomials[i][j])))
                return FNFT__E_INVALID_ARGUMENT(polynomials);
        }
    }

    nfft_init_guru(plan, 1, &bandwidth, (int)M, &fft_length, NFFT_CUTOFF,
            PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F
            | FFTW_INIT | FFT_OUT_OF_PLACE,
            FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    if (plan->x == NULL || plan->f_hat == NULL || plan->f == NULL) {
        nfft_finalize(plan);
        return FNFT__E_NOMEM;
    }
    for (i = 0; i < M; i++)
        plan->x[i] = fnft__poly_eval_nfft_node(angles[i]);

    error_message = nfft_check(plan);
    if (error_message != NULL) {
        /* Form the FNFT error while the NFFT plan is still alive. */
        const FNFT_INT ret_code = FNFT__E_OTHER(error_message);
        nfft_finalize(plan);
        return ret_code;
    }
    if (plan->flags & PRE_ONE_PSI)
        nfft_precompute_one_psi(plan);

    return FNFT_SUCCESS;
}

FNFT_INT fnft__poly_eval_nfft_power(
        const FNFT_UINT degree,
        const FNFT_UINT npoly,
        FNFT_COMPLEX const * const * const polynomials,
        const FNFT_UINT M,
        FNFT_REAL const * const angles,
        FNFT_COMPLEX * const * const values)
{
    nfft_plan plan;
    FNFT_UINT N, half, polynomial_index, i;
    FNFT_INT ret_code;

    if (degree > (FNFT_UINT)INT_MAX - 2U)
        return FNFT__E_INVALID_ARGUMENT(degree);
    N = degree + 1;
    if ((N & 1U) != 0)
        N++;
    half = N/2;

    ret_code = fnft__poly_eval_nfft_prepare_plan(N, degree, npoly,
            polynomials, M, angles, values, &plan);
    if (ret_code != FNFT_SUCCESS)
        return ret_code;

    for (polynomial_index = 0; polynomial_index < npoly;
            polynomial_index++) {
        FNFT_COMPLEX const * const polynomial = polynomials[polynomial_index];
        FNFT_COMPLEX * const output = values[polynomial_index];

        memset(plan.f_hat, 0, (size_t)plan.N_total*sizeof(*plan.f_hat));
        for (i = 0; i <= degree; i++)
            plan.f_hat[degree - i] = polynomial[i];

        nfft_trafo(&plan);
        for (i = 0; i < M; i++) {
            /* NFFT modes run from -N/2 to N/2-1. Undo the common shift,
             * using the same reduced node as the transform itself. */
            const FNFT_REAL phase = -2.0*FNFT_PI
                    * remainder((FNFT_REAL)half*plan.x[i], 1.0);
            output[i] = plan.f[i]*FNFT_CEXP((FNFT_COMPLEX)I*phase);
            if (!isfinite(FNFT_CREAL(output[i]))
                    || !isfinite(FNFT_CIMAG(output[i]))) {
                ret_code = FNFT__E_OTHER("Non-finite NFFT3 polynomial value.");
                goto leave_fun;
            }
        }
    }

leave_fun:
    nfft_finalize(&plan);
    return ret_code;
}

FNFT_INT fnft__poly_eval_nfft_chebyshev(
        const FNFT_UINT degree,
        const FNFT_UINT npoly,
        FNFT_COMPLEX const * const * const polynomials,
        const FNFT_UINT M,
        FNFT_REAL const * const angles,
        FNFT_COMPLEX * const * const values)
{
    nfft_plan plan;
    FNFT_UINT N, half, polynomial_index, i, k;
    FNFT_INT ret_code;

    if (degree > (FNFT_UINT)INT_MAX/2 - 1U)
        return FNFT__E_INVALID_ARGUMENT(degree);
    N = 2*(degree + 1);
    half = N/2;

    ret_code = fnft__poly_eval_nfft_prepare_plan(N, degree, npoly,
            polynomials, M, angles, values, &plan);
    if (ret_code != FNFT_SUCCESS)
        return ret_code;

    for (polynomial_index = 0; polynomial_index < npoly;
            polynomial_index++) {
        FNFT_COMPLEX const * const polynomial = polynomials[polynomial_index];
        FNFT_COMPLEX * const output = values[polynomial_index];

        memset(plan.f_hat, 0, (size_t)plan.N_total*sizeof(*plan.f_hat));
        plan.f_hat[half] = polynomial[0];
        for (k = 1; k <= degree; k++) {
            const FNFT_COMPLEX coefficient = 0.5*polynomial[k];
            plan.f_hat[half + k] = coefficient;
            plan.f_hat[half - k] = coefficient;
        }

        nfft_trafo(&plan);
        for (i = 0; i < M; i++) {
            output[i] = plan.f[i];
            if (!isfinite(FNFT_CREAL(output[i]))
                    || !isfinite(FNFT_CIMAG(output[i]))) {
                ret_code = FNFT__E_OTHER("Non-finite NFFT3 polynomial value.");
                goto leave_fun;
            }
        }
    }

leave_fun:
    nfft_finalize(&plan);
    return ret_code;
}

#else

FNFT_INT fnft__poly_eval_nfft_power(
        const FNFT_UINT degree,
        const FNFT_UINT npoly,
        FNFT_COMPLEX const * const * const polynomials,
        const FNFT_UINT M,
        FNFT_REAL const * const angles,
        FNFT_COMPLEX * const * const values)
{
    (void)degree;
    (void)npoly;
    (void)polynomials;
    (void)M;
    (void)angles;
    (void)values;
    return FNFT__E_NOT_YET_IMPLEMENTED(NFFT3,
            "Rebuild FNFT with -DENABLE_NFFT=ON.");
}

FNFT_INT fnft__poly_eval_nfft_chebyshev(
        const FNFT_UINT degree,
        const FNFT_UINT npoly,
        FNFT_COMPLEX const * const * const polynomials,
        const FNFT_UINT M,
        FNFT_REAL const * const angles,
        FNFT_COMPLEX * const * const values)
{
    (void)degree;
    (void)npoly;
    (void)polynomials;
    (void)M;
    (void)angles;
    (void)values;
    return FNFT__E_NOT_YET_IMPLEMENTED(NFFT3,
            "Rebuild FNFT with -DENABLE_NFFT=ON.");
}

#endif
