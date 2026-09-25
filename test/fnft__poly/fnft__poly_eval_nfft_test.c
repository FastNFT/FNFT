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

#include "fnft_config.h"
#include "fnft_errwarn.h"
#include "fnft__poly_eval_nfft.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_NFFT3

typedef FNFT_INT (*evaluator_t)(FNFT_UINT, FNFT_UINT,
        FNFT_COMPLEX const * const *, FNFT_UINT, FNFT_REAL const *,
        FNFT_COMPLEX * const *);

static long double complex reference_value(const FNFT_UINT degree,
        FNFT_COMPLEX const * const coefficients, const FNFT_REAL angle,
        const int chebyshev)
{
    long double complex result = 0.0L;
    FNFT_UINT i;

    /* Direct sums in extended precision, independent of NFFT mode packing
     * and of the double-precision Horner/Clenshaw fallback. */
    for (i = 0; i <= degree; i++) {
        const FNFT_UINT k = chebyshev ? i : degree - i;
        const long double phase = (long double)k*(long double)angle;
        const long double complex c =
                (long double)FNFT_CREAL(coefficients[i])
                + (long double complex)I
                *(long double)FNFT_CIMAG(coefficients[i]);
        if (c != 0.0L)
            result += c*(chebyshev ? cosl(phase)
                    : cosl(phase) + (long double complex)I*sinl(phase));
    }
    return result;
}

static int check_value(const FNFT_COMPLEX computed,
        const long double complex reference, const long double bound,
        const FNFT_UINT degree, const FNFT_UINT node,
        const int chebyshev, const char * const quantity)
{
    const long double error = cabsl((long double complex)computed - reference);

    if (!isfinite(FNFT_CREAL(computed)) || !isfinite(FNFT_CIMAG(computed))
            || !isfinite(error) || !(error <= bound)) {
        fprintf(stderr,
                "NFFT3 %s failure: chebyshev=%d degree=%lu node=%lu error=%.17g bound=%.17g\n",
                quantity, chebyshev, (unsigned long)degree,
                (unsigned long)node, (double)error, (double)bound);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int test_degree(const FNFT_UINT degree, const int chebyshev)
{
    enum { M = 23, NPOLY = 4 };
    FNFT_COMPLEX *coefficients[NPOLY] = {NULL, NULL, NULL, NULL};
    FNFT_COMPLEX const *polynomials[NPOLY];
    FNFT_COMPLEX output[NPOLY][M];
    FNFT_COMPLEX *values[NPOLY] = {output[0], output[1], output[2], output[3]};
    FNFT_REAL angles[M];
    long double norms[NPOLY] = {0.0L, 0.0L, 0.0L, 0.0L};
    const long double tolerance = FNFT_EPSILON
            *(64.0L + 16.0L*(long double)degree);
    const evaluator_t evaluate = chebyshev ? fnft__poly_eval_nfft_chebyshev
            : fnft__poly_eval_nfft_power;
    FNFT_UINT i, k, j;
    int result = EXIT_FAILURE;

    for (j = 0; j < NPOLY; j++) {
        coefficients[j] = calloc(degree + 1, sizeof(FNFT_COMPLEX));
        if (coefficients[j] == NULL)
            goto leave_fun;
        polynomials[j] = coefficients[j];
    }
    for (k = 0; k <= degree; k++) {
        const FNFT_REAL t = (FNFT_REAL)k;
        const FNFT_UINT index = chebyshev ? k : degree - k;
        coefficients[0][index] = (0.35*cos(0.31*t)
                + (FNFT_COMPLEX)I*0.21*sin(0.17*t))
                /((1.0 + t)*(1.0 + t));
    }
    /* Dense complex series, isolated highest mode, and a denominator whose
     * absolute value is bounded away from zero on the entire circle. */
    coefficients[1][chebyshev ? degree : 0] = 0.75
            - 0.4*(FNFT_COMPLEX)I;
    coefficients[2][chebyshev ? 0 : degree] = 2.0
            + 0.2*(FNFT_COMPLEX)I;
    coefficients[2][chebyshev ? degree : 0] += 0.125
            - 0.05*(FNFT_COMPLEX)I;
    for (j = 0; j < NPOLY; j++) {
        for (k = 0; k <= degree; k++)
            norms[j] += (long double)FNFT_CABS(coefficients[j][k]);
    }
    for (i = 0; i < M; i++) {
        const FNFT_REAL u = ((FNFT_REAL)i + 0.37)/((FNFT_REAL)M + 0.21);
        angles[i] = -0.97*FNFT_PI + 1.94*FNFT_PI*pow(u, 1.31);
    }
    angles[0] = 0.0;
    angles[1] = FNFT_PI;
    angles[2] = -FNFT_PI;
    angles[3] = nextafter(FNFT_PI, 0.0);
    angles[4] = nextafter(-FNFT_PI, 0.0);
    angles[5] = 1e-12;
    angles[6] = -1e-12;
    angles[7] = FNFT_PI/2.0;
    angles[8] = -FNFT_PI/2.0;

    if (evaluate(degree, NPOLY, polynomials, M, angles, values) != FNFT_SUCCESS)
        goto leave_fun;
    for (i = 0; i < M; i++) {
        long double complex reference[NPOLY];
        long double ratio_bound;
        for (j = 0; j < NPOLY; j++) {
            reference[j] = reference_value(degree, coefficients[j], angles[i],
                    chebyshev);
            /* Resolving angles in double precision introduces O(k*epsilon)
             * phase error in mode k, separately from NFFT window error. */
            if (check_value(output[j][i], reference[j], tolerance*norms[j],
                        degree, i, chebyshev, "polynomial") != EXIT_SUCCESS)
                goto leave_fun;
        }
        ratio_bound = 4.0L*tolerance
                *(norms[0] + cabsl(reference[0]/reference[2])*norms[2])
                /cabsl(reference[2]);
        if (check_value(output[0][i]/output[2][i], reference[0]/reference[2],
                    ratio_bound, degree, i, chebyshev, "ratio") != EXIT_SUCCESS)
            goto leave_fun;
    }
    result = EXIT_SUCCESS;

leave_fun:
    for (j = 0; j < NPOLY; j++)
        free(coefficients[j]);
    return result;
}

static int test_invalid_inputs(const evaluator_t evaluate)
{
    FNFT_COMPLEX coefficient = 1.0, output;
    FNFT_COMPLEX const *polynomials[1] = {&coefficient};
    FNFT_COMPLEX *values[1] = {&output};
    FNFT_REAL angle = 0.0;
    const FNFT_REAL invalid[] = {NAN, INFINITY, -INFINITY};
    FNFT_UINT i;

    if (evaluate(0, 0, polynomials, 1, &angle, values) == FNFT_SUCCESS
            || evaluate(0, 1, NULL, 1, &angle, values) == FNFT_SUCCESS
            || evaluate(0, 1, polynomials, 1, &angle, NULL) == FNFT_SUCCESS
            || evaluate(0, 1, polynomials, 0, &angle, values) == FNFT_SUCCESS
            || evaluate(0, 1, polynomials, 1, NULL, values) == FNFT_SUCCESS
            || evaluate((FNFT_UINT)-1, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS
            || evaluate((FNFT_UINT)INT_MAX/2, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS
            || evaluate(0, 1, polynomials, (FNFT_UINT)INT_MAX, &angle, values) == FNFT_SUCCESS)
        return EXIT_FAILURE;
    polynomials[0] = NULL;
    if (evaluate(0, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS)
        return EXIT_FAILURE;
    polynomials[0] = &coefficient;
    values[0] = NULL;
    if (evaluate(0, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS)
        return EXIT_FAILURE;
    values[0] = &output;
    for (i = 0; i < sizeof(invalid)/sizeof(invalid[0]); i++) {
        angle = invalid[i];
        if (evaluate(0, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS)
            return EXIT_FAILURE;
        angle = 0.0;
        coefficient = invalid[i];
        if (evaluate(0, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS)
            return EXIT_FAILURE;
        coefficient = (FNFT_COMPLEX)I*invalid[i];
        if (evaluate(0, 1, polynomials, 1, &angle, values) == FNFT_SUCCESS)
            return EXIT_FAILURE;
        coefficient = 1.0;
    }
    return EXIT_SUCCESS;
}

int main(void)
{
    const FNFT_UINT degrees[] = {0, 1, 2, 3, 7, 31, 32, 33, 64, 257, 4096, 65536};
    FNFT_UINT i;
    int chebyshev;

    for (chebyshev = 0; chebyshev <= 1; chebyshev++) {
        for (i = 0; i < sizeof(degrees)/sizeof(degrees[0]); i++) {
            if (test_degree(degrees[i], chebyshev) != EXIT_SUCCESS)
                return EXIT_FAILURE;
        }
    }
    if (test_invalid_inputs(fnft__poly_eval_nfft_power) != EXIT_SUCCESS
            || test_invalid_inputs(fnft__poly_eval_nfft_chebyshev) != EXIT_SUCCESS) {
        fprintf(stderr, "NFFT3 invalid-input checks failed.\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

#else

int main(void)
{
    return 77;
}

#endif
