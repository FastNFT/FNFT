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

/*
 * Reproduces the continuous-spectrum accuracy and timing sweep for the
 * public slow ES8 and fast Padé schemes used in the high-order comparison.
 * It calls the public fnft_nsev interface only; no code from NFT_fnft is used.
 *
 * Signal: q(t) = A*sech(t)^(1+i*C), A=5.2, C=4, t in [-30,30).
 * Spectral grid: M=2*N points on [-20,20] by default. A fixed M can be
 * selected with --spectrum-size M.
 * Time grids: 2^9,...,2^14 points.
 * Both signs kappa=-1,+1 are evaluated.
 * Each timing is repeated five times with a wall clock. Raw times, their
 * minimum/median/maximum and D*log2(D)^2 normalizations are written out.
 *
 * With no arguments, output:
 *   fes8_paper_continuous_spectrum.dat
 *   fes8_paper_invariants.dat  (for N=4096)
 * A single independent timing task can instead be selected with
 *   --task sigma log2N method [--spectrum-size M]
 * where method is an exact method label below or its numeric id. In task mode,
 * the header and one data row are written to stdout and no files are opened.
 */

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "fnft_config.h"
#include "fnft_nsev.h"

#define SIGNAL_AMPLITUDE 5.2
#define SIGNAL_CHIRP 4.0
#define TIME_LEFT (-30.0)
#define TIME_RIGHT 30.0
#define XI_LEFT (-20.0)
#define XI_RIGHT 20.0
#define SPECTRUM_MULTIPLIER ((FNFT_UINT)2)
#define INVARIANT_TIME_SIZE ((FNFT_UINT)4096)
#define TIMING_REPEATS 5

typedef struct {
    FNFT_UINT id;
    fnft_nse_discretization_t discretization;
    fnft_nsev_pade_representation_t representation;
    FNFT_UINT pade_degree;
    const char *name;
    const char *representation_name;
} method_t;

static const method_t methods[] = {
    {0, fnft_nse_discretization_ES8,
        fnft_nsev_pade_representation_DIRECT_CAYLEY,
        5, "ES8_PADE5", "slow_direct"},
    {1, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_DIRECT_CAYLEY,
        3, "FES8_PADE3", "direct"},
    {2, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_CHEBYSHEV_JOUKOWSKI,
        4, "FES8_PADE4_CHEB", "chebyshev_joukowski"},
    {3, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_DIRECT_CAYLEY,
        4, "FES8_PADE4", "direct"},
    {4, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_CHEBYSHEV_JOUKOWSKI,
        5, "FES8_PADE5_CHEB", "chebyshev_joukowski"},
    {5, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_DIRECT_CAYLEY,
        5, "FES8_PADE5", "direct"},
    {6, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_CHEBYSHEV_JOUKOWSKI,
        6, "FES8_PADE6_CHEB", "chebyshev_joukowski"},
    {7, fnft_nse_discretization_FES8_PADE,
        fnft_nsev_pade_representation_DIRECT_CAYLEY,
        6, "FES8_PADE6", "direct"}
};

static FNFT_REAL wall_time_seconds(void)
{
#if defined(_WIN32)
    LARGE_INTEGER counter, frequency;

    if (!QueryPerformanceFrequency(&frequency)
            || !QueryPerformanceCounter(&counter))
        return FNFT_NAN;
    return (FNFT_REAL)counter.QuadPart/(FNFT_REAL)frequency.QuadPart;
#else
    struct timeval value;

    if (gettimeofday(&value, NULL) != 0)
        return FNFT_NAN;
    return (FNFT_REAL)value.tv_sec + 1e-6*(FNFT_REAL)value.tv_usec;
#endif
}

static int is_finite_complex(const FNFT_COMPLEX z)
{
    return isfinite(FNFT_CREAL(z)) && isfinite(FNFT_CIMAG(z));
}

static int is_gamma_pole(const FNFT_COMPLEX z)
{
    return FNFT_FABS(FNFT_CIMAG(z)) <= 64.0*DBL_EPSILON
        && FNFT_CREAL(z) <= 0.0
        && FNFT_FABS(FNFT_CREAL(z) - FNFT_ROUND(FNFT_CREAL(z)))
            <= 64.0*DBL_EPSILON;
}

static FNFT_COMPLEX log_gamma_right_half_plane(FNFT_COMPLEX z)
{
    static const FNFT_REAL coefficients[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };
    FNFT_COMPLEX x, t;
    FNFT_UINT i;

    z -= 1.0;
    x = coefficients[0];
    for (i = 1; i < sizeof(coefficients)/sizeof(coefficients[0]); i++)
        x += coefficients[i]/(z + (FNFT_REAL)i);
    t = z + 7.5;
    return 0.91893853320467274178 + (z + 0.5)*FNFT_CLOG(t) - t
        + FNFT_CLOG(x);
}

static FNFT_COMPLEX log_gamma_lanczos(FNFT_COMPLEX z)
{
    FNFT_INT shift = 0, k;
    FNFT_COMPLEX log_product = 0.0;

    if (is_gamma_pole(z))
        return DBL_MAX;
    if (FNFT_CREAL(z) < 0.5)
        shift = (FNFT_INT)FNFT_CEIL(0.5 - FNFT_CREAL(z));
    for (k = 0; k < shift; k++) {
        const FNFT_COMPLEX zk = z + (FNFT_REAL)k;

        if (is_gamma_pole(zk))
            return DBL_MAX;
        log_product += FNFT_CLOG(zk);
    }
    return log_gamma_right_half_plane(z + (FNFT_REAL)shift) - log_product;
}

static FNFT_COMPLEX exp_log_gamma_ratio(const FNFT_COMPLEX numerator,
        const FNFT_COMPLEX denominator)
{
    const int numerator_pole = is_gamma_pole(numerator);
    const int denominator_pole = is_gamma_pole(denominator);

    if (denominator_pole) {
        if (!numerator_pole)
            return 0.0;
        {
            const FNFT_INT n = (FNFT_INT)FNFT_ROUND(-FNFT_CREAL(numerator));
            const FNFT_INT m = (FNFT_INT)FNFT_ROUND(-FNFT_CREAL(denominator));
            FNFT_REAL ratio = 1.0;
            FNFT_INT k;

            if (n > m) {
                for (k = m + 1; k <= n; k++)
                    ratio /= -(FNFT_REAL)k;
            } else {
                for (k = n + 1; k <= m; k++)
                    ratio *= -(FNFT_REAL)k;
            }
            return ratio;
        }
    }
    if (numerator_pole)
        return DBL_MAX;
    {
        const FNFT_COMPLEX value = FNFT_CEXP(log_gamma_lanczos(numerator)
                - log_gamma_lanczos(denominator));
        return is_finite_complex(value) ? value : DBL_MAX;
    }
}

static void chirped_secant_parameters(const FNFT_REAL sigma,
        FNFT_COMPLEX * const alpha, FNFT_COMPLEX * const beta)
{
    const FNFT_COMPLEX root = FNFT_CSQRT(sigma*SIGNAL_AMPLITUDE
            *SIGNAL_AMPLITUDE - 0.25*SIGNAL_CHIRP*SIGNAL_CHIRP);

    *alpha = root - 0.5*(FNFT_COMPLEX)I*SIGNAL_CHIRP;
    *beta = -root - 0.5*(FNFT_COMPLEX)I*SIGNAL_CHIRP;
}

static FNFT_COMPLEX exact_a(const FNFT_REAL xi, const FNFT_REAL sigma)
{
    FNFT_COMPLEX alpha, beta;
    const FNFT_COMPLEX gamma = 0.5
            - (FNFT_COMPLEX)I*(xi + 0.5*SIGNAL_CHIRP);

    chirped_secant_parameters(sigma, &alpha, &beta);
    return exp_log_gamma_ratio(gamma, gamma - beta)
        * exp_log_gamma_ratio(gamma - alpha - beta, gamma - alpha);
}

static FNFT_COMPLEX exact_b(const FNFT_REAL xi, const FNFT_REAL sigma)
{
    FNFT_COMPLEX alpha, beta;
    const FNFT_COMPLEX gamma = 0.5
            - (FNFT_COMPLEX)I*(xi + 0.5*SIGNAL_CHIRP);

    chirped_secant_parameters(sigma, &alpha, &beta);
    return FNFT_CPOW(2.0, -(FNFT_COMPLEX)I*SIGNAL_CHIRP)/SIGNAL_AMPLITUDE
        * exp_log_gamma_ratio(gamma + 2.0*(FNFT_COMPLEX)I*xi, beta)
        * exp_log_gamma_ratio(gamma, alpha);
}

static FNFT_REAL msre(const FNFT_COMPLEX * const exact,
        const FNFT_COMPLEX * const numerical, const FNFT_UINT M)
{
    FNFT_REAL error = 0.0;
    FNFT_UINT i;

    if (M == 0)
        return FNFT_NAN;
    for (i = 0; i < M; i++) {
        const FNFT_REAL exact_abs = FNFT_CABS(exact[i]);
        const FNFT_COMPLEX difference = numerical[i] - exact[i];
        const FNFT_REAL difference_abs = FNFT_CABS(difference);
        const FNFT_REAL denominator = exact_abs > 1.0
            ? exact_abs*exact_abs : 1.0;
        const FNFT_REAL term = difference_abs*difference_abs/denominator;

        if (!is_finite_complex(exact[i])
                || !is_finite_complex(numerical[i])
                || !is_finite_complex(difference) || !isfinite(exact_abs)
                || !isfinite(difference_abs) || !isfinite(denominator)
                || !isfinite(term))
            return FNFT_NAN;
        error += term;
        if (!isfinite(error))
            return FNFT_NAN;
    }
    error /= (FNFT_REAL)M;
    return isfinite(error) ? error : FNFT_NAN;
}

static FNFT_REAL abs_msre(const FNFT_COMPLEX * const exact,
        const FNFT_COMPLEX * const numerical, const FNFT_UINT M)
{
    FNFT_REAL error = 0.0;
    FNFT_UINT i;

    if (M == 0)
        return FNFT_NAN;
    for (i = 0; i < M; i++) {
        const FNFT_REAL exact_abs = FNFT_CABS(exact[i]);
        const FNFT_REAL numerical_abs = FNFT_CABS(numerical[i]);
        const FNFT_REAL difference = numerical_abs - exact_abs;
        const FNFT_REAL denominator = exact_abs > 1.0
            ? exact_abs*exact_abs : 1.0;
        const FNFT_REAL term = difference*difference/denominator;

        if (!is_finite_complex(exact[i])
                || !is_finite_complex(numerical[i]) || !isfinite(exact_abs)
                || !isfinite(numerical_abs) || !isfinite(difference)
                || !isfinite(denominator) || !isfinite(term))
            return FNFT_NAN;
        error += term;
        if (!isfinite(error))
            return FNFT_NAN;
    }
    error /= (FNFT_REAL)M;
    return isfinite(error) ? error : FNFT_NAN;
}

static void timing_summary(FNFT_REAL const times[TIMING_REPEATS],
        FNFT_REAL * const minimum, FNFT_REAL * const median,
        FNFT_REAL * const maximum)
{
    FNFT_REAL sorted[TIMING_REPEATS];
    FNFT_UINT i, j;

    for (i = 0; i < TIMING_REPEATS; i++) {
        sorted[i] = times[i];
        for (j = i; j > 0 && sorted[j] < sorted[j - 1]; j--) {
            const FNFT_REAL temporary = sorted[j];
            sorted[j] = sorted[j - 1];
            sorted[j - 1] = temporary;
        }
    }
    *minimum = sorted[0];
    *median = sorted[TIMING_REPEATS/2];
    *maximum = sorted[TIMING_REPEATS - 1];
}

static int build_exact_spectrum(const FNFT_INT kappa,
        const FNFT_UINT M,
        FNFT_COMPLEX * const a_exact, FNFT_COMPLEX * const b_exact,
        FNFT_COMPLEX * const r_exact)
{
    FNFT_UINT i;

    for (i = 0; i < M; i++) {
        const FNFT_REAL xi = XI_LEFT + (XI_RIGHT - XI_LEFT)*(FNFT_REAL)i
            /((FNFT_REAL)M - 1.0);

        a_exact[i] = exact_a(xi, (FNFT_REAL)kappa);
        b_exact[i] = exact_b(xi, (FNFT_REAL)kappa);
        if (!is_finite_complex(a_exact[i])
                || !is_finite_complex(b_exact[i])
                || !isfinite(FNFT_CABS(a_exact[i]))
                || FNFT_CABS(a_exact[i]) == 0.0)
            return FNFT_EC_OTHER;
        r_exact[i] = b_exact[i]/a_exact[i];
        if (!is_finite_complex(r_exact[i]))
            return FNFT_EC_OTHER;
    }
    return FNFT_SUCCESS;
}

static int flush_output(FILE * const stream)
{
    if (fflush(stream) == EOF || ferror(stream))
        return FNFT_EC_OTHER;
    return FNFT_SUCCESS;
}

static int run_method(FILE * const summary, FILE * const invariant_file,
        const method_t * const method, const FNFT_INT kappa,
        const FNFT_UINT D, const FNFT_UINT M,
        FNFT_COMPLEX const * const q,
        FNFT_COMPLEX const * const a_exact, FNFT_COMPLEX const * const b_exact,
        FNFT_COMPLEX const * const r_exact)
{
    const FNFT_REAL D_real = (FNFT_REAL)D;
    const FNFT_REAL M_real = (FNFT_REAL)M;
    const FNFT_REAL log2_D = FNFT_LOG2(D_real);
    const FNFT_REAL dt = (TIME_RIGHT - TIME_LEFT)/D_real;
    FNFT_REAL T[2] = {TIME_LEFT, TIME_RIGHT - dt};
    FNFT_REAL XI[2] = {XI_LEFT, XI_RIGHT};
    FNFT_COMPLEX *contspec = malloc(3*M*sizeof(FNFT_COMPLEX));
    FNFT_COMPLEX * const rho = contspec;
    FNFT_COMPLEX * const a = contspec == NULL ? NULL : contspec + M;
    FNFT_COMPLEX * const b = contspec == NULL ? NULL : contspec + 2*M;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    FNFT_REAL times[TIMING_REPEATS];
    FNFT_REAL time_minimum, time_median, time_maximum;
    FNFT_REAL time_median_per_D, time_median_per_D_log2_D2;
    FNFT_REAL max_invariant_error = 0.0, l2_invariant_error = 0.0;
    FNFT_REAL a_rmsre, a_abs_rmsre, b_rmsre, b_abs_rmsre, r_rmsre;
    FNFT_UINT launch, i;
    int ret_code = FNFT_SUCCESS;

    if (contspec == NULL)
        return FNFT_EC_NOMEM;
    opts.discretization = method->discretization;
    opts.pade_degree = method->pade_degree;
    opts.pade_representation = method->representation;
    opts.pade_h = 0.0;
    opts.contspec_type = fnft_nsev_cstype_BOTH;
    opts.normalization_flag = 1;

    for (launch = 0; launch < TIMING_REPEATS; launch++) {
        const FNFT_REAL begin = wall_time_seconds();
        FNFT_REAL end;

        if (!isfinite(begin)) {
            ret_code = FNFT_EC_OTHER;
            goto release_mem;
        }
        ret_code = fnft_nsev(D, q, T, M, contspec, XI,
                NULL, NULL, NULL, kappa, &opts);
        if (ret_code != FNFT_SUCCESS)
            goto release_mem;
        end = wall_time_seconds();
        if (!isfinite(end) || end < begin) {
            ret_code = FNFT_EC_OTHER;
            goto release_mem;
        }
        times[launch] = end - begin;
    }
    timing_summary(times, &time_minimum, &time_median, &time_maximum);

    for (i = 0; i < M; i++) {
        const FNFT_REAL invariant = FNFT_CABS(a[i])*FNFT_CABS(a[i])
            + (FNFT_REAL)kappa*FNFT_CABS(b[i])*FNFT_CABS(b[i]) - 1.0;
        const FNFT_REAL invariant_error = FNFT_FABS(invariant);

        if (!is_finite_complex(rho[i]) || !is_finite_complex(a[i])
                || !is_finite_complex(b[i]) || !isfinite(invariant)
                || !isfinite(invariant_error)) {
            ret_code = FNFT_EC_OTHER;
            goto release_mem;
        }
        if (invariant_error > max_invariant_error)
            max_invariant_error = invariant_error;
        l2_invariant_error += invariant*invariant;
        if (invariant_file != NULL && D == INVARIANT_TIME_SIZE) {
            const FNFT_REAL xi = XI_LEFT + (XI_RIGHT - XI_LEFT)
                *(FNFT_REAL)i/(M_real - 1.0);
            if (fprintf(invariant_file,
                "%d\t%lu\t%s\t%s\t%lu\t%lu\t%lu\t%.17g\t%.17g\t%.17g\t%.17g\n",
                (int)kappa, (unsigned long)method->id, method->name,
                method->representation_name,
                (unsigned long)method->pade_degree, (unsigned long)D,
                (unsigned long)M, (double)xi, (double)invariant_error,
                (double)FNFT_CABS(a[i]), (double)FNFT_CABS(b[i])) < 0
                    || flush_output(invariant_file) != FNFT_SUCCESS) {
                ret_code = FNFT_EC_OTHER;
                goto release_mem;
            }
        }
    }
    if (invariant_file != NULL && flush_output(invariant_file)
            != FNFT_SUCCESS) {
        ret_code = FNFT_EC_OTHER;
        goto release_mem;
    }
    l2_invariant_error = FNFT_SQRT(l2_invariant_error/M_real);
    a_rmsre = FNFT_SQRT(msre(a_exact, a, M));
    a_abs_rmsre = FNFT_SQRT(abs_msre(a_exact, a, M));
    b_rmsre = FNFT_SQRT(msre(b_exact, b, M));
    b_abs_rmsre = FNFT_SQRT(abs_msre(b_exact, b, M));
    r_rmsre = FNFT_SQRT(msre(r_exact, rho, M));
    time_median_per_D = time_median/D_real;
    time_median_per_D_log2_D2 = time_median
        /(D_real*log2_D*log2_D);
    if (!isfinite(time_minimum) || !isfinite(time_median)
            || !isfinite(time_maximum) || !isfinite(max_invariant_error)
            || !isfinite(l2_invariant_error) || !isfinite(a_rmsre)
            || !isfinite(a_abs_rmsre) || !isfinite(b_rmsre)
            || !isfinite(b_abs_rmsre) || !isfinite(r_rmsre)
            || !isfinite(time_median_per_D)
            || !isfinite(time_median_per_D_log2_D2)) {
        ret_code = FNFT_EC_OTHER;
        goto release_mem;
    }

    if (fprintf(summary,
        "%d\t%lu\t%s\t%s\t%lu\t%.17g\t%.17g\t%lu\t%u\t%lu\t%.17g\t%u"
        "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
        "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g"
        "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
        (int)kappa, (unsigned long)method->id, method->name,
        method->representation_name,
        (unsigned long)method->pade_degree,
        (double)SIGNAL_AMPLITUDE, (double)SIGNAL_CHIRP,
        (unsigned long)D, (unsigned int)FNFT_ROUND(log2_D),
        (unsigned long)M, (double)(M_real/D_real),
        (unsigned int)TIMING_REPEATS,
        (double)times[0], (double)times[1], (double)times[2],
        (double)times[3], (double)times[4],
        (double)time_minimum, (double)time_median, (double)time_maximum,
        (double)time_median_per_D, (double)time_median_per_D_log2_D2,
        (double)a_rmsre, (double)a_abs_rmsre,
        (double)b_rmsre, (double)b_abs_rmsre, (double)r_rmsre,
        (double)max_invariant_error, (double)l2_invariant_error) < 0
            || flush_output(summary) != FNFT_SUCCESS) {
        ret_code = FNFT_EC_OTHER;
        goto release_mem;
    }

release_mem:
    free(contspec);
    return ret_code;
}

static int write_summary_header(FILE * const stream)
{
    if (fprintf(stream,
        "sigma\tmethod_id\tmethod\trepresentation\tpade_degree\tamplitude"
        "\tchirp\tN\tlog2N\tM\tM_over_N\ttiming_repeats"
        "\ttime_0_sec\ttime_1_sec\ttime_2_sec\ttime_3_sec\ttime_4_sec"
        "\ttime_min_sec\ttime_median_sec\ttime_max_sec"
        "\ttime_median_per_N\ttime_median_per_N_log2N2"
        "\ta_rmsre\ta_abs_rmsre"
        "\tb_rmsre\tb_abs_rmsre\tr_rmsre\tmax_invariant_error"
        "\tl2_invariant_error\n") < 0)
        return FNFT_EC_OTHER;
    return flush_output(stream);
}

static int write_usage(FILE * const stream, const char * const program)
{
    FNFT_UINT i;

    if (fprintf(stream, "Usage: %s\n", program) < 0
            || fprintf(stream, "       %s --spectrum-size M\n",
                program) < 0
            || fprintf(stream,
                "       %s --task {-1|+1} {9..14} {method-label|method-id}\n",
                program) < 0
            || fprintf(stream,
                "       %s --task {-1|+1} {9..14} {method-label|method-id}"
                " --spectrum-size M\n", program) < 0
            || fprintf(stream, "Methods:\n") < 0)
        return FNFT_EC_OTHER;
    for (i = 0; i < sizeof(methods)/sizeof(methods[0]); i++) {
        if (fprintf(stream, "  %lu  %s\n", (unsigned long)methods[i].id,
                methods[i].name) < 0)
            return FNFT_EC_OTHER;
    }
    return flush_output(stream);
}

static int parse_integer(const char * const text, long * const value)
{
    char *end;

    errno = 0;
    *value = strtol(text, &end, 10);
    return errno == 0 && end != text && *end == '\0';
}

static int spectrum_size_is_valid(const FNFT_UINT M)
{
    const FNFT_UINT maximum = (FNFT_UINT)-1;

    /* Match the NFFT evaluator's oversampling guard for cutoff 12. */
    return M >= 2 && M <= (FNFT_UINT)INT_MAX/26U
        && M <= maximum/(26U*sizeof(FNFT_REAL))
        && M <= maximum/(3U*sizeof(FNFT_COMPLEX));
}

static int parse_spectrum_size(const char * const text, FNFT_UINT * const M)
{
    char *end;
    unsigned long long parsed;

    if (text == NULL || text[0] < '0' || text[0] > '9')
        return 0;
    errno = 0;
    parsed = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0'
            || parsed > (unsigned long long)((FNFT_UINT)-1))
        return 0;
    *M = (FNFT_UINT)parsed;
    return spectrum_size_is_valid(*M);
}

static const method_t *find_method(const char * const selector)
{
    FNFT_UINT i;
    long id;

    for (i = 0; i < sizeof(methods)/sizeof(methods[0]); i++) {
        if (strcmp(selector, methods[i].name) == 0)
            return &methods[i];
    }
    if (!parse_integer(selector, &id) || id < 0)
        return NULL;
    for (i = 0; i < sizeof(methods)/sizeof(methods[0]); i++) {
        if ((unsigned long)id == (unsigned long)methods[i].id)
            return &methods[i];
    }
    return NULL;
}

static int run_case(FILE * const summary, FILE * const invariant_file,
        const method_t * const method, const FNFT_INT kappa,
        const FNFT_UINT log2_D, const FNFT_UINT M)
{
    const FNFT_UINT D = (FNFT_UINT)1U << log2_D;
    const FNFT_REAL dt = (TIME_RIGHT - TIME_LEFT)/(FNFT_REAL)D;
    FNFT_COMPLEX *q = NULL;
    FNFT_COMPLEX *a_exact = NULL;
    FNFT_COMPLEX *b_exact = NULL;
    FNFT_COMPLEX *r_exact = NULL;
    FNFT_UINT i;
    int ret_code;

    if (!spectrum_size_is_valid(M)) {
        ret_code = FNFT_EC_INVALID_ARGUMENT;
        goto release_mem;
    }
    q = malloc(D*sizeof(FNFT_COMPLEX));
    a_exact = malloc(M*sizeof(FNFT_COMPLEX));
    b_exact = malloc(M*sizeof(FNFT_COMPLEX));
    r_exact = malloc(M*sizeof(FNFT_COMPLEX));
    if (q == NULL || a_exact == NULL || b_exact == NULL || r_exact == NULL) {
        ret_code = FNFT_EC_NOMEM;
        goto release_mem;
    }
    ret_code = build_exact_spectrum(kappa, M, a_exact, b_exact, r_exact);
    if (ret_code != FNFT_SUCCESS)
        goto release_mem;
    for (i = 0; i < D; i++) {
        const FNFT_REAL t = TIME_LEFT + dt*(FNFT_REAL)i;
        const FNFT_REAL sech = 1.0/FNFT_COSH(t);

        q[i] = SIGNAL_AMPLITUDE
            * FNFT_CPOW(sech, 1.0 + (FNFT_COMPLEX)I*SIGNAL_CHIRP);
        if (!is_finite_complex(q[i])) {
            ret_code = FNFT_EC_OTHER;
            goto release_mem;
        }
    }
    ret_code = run_method(summary, invariant_file, method, kappa, D, M, q,
        a_exact, b_exact, r_exact);

release_mem:
    free(q);
    free(a_exact);
    free(b_exact);
    free(r_exact);
    return ret_code;
}

int main(int argc, char **argv)
{
#ifndef HAVE_NFFT3
    if (fputs("This timing benchmark requires FNFT built with -DENABLE_NFFT=ON.\n",
            stderr) == EOF)
        return EXIT_FAILURE;
    (void)flush_output(stderr);
    return EXIT_FAILURE;
#endif
    FILE *summary = NULL, *invariant_file = NULL;
    FNFT_UINT k, method_index;
    FNFT_UINT fixed_M = 0;
    FNFT_INT kappa;
    int ret_code = EXIT_SUCCESS;

    if (argc == 3 && strcmp(argv[1], "--spectrum-size") == 0) {
        if (!parse_spectrum_size(argv[2], &fixed_M)) {
            (void)write_usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }
    } else if (argc != 1) {
        const method_t *method;
        long parsed_kappa, parsed_log2_D;
        FNFT_UINT M;

        if ((argc != 5 && argc != 7) || strcmp(argv[1], "--task") != 0
                || !parse_integer(argv[2], &parsed_kappa)
                || (parsed_kappa != -1 && parsed_kappa != 1)
                || !parse_integer(argv[3], &parsed_log2_D)
                || parsed_log2_D < 9 || parsed_log2_D > 14
                || (method = find_method(argv[4])) == NULL
                || (argc == 7 && (strcmp(argv[5], "--spectrum-size") != 0
                    || !parse_spectrum_size(argv[6], &fixed_M)))) {
            (void)write_usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }
        M = fixed_M == 0
            ? SPECTRUM_MULTIPLIER*((FNFT_UINT)1U << parsed_log2_D)
            : fixed_M;
        ret_code = write_summary_header(stdout);
        if (ret_code != FNFT_SUCCESS)
            return EXIT_FAILURE;
        ret_code = run_case(stdout, NULL, method, (FNFT_INT)parsed_kappa,
            (FNFT_UINT)parsed_log2_D, M);
        if (ret_code == FNFT_SUCCESS
                && flush_output(stdout) != FNFT_SUCCESS)
            ret_code = FNFT_EC_OTHER;
        return ret_code == FNFT_SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    summary = fopen("fes8_paper_continuous_spectrum.dat", "w");
    invariant_file = fopen("fes8_paper_invariants.dat", "w");
    if (summary == NULL || invariant_file == NULL) {
        ret_code = EXIT_FAILURE;
        goto release_mem;
    }
    ret_code = write_summary_header(summary);
    if (ret_code != FNFT_SUCCESS)
        goto release_mem;
    if (fprintf(invariant_file,
        "sigma\tmethod_id\tmethod\trepresentation\tpade_degree\tN\tM\txi"
        "\tinvariant_error\tabs_a\tabs_b\n") < 0
            || flush_output(invariant_file) != FNFT_SUCCESS) {
        ret_code = FNFT_EC_OTHER;
        goto release_mem;
    }

    for (kappa = -1; kappa <= 1; kappa += 2) {
        for (k = 9; k <= 14; k++) {
            const FNFT_UINT D = (FNFT_UINT)1U << k;
            const FNFT_UINT M = fixed_M == 0
                ? SPECTRUM_MULTIPLIER*D : fixed_M;
            const FNFT_REAL dt = (TIME_RIGHT - TIME_LEFT)/(FNFT_REAL)D;
            FNFT_COMPLEX *q = NULL;
            FNFT_COMPLEX *a_exact = NULL;
            FNFT_COMPLEX *b_exact = NULL;
            FNFT_COMPLEX *r_exact = NULL;
            FNFT_UINT i;

            if (!spectrum_size_is_valid(M)) {
                free(q);
                free(a_exact);
                free(b_exact);
                free(r_exact);
                ret_code = FNFT_EC_INVALID_ARGUMENT;
                goto release_mem;
            }
            q = malloc(D*sizeof(FNFT_COMPLEX));
            a_exact = malloc(M*sizeof(FNFT_COMPLEX));
            b_exact = malloc(M*sizeof(FNFT_COMPLEX));
            r_exact = malloc(M*sizeof(FNFT_COMPLEX));
            if (q == NULL || a_exact == NULL || b_exact == NULL
                    || r_exact == NULL) {
                free(q);
                free(a_exact);
                free(b_exact);
                free(r_exact);
                ret_code = FNFT_EC_NOMEM;
                goto release_mem;
            }
            ret_code = build_exact_spectrum(kappa, M, a_exact, b_exact,
                r_exact);
            if (ret_code != FNFT_SUCCESS) {
                free(q);
                free(a_exact);
                free(b_exact);
                free(r_exact);
                goto release_mem;
            }
            for (i = 0; i < D; i++) {
                const FNFT_REAL t = TIME_LEFT + dt*(FNFT_REAL)i;
                const FNFT_REAL sech = 1.0/FNFT_COSH(t);

                q[i] = SIGNAL_AMPLITUDE
                    * FNFT_CPOW(sech,
                        1.0 + (FNFT_COMPLEX)I*SIGNAL_CHIRP);
                if (!is_finite_complex(q[i])) {
                    free(q);
                    free(a_exact);
                    free(b_exact);
                    free(r_exact);
                    ret_code = FNFT_EC_OTHER;
                    goto release_mem;
                }
            }

            for (method_index = 0;
                    method_index < sizeof(methods)/sizeof(methods[0]);
                    method_index++) {
                if (printf("sigma=%d N=%lu %s degree=%lu %s\n",
                    (int)kappa, (unsigned long)D,
                    methods[method_index].name,
                    (unsigned long)methods[method_index].pade_degree,
                    methods[method_index].representation_name) < 0
                        || flush_output(stdout) != FNFT_SUCCESS) {
                    free(q);
                    free(a_exact);
                    free(b_exact);
                    free(r_exact);
                    ret_code = FNFT_EC_OTHER;
                    goto release_mem;
                }
                ret_code = run_method(summary, invariant_file,
                    &methods[method_index], kappa, D, M, q,
                    a_exact, b_exact, r_exact);
                if (ret_code != FNFT_SUCCESS) {
                    free(q);
                    free(a_exact);
                    free(b_exact);
                    free(r_exact);
                    goto release_mem;
                }
            }
            free(q);
            free(a_exact);
            free(b_exact);
            free(r_exact);
        }
    }

release_mem:
    if (summary != NULL) {
        if (ferror(summary))
            ret_code = FNFT_EC_OTHER;
        if (fclose(summary) == EOF)
            ret_code = FNFT_EC_OTHER;
    }
    if (invariant_file != NULL) {
        if (ferror(invariant_file))
            ret_code = FNFT_EC_OTHER;
        if (fclose(invariant_file) == EOF)
            ret_code = FNFT_EC_OTHER;
    }
    return ret_code == FNFT_SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
