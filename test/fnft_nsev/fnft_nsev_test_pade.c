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
 * Igor Chekhovskoy 2026.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__akns_discretization.h"
#include "fnft__nse_discretization.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    COMPLEX entry[4];
} matrix_t;

static matrix_t matrix_identity(void)
{
    const matrix_t result = {{1.0, 0.0, 0.0, 1.0}};
    return result;
}

static matrix_t matrix_multiply(matrix_t const * const left,
        matrix_t const * const right)
{
    matrix_t result;

    result.entry[0] = left->entry[0]*right->entry[0]
            + left->entry[1]*right->entry[2];
    result.entry[1] = left->entry[0]*right->entry[1]
            + left->entry[1]*right->entry[3];
    result.entry[2] = left->entry[2]*right->entry[0]
            + left->entry[3]*right->entry[2];
    result.entry[3] = left->entry[2]*right->entry[1]
            + left->entry[3]*right->entry[3];
    return result;
}

static matrix_t pointwise_es8_z(const UINT D, const UINT index,
        COMPLEX const * const q, COMPLEX const * const r,
        const REAL eps_t, const COMPLEX lambda)
{
    const UINT im1 = (index + D - 1) % D;
    const UINT ip1 = (index + 1) % D;
    const UINT im2 = (index + D - 2) % D;
    const UINT ip2 = (index + 2) % D;
    const UINT im3 = (index + D - 3) % D;
    const UINT ip3 = (index + 3) % D;
    const COMPLEX q_samples[7] = {
        q[im3], q[im2], q[im1], q[index], q[ip1], q[ip2], q[ip3]
    };
    const COMPLEX r_samples[7] = {
        r[im3], r[im2], r[im1], r[index], r[ip1], r[ip2], r[ip3]
    };
    fnft__akns_es8_stencil_t qs, rs;
    COMPLEX q_values[7], r_values[7], coefficients[24];
    const COMPLEX x = eps_t*lambda;
    matrix_t result = {{0.0, 0.0, 0.0, 0.0}};
    UINT i, j;

    fnft__akns_es8_stencil(q_samples, eps_t, &qs);
    fnft__akns_es8_stencil(r_samples, eps_t, &rs);
    q_values[0] = qs.value;
    q_values[1] = qs.first;
    q_values[2] = qs.second;
    q_values[3] = qs.third;
    q_values[4] = qs.fourth;
    q_values[5] = qs.fifth;
    q_values[6] = qs.sixth;
    r_values[0] = rs.value;
    r_values[1] = rs.first;
    r_values[2] = rs.second;
    r_values[3] = rs.third;
    r_values[4] = rs.fourth;
    r_values[5] = rs.fifth;
    r_values[6] = rs.sixth;
    fnft__akns_es8_z_coefficients(q_values, r_values, coefficients);
    for (i = 0; i < 4; i++) {
        for (j = 6; j-- > 0; )
            result.entry[i] = result.entry[i]*x + coefficients[4*j + i];
    }
    return result;
}

static matrix_t pointwise_pade(matrix_t const * const z,
        const UINT pade_degree)
{
    matrix_t p_plus = matrix_identity();
    matrix_t p_minus = matrix_identity();
    matrix_t power = matrix_identity();
    matrix_t inverse_minus, result;
    REAL coefficient = 1.0;
    COMPLEX determinant;
    UINT i, j;

    for (j = 1; j <= pade_degree; j++) {
        power = matrix_multiply(&power, z);
        coefficient *= (REAL)(pade_degree + 1 - j)
                / ((REAL)j*(REAL)(2*pade_degree + 1 - j));
        for (i = 0; i < 4; i++) {
            p_plus.entry[i] += coefficient*power.entry[i];
            p_minus.entry[i] += ((j & 1U) == 0 ? coefficient : -coefficient)
                    *power.entry[i];
        }
    }
    determinant = p_minus.entry[0]*p_minus.entry[3]
            - p_minus.entry[1]*p_minus.entry[2];
    inverse_minus.entry[0] = p_minus.entry[3]/determinant;
    inverse_minus.entry[1] = -p_minus.entry[1]/determinant;
    inverse_minus.entry[2] = -p_minus.entry[2]/determinant;
    inverse_minus.entry[3] = p_minus.entry[0]/determinant;
    result = matrix_multiply(&inverse_minus, &p_plus);
    return result;
}

static COMPLEX pointwise_reflection(const UINT D,
        COMPLEX const * const q, const REAL T[2], const COMPLEX lambda,
        const INT kappa, const UINT pade_degree)
{
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    COMPLEX *r = malloc(D*sizeof(COMPLEX));
    matrix_t product = matrix_identity();
    REAL phase_factor;
    UINT i;

    if (r == NULL)
        return NAN;
    for (i = 0; i < D; i++)
        r[i] = -(REAL)kappa*CONJ(q[i]);
    for (i = 0; i < D; i++) {
        const matrix_t z = pointwise_es8_z(D, i, q, r, eps_t, lambda);
        const matrix_t step = pointwise_pade(&z, pade_degree);
        product = matrix_multiply(&step, &product);
    }
    free(r);
    if (nse_discretization_phase_factor_rho(eps_t, T[1], &phase_factor,
                nse_discretization_FES8_PADE) != SUCCESS)
        return NAN;
    return product.entry[2]*CEXP(I*lambda*phase_factor)/product.entry[0];
}

static INT test_zero_signal(const nse_discretization_t discretization,
        const UINT pade_degree, const INT kappa, const INT normalization_flag,
        const fnft_nsev_pade_representation_t pade_representation)
{
    const UINT D = 16, M = 5;
    const REAL T[2] = {-1.0, 1.0};
    const REAL XI[2] = {-2.0, 2.0};
    COMPLEX q[16] = {0.0};
    COMPLEX spectrum[15];
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    UINT i;
    INT ret_code;

    opts.discretization = discretization;
    opts.pade_degree = pade_degree;
    opts.pade_representation = pade_representation;
    opts.contspec_type = nsev_cstype_BOTH;
    opts.normalization_flag = normalization_flag;
    ret_code = fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL,
            kappa, &opts);
    if (ret_code != SUCCESS)
        return ret_code;
    for (i = 0; i < M; i++) {
        if (FNFT_CABS(spectrum[i]) > 2e-11
                || FNFT_FABS(FNFT_CABS(spectrum[M + i]) - 1.0) > 2e-11
                || FNFT_CABS(spectrum[2*M + i]) > 2e-11) {
            fprintf(stderr, "Padé zero-signal failure: degree=%lu kappa=%d norm=%d i=%lu rho=%.3e a=(%.6e,%.6e) b=%.3e\n",
                    (unsigned long)pade_degree, (int)kappa,
                    (int)normalization_flag, (unsigned long)i,
                    FNFT_CABS(spectrum[i]), FNFT_CREAL(spectrum[M + i]),
                    FNFT_CIMAG(spectrum[M + i]),
                    FNFT_CABS(spectrum[2*M + i]));
            return E_TEST_FAILED;
        }
    }
    return SUCCESS;
}

static REAL continuous_spectrum_error(const nsev_testcases_t testcase,
        const UINT D, const nse_discretization_t discretization,
        const UINT pade_degree,
        const fnft_nsev_pade_representation_t pade_representation,
        INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL, *computed = NULL;
    REAL T[2], XI[2], error = 0.0, norm = 0.0;
    UINT M, K, i;
    INT kappa;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();

    *ret_code = nsev_testcases(testcase, D, &q, T, &M, &exact, &ab, XI,
            &K, &bound_states, &normconsts, &residues, &kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    computed = malloc(M*sizeof(COMPLEX));
    if (computed == NULL) {
        *ret_code = E_NOMEM;
        goto release_mem;
    }
    opts.discretization = discretization;
    opts.pade_degree = pade_degree;
    opts.pade_representation = pade_representation;
    opts.contspec_type = nsev_cstype_REFLECTION_COEFFICIENT;
    *ret_code = fnft_nsev(D, q, T, M, computed, XI, NULL, NULL, NULL,
            kappa, &opts);
    if (*ret_code != SUCCESS)
        goto release_mem;
    for (i = 0; i < M; i++) {
        error += FNFT_CABS(computed[i] - exact[i]);
        norm += FNFT_CABS(exact[i]);
    }
    error /= 1.0 + norm;

release_mem:
    free(q);
    free(exact);
    free(ab);
    free(bound_states);
    free(normconsts);
    free(residues);
    free(computed);
    return error;
}

static REAL pointwise_spectrum_error(const nsev_testcases_t testcase,
        const UINT D, const UINT pade_degree, INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL;
    REAL T[2], XI[2], error = 0.0, norm = 0.0;
    UINT M, K, i;
    INT kappa;

    *ret_code = nsev_testcases(testcase, D, &q, T, &M, &exact, &ab, XI,
            &K, &bound_states, &normconsts, &residues, &kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    for (i = 0; i < M; i++) {
        const COMPLEX lambda = XI[0] + (XI[1] - XI[0])*i/(M - 1);
        const COMPLEX computed = pointwise_reflection(D, q, T, lambda,
                kappa, pade_degree);

        if (isnan(CREAL(computed)) || isnan(CIMAG(computed))) {
            *ret_code = E_TEST_FAILED;
            goto release_mem;
        }
        error += CABS(computed - exact[i]);
        norm += CABS(exact[i]);
    }
    error /= 1.0 + norm;

release_mem:
    free(q);
    free(exact);
    free(ab);
    free(bound_states);
    free(normconsts);
    free(residues);
    return error;
}

static INT test_public_es8_small_grid(void)
{
    const UINT D = 7, M = 2;
    const REAL T[2] = {-0.21, 0.21};
    const REAL XI[2] = {0.63, 0.64};
    const UINT degrees[] = {3, 4, 7};
    COMPLEX q[7], computed[2];
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    UINT degree_index, i;
    INT kappa, normalization_flag, ret_code;

    for (i = 0; i < D; i++)
        q[i] = 0.18 + 0.07*CEXP(6.2831853071795864769*I*i/D);
    opts.discretization = nse_discretization_FES8_PADE;
    opts.contspec_type = nsev_cstype_REFLECTION_COEFFICIENT;
    for (normalization_flag = 0; normalization_flag <= 1;
            normalization_flag++) {
        opts.normalization_flag = normalization_flag;
        for (kappa = -1; kappa <= 1; kappa += 2) {
            for (degree_index = 0;
                    degree_index < sizeof(degrees)/sizeof(degrees[0]);
                    degree_index++) {
                opts.pade_degree = degrees[degree_index];
                /* Keep the degree-3 family default under test. For the
                 * higher-degree characterization, a smaller allowed scale
                 * keeps w farther from one and avoids the documented global
                 * power-basis cancellation on this tiny grid. */
                opts.pade_h = opts.pade_degree == 3 ? 0.0 : 5.0;
                ret_code = fnft_nsev(D, q, T, M, computed, XI, NULL, NULL,
                        NULL, kappa, &opts);
                if (ret_code != SUCCESS)
                    return ret_code;
                for (i = 0; i < M; i++) {
                    const COMPLEX expected = pointwise_reflection(D, q, T,
                            XI[0] + (XI[1] - XI[0])*i/(M - 1), kappa,
                            opts.pade_degree);
                    if (CABS(computed[i] - expected)
                            > 2e-9*(1.0 + CABS(expected))) {
                        fprintf(stderr, "FES8 public small-grid mismatch: degree=%lu kappa=%d norm=%d i=%lu error=%.3e\n",
                                (unsigned long)opts.pade_degree, (int)kappa,
                                (int)normalization_flag, (unsigned long)i,
                                CABS(computed[i] - expected));
                        return E_TEST_FAILED;
                    }
                }
            }
        }
    }
    return SUCCESS;
}

static INT test_public_es8_chebyshev(void)
{
    const UINT D = 7, M = 5;
    const REAL T[2] = {-0.21, 0.21};
    const REAL XI[2] = {-0.83, 1.27};
    const UINT degrees[] = {3, 4, 5, 6, 7};
    COMPLEX q[7], computed[5];
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    UINT degree_index, i;
    INT kappa, normalization_flag, ret_code;

    for (i = 0; i < D; i++)
        q[i] = 0.18 + 0.07*CEXP(6.2831853071795864769*I*i/D);
    opts.discretization = nse_discretization_FES8_PADE;
    opts.contspec_type = nsev_cstype_REFLECTION_COEFFICIENT;
    opts.pade_representation =
            nsev_pade_representation_CHEBYSHEV_JOUKOWSKI;
    for (normalization_flag = 0; normalization_flag <= 1;
            normalization_flag++) {
        opts.normalization_flag = normalization_flag;
        for (kappa = -1; kappa <= 1; kappa += 2) {
            for (degree_index = 0;
                    degree_index < sizeof(degrees)/sizeof(degrees[0]);
                    degree_index++) {
                opts.pade_degree = degrees[degree_index];
                ret_code = fnft_nsev(D, q, T, M, computed, XI, NULL, NULL,
                        NULL, kappa, &opts);
                if (ret_code != SUCCESS)
                    return ret_code;
                for (i = 0; i < M; i++) {
                    const COMPLEX lambda = XI[0]
                            + (XI[1] - XI[0])*i/(M - 1);
                    const COMPLEX expected = pointwise_reflection(D, q, T,
                            lambda, kappa, opts.pade_degree);
                    const REAL tolerance = 4e-8*(1.0 + CABS(expected));

                    if (CABS(computed[i] - expected) > tolerance) {
                        fprintf(stderr, "FES8 public Chebyshev mismatch: degree=%lu kappa=%d norm=%d i=%lu error=%.3e\n",
                                (unsigned long)opts.pade_degree, (int)kappa,
                                (int)normalization_flag, (unsigned long)i,
                                CABS(computed[i] - expected));
                        return E_TEST_FAILED;
                    }
                }
            }
        }
    }
    return SUCCESS;
}

static INT test_public_es8_rejections(void)
{
    const UINT D = 8, M = 2;
    const REAL T[2] = {-0.5, 0.5};
    const REAL XI[2] = {-1.0, 1.0};
    COMPLEX q[8] = {0.0}, q_short[6] = {0.0};
    COMPLEX spectrum[2], bound_state[1];
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    UINT K = 1;

    opts.discretization = nse_discretization_FES8_PADE;
    opts.pade_degree = 2;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_degree = 8;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_degree = 3;
    opts.pade_h = -1.0;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_h = NAN;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_h = 0.0;
    if (fnft_nsev(6, q_short, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    K = 1;
    if (fnft_nsev(D, q, T, 0, NULL, NULL, &K, NULL, NULL, +1,
                &opts) != SUCCESS || K != 0)
        return E_TEST_FAILED;
    K = 1;
    if (fnft_nsev(D, q, T, 0, NULL, NULL, &K, bound_state, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_representation =
            nsev_pade_representation_CHEBYSHEV_JOUKOWSKI;
    opts.pade_h = 1.0;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_h = 0.0;
    opts.discretization = nse_discretization_FES6_PADE;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    opts.pade_representation = (fnft_nsev_pade_representation_t)42;
    if (fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL, +1,
                &opts) == SUCCESS)
        return E_TEST_FAILED;
    return SUCCESS;
}

static INT test_es8_metadata(void)
{
    if (nse_discretization_method_order(nse_discretization_FES8_PADE) != 8
            || nse_discretization_effective_order(
                nse_discretization_FES8_PADE, 3) != 6
            || nse_discretization_effective_order(
                nse_discretization_FES8_PADE, 4) != 8
            || nse_discretization_effective_order(
                nse_discretization_FES8_PADE, 7) != 8
            || nse_discretization_effective_order(
                nse_discretization_FES8_PADE, 2) != 0
            || FABS(nse_discretization_pade_h(
                nse_discretization_FES8_PADE, 3, 0.0) - 14.9) > EPSILON
            || FABS(nse_discretization_pade_h(
                nse_discretization_FES8_PADE, 4, 0.0) - 19.4) > EPSILON
            || FABS(nse_discretization_pade_h(
                nse_discretization_FES8_PADE, 7, 0.0) - 21.8) > EPSILON)
        return E_TEST_FAILED;
    {
        fnft_nsev_opts_t opts = fnft_nsev_default_opts();
        if (opts.pade_representation
                != nsev_pade_representation_DIRECT_CAYLEY)
            return E_TEST_FAILED;
        opts.discretization = nse_discretization_FES8_PADE;
        if (fnft_nsev_max_K(128, &opts) != 0)
            return E_TEST_FAILED;
    }
    {
        const REAL eps_t = 0.17, h = 19.4;
        const COMPLEX lambda = 0.63 + 0.21*I;
        COMPLEX generic_value = lambda, cayley_value = lambda;
        akns_discretization_t akns_discretization;

        if (nse_discretization_to_akns_discretization(
                    nse_discretization_FES8_PADE,
                    &akns_discretization) == SUCCESS
                || nse_discretization_lambda_to_z(1, eps_t, &generic_value,
                    nse_discretization_FES8_PADE) == SUCCESS
                || generic_value != lambda)
            return E_TEST_FAILED;
        generic_value = lambda;
        if (nse_discretization_z_to_lambda(1, eps_t, &generic_value,
                    nse_discretization_FES8_PADE) == SUCCESS
                || generic_value != lambda)
            return E_TEST_FAILED;
        if (nse_discretization_pade_lambda_to_z(1, eps_t, &cayley_value, h)
                    != SUCCESS
                || nse_discretization_pade_z_to_lambda(1, eps_t,
                    &cayley_value, h) != SUCCESS
                /* The inverse subtracts w from one; for this w close to one,
                 * the measured roundtrip error is 18 scaled ulps. */
                || CABS(cayley_value - lambda)
                    > 64.0*EPSILON*(1.0 + CABS(lambda)))
            return E_TEST_FAILED;
    }
    return SUCCESS;
}

static INT test_es8_pointwise_convergence(const nsev_testcases_t testcase)
{
    REAL errors[2][3], slow_error;
    const UINT degrees[2] = {3, 4};
    const UINT grids[3] = {128, 256, 512};
    UINT degree_index, grid_index;
    INT ret_code;

    for (degree_index = 0; degree_index < 2; degree_index++) {
        for (grid_index = 0; grid_index < 3; grid_index++) {
            errors[degree_index][grid_index] = pointwise_spectrum_error(
                    testcase, grids[grid_index], degrees[degree_index],
                    &ret_code);
            if (ret_code != SUCCESS)
                return ret_code;
        }
    }
    if (!(errors[0][0]/errors[0][1] >= 60.0
            && errors[0][1]/errors[0][2] >= 60.0
            && errors[0][2] < 1e-5)) {
        fprintf(stderr, "FES8 pointwise order-six failure: errors=%.3e,%.3e,%.3e ratios=%.2f,%.2f\n",
                errors[0][0], errors[0][1], errors[0][2],
                errors[0][0]/errors[0][1], errors[0][1]/errors[0][2]);
        return E_TEST_FAILED;
    }
    if (!(errors[1][0]/errors[1][1] >= 90.0
            && errors[1][1]/errors[1][2] >= 140.0
            && errors[1][2] < 1e-5)) {
        fprintf(stderr, "FES8 pointwise order-eight failure: errors=%.3e,%.3e,%.3e ratios=%.2f,%.2f\n",
                errors[1][0], errors[1][1], errors[1][2],
                errors[1][0]/errors[1][1], errors[1][1]/errors[1][2]);
        return E_TEST_FAILED;
    }
    slow_error = continuous_spectrum_error(testcase, grids[2],
            nse_discretization_ES8, 0,
            nsev_pade_representation_DIRECT_CAYLEY, &ret_code);
    if (ret_code != SUCCESS)
        return ret_code;
    if (FABS(errors[1][2] - slow_error)
            > 0.02*slow_error + 1e-12) {
        fprintf(stderr, "FES8 pointwise/slow mismatch: Pade=%.3e slow=%.3e\n",
                errors[1][2], slow_error);
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

static INT test_convergence(const nsev_testcases_t testcase,
        const nse_discretization_t discretization, const UINT pade_degree,
        const UINT coarse_D, const REAL minimum_ratio)
{
    INT ret_code;
    const REAL coarse = continuous_spectrum_error(testcase, coarse_D,
            discretization, pade_degree,
            nsev_pade_representation_DIRECT_CAYLEY, &ret_code);
    REAL fine;

    if (ret_code != SUCCESS)
        return ret_code;
    fine = continuous_spectrum_error(testcase, 2*coarse_D, discretization,
            pade_degree, nsev_pade_representation_DIRECT_CAYLEY, &ret_code);
    if (ret_code != SUCCESS)
        return ret_code;
    if (!(fine > 0.0 && coarse/fine >= minimum_ratio)) {
        fprintf(stderr, "Padé convergence failure: degree=%lu coarse=%.3e fine=%.3e ratio=%.2f\n",
                (unsigned long)pade_degree, coarse, fine, coarse/fine);
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

static INT test_chebyshev_convergence(const nsev_testcases_t testcase)
{
    const UINT degrees[] = {3, 4};
    const REAL minimum_ratios[] = {45.0, 90.0};
    UINT degree_index;
    INT ret_code;

    for (degree_index = 0;
            degree_index < sizeof(degrees)/sizeof(degrees[0]);
            degree_index++) {
        const REAL coarse = continuous_spectrum_error(testcase, 128,
                nse_discretization_FES8_PADE, degrees[degree_index],
                nsev_pade_representation_CHEBYSHEV_JOUKOWSKI, &ret_code);
        REAL fine;

        if (ret_code != SUCCESS)
            return ret_code;
        fine = continuous_spectrum_error(testcase, 256,
                nse_discretization_FES8_PADE, degrees[degree_index],
                nsev_pade_representation_CHEBYSHEV_JOUKOWSKI, &ret_code);
        if (ret_code != SUCCESS)
            return ret_code;
        if (!(fine > 0.0 && coarse/fine >= minimum_ratios[degree_index])) {
            fprintf(stderr, "FES8 Chebyshev convergence failure: degree=%lu coarse=%.3e fine=%.3e ratio=%.2f\n",
                    (unsigned long)degrees[degree_index], coarse, fine,
                    coarse/fine);
            return E_TEST_FAILED;
        }
    }
    return SUCCESS;
}

static REAL bound_state_error(const UINT D, INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL, *computed = NULL;
    REAL T[2], XI[2], error = 0.0;
    UINT M, K, K_exact, i, j;
    INT kappa;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();

    *ret_code = nsev_testcases(nsev_testcases_SECH_FOCUSING, D, &q, T,
            &M, &exact, &ab, XI, &K_exact, &bound_states, &normconsts,
            &residues, &kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    opts.discretization = nse_discretization_FES4_PADE;
    opts.pade_degree = 2;
    opts.bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE;
    K = fnft_nsev_max_K(D, &opts);
    computed = malloc(K*sizeof(COMPLEX));
    if (computed == NULL) {
        *ret_code = E_NOMEM;
        goto release_mem;
    }
    *ret_code = fnft_nsev(D, q, T, 0, NULL, XI, &K, computed, NULL,
            kappa, &opts);
    if (*ret_code != SUCCESS)
        goto release_mem;
    if (K != K_exact) {
        fprintf(stderr, "Padé bound-state count failure: D=%lu computed=%lu exact=%lu\n",
                (unsigned long)D, (unsigned long)K, (unsigned long)K_exact);
        *ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    for (i = 0; i < K; i++) {
        REAL minimum = INFINITY;
        for (j = 0; j < K_exact; j++) {
            const REAL distance = FNFT_CABS(computed[i] - bound_states[j]);
            if (distance < minimum)
                minimum = distance;
        }
        if (minimum > error)
            error = minimum;
    }

release_mem:
    free(q);
    free(exact);
    free(ab);
    free(bound_states);
    free(normconsts);
    free(residues);
    free(computed);
    return error;
}

static INT test_bound_state_localization(void)
{
    INT ret_code;
    const REAL coarse = bound_state_error(256, &ret_code);
    REAL fine;

    if (ret_code != SUCCESS)
        return ret_code;
    fine = bound_state_error(512, &ret_code);
    if (ret_code != SUCCESS)
        return ret_code;
    if (!(fine > 0.0 && fine < 1e-3 && coarse/fine >= 3.0)) {
        fprintf(stderr, "Padé bound-state localization failure: coarse=%.3e fine=%.3e ratio=%.2f\n",
                coarse, fine, coarse/fine);
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

INT main(void)
{
    INT ret_code, kappa, normalization_flag;

    ret_code = test_es8_metadata();
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_public_es8_small_grid();
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_public_es8_chebyshev();
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_public_es8_rejections();
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_es8_pointwise_convergence(
            nsev_testcases_SECH_FOCUSING_CONTSPEC);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_es8_pointwise_convergence(
            nsev_testcases_SECH_DEFOCUSING);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_chebyshev_convergence(
            nsev_testcases_SECH_FOCUSING_CONTSPEC);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_chebyshev_convergence(nsev_testcases_SECH_DEFOCUSING);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    for (normalization_flag = 0; normalization_flag <= 1;
            normalization_flag++) {
        for (kappa = -1; kappa <= 1; kappa += 2) {
            ret_code = test_zero_signal(nse_discretization_FES4_PADE, 2,
                    kappa, normalization_flag,
                    nsev_pade_representation_DIRECT_CAYLEY);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            ret_code = test_zero_signal(nse_discretization_FES6_PADE, 4,
                    kappa, normalization_flag,
                    nsev_pade_representation_DIRECT_CAYLEY);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            ret_code = test_zero_signal(nse_discretization_FES8_PADE, 3,
                    kappa, normalization_flag,
                    nsev_pade_representation_DIRECT_CAYLEY);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            ret_code = test_zero_signal(nse_discretization_FES8_PADE, 7,
                    kappa, normalization_flag,
                    nsev_pade_representation_CHEBYSHEV_JOUKOWSKI);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            if (normalization_flag == 0 && kappa == 1) {
                ret_code = test_zero_signal(nse_discretization_FES6_PADE, 7,
                        kappa, normalization_flag,
                        nsev_pade_representation_DIRECT_CAYLEY);
                if (ret_code != SUCCESS)
                    return EXIT_FAILURE;
            }
        }
    }

    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,
            nse_discretization_FES4_PADE, 2, 512, 8.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_convergence(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            nse_discretization_FES6_PADE, 3, 128, 24.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,
            nse_discretization_FES6_PADE, 4, 512, 24.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_bound_state_localization();
    return ret_code == SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
