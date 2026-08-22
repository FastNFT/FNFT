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

#include "fnft__akns_fscatter_pade.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    FNFT_COMPLEX entry[4];
} matrix_t;

static matrix_t matrix_zero(void)
{
    const matrix_t result = {{0.0, 0.0, 0.0, 0.0}};
    return result;
}

static matrix_t matrix_identity(void)
{
    const matrix_t result = {{1.0, 0.0, 0.0, 1.0}};
    return result;
}

static void matrix_add_scaled(matrix_t * const result,
        matrix_t const * const term, const FNFT_COMPLEX scale)
{
    FNFT_UINT i;

    for (i = 0; i < 4; i++)
        result->entry[i] += scale*term->entry[i];
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

static matrix_t matrix_commutator(matrix_t const * const left,
        matrix_t const * const right)
{
    matrix_t result = matrix_multiply(left, right);
    const matrix_t reverse = matrix_multiply(right, left);

    matrix_add_scaled(&result, &reverse, -1.0);
    return result;
}

static matrix_t derivative_matrix(const FNFT_COMPLEX q,
        const FNFT_COMPLEX r)
{
    matrix_t result = matrix_zero();
    result.entry[1] = q;
    result.entry[2] = r;
    return result;
}

static matrix_t direct_z(const FNFT_UINT D, const FNFT_UINT index,
        FNFT_COMPLEX const * const q, FNFT_COMPLEX const * const r,
        const FNFT_REAL eps_t, const FNFT_COMPLEX lambda,
        const FNFT_UINT method_order)
{
    const FNFT_UINT im1 = (index + D - 1) % D;
    const FNFT_UINT ip1 = (index + 1) % D;
    const FNFT_UINT im2 = (index + D - 2) % D;
    const FNFT_UINT ip2 = (index + 2) % D;
    matrix_t a0 = derivative_matrix(eps_t*q[index], eps_t*r[index]);
    matrix_t d1, d2, term;
    matrix_t result;

    a0.entry[0] = -I*eps_t*lambda;
    a0.entry[3] = I*eps_t*lambda;

    if (method_order == 4) {
        d1 = derivative_matrix(eps_t*(q[ip1] - q[im1])/2.0,
                eps_t*(r[ip1] - r[im1])/2.0);
        d2 = derivative_matrix(eps_t*(q[ip1] - 2.0*q[index] + q[im1]),
                eps_t*(r[ip1] - 2.0*r[index] + r[im1]));
        result = a0;
        matrix_add_scaled(&result, &d2, 1.0/24.0);
        term = matrix_commutator(&d1, &a0);
        matrix_add_scaled(&result, &term, 1.0/12.0);
    } else {
        matrix_t d1_low, d2_low, d3, d4, t1, t2, a0_squared, a0_cubed;

        d1 = derivative_matrix(eps_t*(-q[ip2] + 8.0*q[ip1]
                    - 8.0*q[im1] + q[im2])/12.0,
                eps_t*(-r[ip2] + 8.0*r[ip1]
                    - 8.0*r[im1] + r[im2])/12.0);
        d2 = derivative_matrix(eps_t*(-q[ip2] + 16.0*q[ip1]
                    - 30.0*q[index] + 16.0*q[im1] - q[im2])/12.0,
                eps_t*(-r[ip2] + 16.0*r[ip1]
                    - 30.0*r[index] + 16.0*r[im1] - r[im2])/12.0);
        result = a0;
        matrix_add_scaled(&result, &d2, 1.0/24.0);
        term = matrix_commutator(&d1, &a0);
        matrix_add_scaled(&result, &term, 1.0/12.0);

        d1_low = derivative_matrix(eps_t*(q[ip1] - q[im1])/2.0,
                eps_t*(r[ip1] - r[im1])/2.0);
        d2_low = derivative_matrix(eps_t*(q[ip1] - 2.0*q[index] + q[im1]),
                eps_t*(r[ip1] - 2.0*r[index] + r[im1]));
        d3 = derivative_matrix(eps_t*(q[ip2] - 2.0*q[ip1]
                    + 2.0*q[im1] - q[im2])/2.0,
                eps_t*(r[ip2] - 2.0*r[ip1]
                    + 2.0*r[im1] - r[im2])/2.0);
        d4 = derivative_matrix(eps_t*(q[ip2] - 4.0*q[ip1]
                    + 6.0*q[index] - 4.0*q[im1] + q[im2]),
                eps_t*(r[ip2] - 4.0*r[ip1]
                    + 6.0*r[index] - 4.0*r[im1] + r[im2]));
        matrix_add_scaled(&result, &d4, 1.0/1920.0);
        term = matrix_commutator(&d3, &a0);
        matrix_add_scaled(&result, &term, 1.0/480.0);
        term = matrix_commutator(&d1_low, &d2_low);
        matrix_add_scaled(&result, &term, 1.0/480.0);
        t1 = matrix_commutator(&d2_low, &a0);
        term = matrix_commutator(&t1, &a0);
        matrix_add_scaled(&result, &term, 1.0/720.0);
        t1 = matrix_commutator(&a0, &d1_low);
        term = matrix_commutator(&t1, &d1_low);
        matrix_add_scaled(&result, &term, 1.0/240.0);
        a0_squared = matrix_multiply(&a0, &a0);
        a0_cubed = matrix_multiply(&a0_squared, &a0);
        term = matrix_commutator(&a0_cubed, &d1_low);
        matrix_add_scaled(&result, &term, 1.0/720.0);
        t1 = matrix_multiply(&a0, &d1_low);
        t2 = matrix_multiply(&t1, &a0);
        term = matrix_commutator(&t2, &a0);
        matrix_add_scaled(&result, &term, 1.0/240.0);
    }
    return result;
}

static matrix_t direct_pade(matrix_t const * const z,
        const FNFT_UINT pade_degree)
{
    matrix_t p_plus = matrix_identity();
    matrix_t p_minus = matrix_identity();
    matrix_t power = matrix_identity();
    matrix_t inverse_minus, result;
    FNFT_REAL coefficient = 1.0;
    FNFT_COMPLEX determinant;
    FNFT_UINT j;

    for (j = 1; j <= pade_degree; j++) {
        power = matrix_multiply(&power, z);
        coefficient *= (FNFT_REAL)(pade_degree + 1 - j)
                / ((FNFT_REAL)j*(FNFT_REAL)(2*pade_degree + 1 - j));
        matrix_add_scaled(&p_plus, &power, coefficient);
        matrix_add_scaled(&p_minus, &power,
                (j & 1U) == 0 ? coefficient : -coefficient);
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

static FNFT_COMPLEX evaluate_polynomial(FNFT_COMPLEX const * const p,
        const FNFT_UINT degree, const FNFT_COMPLEX z)
{
    FNFT_COMPLEX result = p[0];
    FNFT_UINT i;

    for (i = 1; i <= degree; i++)
        result = result*z + p[i];
    return result;
}

static int run_case(const FNFT_UINT method_order, const FNFT_UINT pade_degree,
        const FNFT_INT kappa, const FNFT_INT normalization_flag)
{
    const FNFT_UINT D = 5;
    const FNFT_REAL eps_t = 0.07;
    const FNFT_REAL h = method_order == 4 ? 5.0 : 11.0;
    const FNFT_COMPLEX lambda = 0.63 + 0.17*I;
    const FNFT_COMPLEX q[5] = {
        0.31 + 0.08*I, -0.17 + 0.23*I, 0.42 - 0.19*I,
        -0.28 - 0.11*I, 0.09 + 0.37*I
    };
    FNFT_COMPLEX r[5];
    FNFT_COMPLEX *numerator;
    FNFT_COMPLEX *denominator;
    FNFT_UINT numerator_degree, denominator_degree, i, component;
    FNFT_INT numerator_exponent = 0, denominator_exponent = 0;
    FNFT_INT ret_code;
    matrix_t direct = matrix_identity();
    matrix_t computed;
    FNFT_COMPLEX w, denominator_value;
    FNFT_REAL scale, error = 0.0, reference_norm = 0.0;

    numerator = malloc(fnft__akns_fscatter_pade_numel(D, method_order,
                pade_degree)*sizeof(FNFT_COMPLEX));
    denominator = malloc(fnft__akns_fscatter_pade_den_numel(D, method_order,
                pade_degree)*sizeof(FNFT_COMPLEX));
    if (numerator == NULL || denominator == NULL)
        return EXIT_FAILURE;

    for (i = 0; i < D; i++)
        r[i] = -(FNFT_REAL)kappa*FNFT_CONJ(q[i]);

    ret_code = fnft__akns_fscatter_pade(D, q, r, eps_t, method_order,
            pade_degree, h, numerator, &numerator_degree,
            normalization_flag ? &numerator_exponent : NULL,
            denominator, &denominator_degree,
            normalization_flag ? &denominator_exponent : NULL);
    if (ret_code != FNFT_SUCCESS) {
        free(numerator);
        free(denominator);
        return EXIT_FAILURE;
    }

    for (i = 0; i < D; i++) {
        const matrix_t z = direct_z(D, i, q, r, eps_t, lambda, method_order);
        const matrix_t step = direct_pade(&z, pade_degree);
        direct = matrix_multiply(&step, &direct);
    }

    w = (I*h - eps_t*lambda)/(I*h + eps_t*lambda);
    denominator_value = evaluate_polynomial(denominator,
            denominator_degree, w);
    scale = ldexp(1.0, numerator_exponent - denominator_exponent);
    for (component = 0; component < 4; component++) {
        computed.entry[component] = scale*evaluate_polynomial(
                numerator + component*(numerator_degree + 1),
                numerator_degree, w)/denominator_value;
        error += FNFT_CABS(computed.entry[component] - direct.entry[component]);
        reference_norm += FNFT_CABS(direct.entry[component]);
    }

    free(numerator);
    free(denominator);
    if (error > 2e-9*(1.0 + reference_norm)) {
        fprintf(stderr, "Padé mismatch: order=%lu degree=%lu kappa=%d norm=%d error=%.3e\n",
                (unsigned long)method_order, (unsigned long)pade_degree,
                (int)kappa, (int)normalization_flag, error);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int run_zero_case(const FNFT_INT normalization_flag)
{
    const FNFT_UINT D = 16, method_order = 4, pade_degree = 2;
    const FNFT_REAL eps_t = 2.0/15.0, h = 3.4641016151377546;
    const FNFT_COMPLEX lambda = -2.0;
    FNFT_COMPLEX q[16] = {0.0}, r[16] = {0.0};
    FNFT_COMPLEX *numerator, *denominator, w, numerator_value;
    FNFT_UINT numerator_degree, denominator_degree;
    FNFT_INT numerator_exponent = 0, denominator_exponent = 0, ret_code;
    FNFT_REAL scale;

    numerator = malloc(fnft__akns_fscatter_pade_numel(D, method_order,
                pade_degree)*sizeof(FNFT_COMPLEX));
    denominator = malloc(fnft__akns_fscatter_pade_den_numel(D, method_order,
                pade_degree)*sizeof(FNFT_COMPLEX));
    if (numerator == NULL || denominator == NULL)
        return EXIT_FAILURE;
    ret_code = fnft__akns_fscatter_pade(D, q, r, eps_t, method_order,
            pade_degree, h, numerator, &numerator_degree,
            normalization_flag ? &numerator_exponent : NULL, denominator,
            &denominator_degree,
            normalization_flag ? &denominator_exponent : NULL);
    if (ret_code != FNFT_SUCCESS)
        return EXIT_FAILURE;
    w = (I*h - eps_t*lambda)/(I*h + eps_t*lambda);
    numerator_value = evaluate_polynomial(numerator, numerator_degree, w);
    scale = ldexp(1.0, numerator_exponent - denominator_exponent);
    numerator_value *= scale/evaluate_polynomial(denominator,
            denominator_degree, w);
    free(numerator);
    free(denominator);
    if (FNFT_CABS(numerator_value) < 0.9
            || FNFT_CABS(numerator_value) > 1.1) {
        fprintf(stderr, "Padé zero core failure: norm=%d abs=%.3e\n",
                (int)normalization_flag, FNFT_CABS(numerator_value));
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int main(void)
{
    static const FNFT_UINT fourth_degrees[] = {2, 3, 4, 7};
    static const FNFT_UINT sixth_degrees[] = {3, 4, 7};
    FNFT_UINT i;
    FNFT_INT kappa, normalization_flag;

    if (run_zero_case(0) != EXIT_SUCCESS
            || run_zero_case(1) != EXIT_SUCCESS)
        return EXIT_FAILURE;

    for (normalization_flag = 0; normalization_flag <= 1; normalization_flag++) {
        for (kappa = -1; kappa <= 1; kappa += 2) {
            for (i = 0; i < sizeof(fourth_degrees)/sizeof(fourth_degrees[0]); i++) {
                if (run_case(4, fourth_degrees[i], kappa,
                            normalization_flag) != EXIT_SUCCESS)
                    return EXIT_FAILURE;
            }
            for (i = 0; i < sizeof(sixth_degrees)/sizeof(sixth_degrees[0]); i++) {
                if (run_case(6, sixth_degrees[i], kappa,
                            normalization_flag) != EXIT_SUCCESS)
                    return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
