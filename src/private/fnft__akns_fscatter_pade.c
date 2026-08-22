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

#include "fnft__akns_fscatter_pade.h"
#include "fnft__akns_discretization.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

#include <limits.h>
#include <stdlib.h>

enum {
    PADE_MAX_DEGREE = 7,
    Z_MAX_DEGREE = 5,
    LOCAL_MAX_DEGREE = 2*Z_MAX_DEGREE*PADE_MAX_DEGREE
};

typedef struct {
    UINT degree;
    COMPLEX coefficient[Z_MAX_DEGREE + 1];
} small_poly_t;

typedef struct {
    small_poly_t entry[4];
} small_matrix_t;

static void small_poly_zero(small_poly_t * const p)
{
    UINT i;

    p->degree = 0;
    for (i = 0; i <= Z_MAX_DEGREE; i++)
        p->coefficient[i] = 0.0;
}

static void small_poly_trim(small_poly_t * const p)
{
    while (p->degree > 0 && p->coefficient[p->degree] == 0.0)
        p->degree--;
}

static void small_poly_add_scaled(small_poly_t * const result,
        small_poly_t const * const term, const COMPLEX scale)
{
    UINT i;

    if (result->degree < term->degree)
        result->degree = term->degree;
    for (i = 0; i <= term->degree; i++)
        result->coefficient[i] += scale*term->coefficient[i];
    small_poly_trim(result);
}

static void small_poly_multiply(small_poly_t const * const left,
        small_poly_t const * const right, small_poly_t * const result)
{
    UINT i, j;

    small_poly_zero(result);
    result->degree = left->degree + right->degree;
    for (i = 0; i <= left->degree; i++) {
        for (j = 0; j <= right->degree; j++)
            result->coefficient[i + j] +=
                    left->coefficient[i]*right->coefficient[j];
    }
    small_poly_trim(result);
}

static void small_matrix_zero(small_matrix_t * const matrix)
{
    UINT i;

    for (i = 0; i < 4; i++)
        small_poly_zero(&matrix->entry[i]);
}

static void small_matrix_add_scaled(small_matrix_t * const result,
        small_matrix_t const * const term, const COMPLEX scale)
{
    UINT i;

    for (i = 0; i < 4; i++)
        small_poly_add_scaled(&result->entry[i], &term->entry[i], scale);
}

static void small_matrix_multiply(small_matrix_t const * const left,
        small_matrix_t const * const right, small_matrix_t * const result)
{
    small_poly_t t1, t2;
    UINT row, column;

    small_matrix_zero(result);
    for (row = 0; row < 2; row++) {
        for (column = 0; column < 2; column++) {
            small_poly_multiply(&left->entry[2*row],
                    &right->entry[column], &t1);
            small_poly_multiply(&left->entry[2*row + 1],
                    &right->entry[2 + column], &t2);
            small_poly_add_scaled(&result->entry[2*row + column], &t1, 1.0);
            small_poly_add_scaled(&result->entry[2*row + column], &t2, 1.0);
        }
    }
}

static void small_matrix_commutator(small_matrix_t const * const left,
        small_matrix_t const * const right, small_matrix_t * const result)
{
    small_matrix_t lr, rl;

    small_matrix_multiply(left, right, &lr);
    small_matrix_multiply(right, left, &rl);
    *result = lr;
    small_matrix_add_scaled(result, &rl, -1.0);
}

static void constant_derivative_matrix(const COMPLEX q_value,
        const COMPLEX r_value, small_matrix_t * const matrix)
{
    small_matrix_zero(matrix);
    matrix->entry[1].coefficient[0] = q_value;
    matrix->entry[2].coefficient[0] = r_value;
}

static void build_z_polynomial(const UINT D, const UINT index,
        COMPLEX const * const q, COMPLEX const * const r,
        const REAL eps_t, const UINT method_order,
        small_matrix_t * const z_matrix)
{
    const UINT im1 = (index + D - 1) % D;
    const UINT ip1 = (index + 1) % D;
    const UINT im2 = (index + D - 2) % D;
    const UINT ip2 = (index + 2) % D;
    small_matrix_t a0, d1, d2, comm;

    constant_derivative_matrix(eps_t*q[index], eps_t*r[index], &a0);
    a0.entry[0].degree = 1;
    a0.entry[0].coefficient[1] = -I;
    a0.entry[3].degree = 1;
    a0.entry[3].coefficient[1] = I;

    if (method_order == 4) {
        constant_derivative_matrix(eps_t*(q[ip1] - q[im1])/2.0,
                eps_t*(r[ip1] - r[im1])/2.0, &d1);
        constant_derivative_matrix(eps_t*(q[ip1] - 2.0*q[index] + q[im1]),
                eps_t*(r[ip1] - 2.0*r[index] + r[im1]), &d2);
        small_matrix_commutator(&d1, &a0, &comm);

        *z_matrix = a0;
        small_matrix_add_scaled(z_matrix, &d2, 1.0/24.0);
        small_matrix_add_scaled(z_matrix, &comm, 1.0/12.0);
    } else if (method_order == 6) {
        const COMPLEX q_samples[5] = {
            q[im2], q[im1], q[index], q[ip1], q[ip2]
        };
        const COMPLEX r_samples[5] = {
            r[im2], r[im1], r[index], r[ip1], r[ip2]
        };
        fnft__akns_es6_stencil_t qs, rs;
        small_matrix_t d1_low, d2_low, d3, d4;
        small_matrix_t term, t1, t2, a0_squared, a0_cubed;

        akns_discretization_es6_stencil(q_samples, eps_t, &qs);
        akns_discretization_es6_stencil(r_samples, eps_t, &rs);
        constant_derivative_matrix(qs.first, rs.first, &d1);
        constant_derivative_matrix(qs.second, rs.second, &d2);
        small_matrix_commutator(&d1, &a0, &comm);
        *z_matrix = a0;
        small_matrix_add_scaled(z_matrix, &d2, 1.0/24.0);
        small_matrix_add_scaled(z_matrix, &comm, 1.0/12.0);

        constant_derivative_matrix(qs.first_low, rs.first_low, &d1_low);
        constant_derivative_matrix(qs.second_low, rs.second_low, &d2_low);
        constant_derivative_matrix(qs.third, rs.third, &d3);
        constant_derivative_matrix(qs.fourth, rs.fourth, &d4);

        small_matrix_add_scaled(z_matrix, &d4, 1.0/1920.0);
        small_matrix_commutator(&d3, &a0, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/480.0);
        small_matrix_commutator(&d1_low, &d2_low, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/480.0);
        small_matrix_commutator(&d2_low, &a0, &t1);
        small_matrix_commutator(&t1, &a0, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/720.0);
        small_matrix_commutator(&a0, &d1_low, &t1);
        small_matrix_commutator(&t1, &d1_low, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/240.0);
        small_matrix_multiply(&a0, &a0, &a0_squared);
        small_matrix_multiply(&a0_squared, &a0, &a0_cubed);
        small_matrix_commutator(&a0_cubed, &d1_low, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/720.0);
        small_matrix_multiply(&a0, &d1_low, &t1);
        small_matrix_multiply(&t1, &a0, &t2);
        small_matrix_commutator(&t2, &a0, &term);
        small_matrix_add_scaled(z_matrix, &term, 1.0/240.0);
    } else {
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

        small_matrix_zero(z_matrix);
        for (i = 0; i < 4; i++) {
            z_matrix->entry[i].degree = Z_MAX_DEGREE;
            for (j = 0; j <= Z_MAX_DEGREE; j++)
                z_matrix->entry[i].coefficient[j] = coefficients[4*j + i];
            small_poly_trim(&z_matrix->entry[i]);
        }
    }
}

static REAL binomial_coefficient(const UINT n, const UINT k)
{
    REAL value = 1.0;
    UINT i;

    for (i = 1; i <= k; i++)
        value *= (REAL)(n + 1 - i)/(REAL)i;
    return value;
}

static void transform_z_to_w(small_poly_t const * const z_polynomial,
        const UINT z_degree, const REAL h,
        COMPLEX w_polynomial[Z_MAX_DEGREE + 1])
{
    UINT j, left, right;

    for (j = 0; j <= Z_MAX_DEGREE; j++)
        w_polynomial[j] = 0.0;

    for (j = 0; j <= z_polynomial->degree; j++) {
        const COMPLEX scale = z_polynomial->coefficient[j]*CPOW(I*h, j);
        for (left = 0; left <= j; left++) {
            const REAL c_left = binomial_coefficient(j, left)
                    * ((left & 1U) == 0 ? 1.0 : -1.0);
            for (right = 0; right <= z_degree - j; right++) {
                w_polynomial[left + right] += scale*c_left
                        * binomial_coefficient(z_degree - j, right);
            }
        }
    }
}

static void polynomial_zero(COMPLEX p[LOCAL_MAX_DEGREE + 1])
{
    UINT i;

    for (i = 0; i <= LOCAL_MAX_DEGREE; i++)
        p[i] = 0.0;
}

static void polynomial_multiply(COMPLEX const a[LOCAL_MAX_DEGREE + 1],
        const UINT degree_a, COMPLEX const b[LOCAL_MAX_DEGREE + 1],
        const UINT degree_b, COMPLEX result[LOCAL_MAX_DEGREE + 1])
{
    UINT i, j;

    polynomial_zero(result);
    for (i = 0; i <= degree_a; i++) {
        for (j = 0; j <= degree_b; j++)
            result[i + j] += a[i]*b[j];
    }
}

static void polynomial_add_padded_power(
        COMPLEX result[LOCAL_MAX_DEGREE + 1],
        COMPLEX const lambda_power[LOCAL_MAX_DEGREE + 1],
        const UINT lambda_degree, const UINT padding_power,
        const COMPLEX scale)
{
    UINT i, j;

    for (i = 0; i <= lambda_degree; i++) {
        for (j = 0; j <= padding_power; j++)
            result[i + j] += scale*lambda_power[i]
                    * binomial_coefficient(padding_power, j);
    }
}

static void pade_coefficients(const UINT degree,
        REAL f0[PADE_MAX_DEGREE + 1],
        REAL temp[PADE_MAX_DEGREE],
        REAL denominator[PADE_MAX_DEGREE + 1])
{
    REAL p[PADE_MAX_DEGREE + 1] = {0.0};
    REAL square[2*PADE_MAX_DEGREE + 1] = {0.0};
    REAL product[2*PADE_MAX_DEGREE + 1] = {0.0};
    UINT i, j;

    p[0] = 1.0;
    for (i = 0; i < degree; i++)
        p[i + 1] = p[i]*(REAL)(degree - i)
                / ((REAL)(i + 1)*(REAL)(2*degree - i));

    for (i = 0; i <= degree; i++) {
        for (j = 0; j <= degree; j++) {
            square[i + j] += p[i]*p[j];
            product[i + j] += p[i]*p[j]
                    * ((j & 1U) == 0 ? 1.0 : -1.0);
        }
    }
    for (i = 0; i <= degree; i++) {
        f0[i] = square[2*i];
        denominator[i] = product[2*i];
        if (i < degree)
            temp[i] = square[2*i + 1];
    }
}

static void build_local_transition(const UINT D, const UINT index,
        COMPLEX const * const q, COMPLEX const * const r,
        const REAL eps_t, const UINT method_order,
        const UINT pade_degree, const REAL h,
        COMPLEX numerator[4][LOCAL_MAX_DEGREE + 1],
        COMPLEX denominator_polynomial[LOCAL_MAX_DEGREE + 1])
{
    const UINT z_degree = method_order == 4 ? 1
            : (method_order == 6 ? 3 : 5);
    const UINT common_degree = 2*z_degree*pade_degree;
    small_matrix_t z_matrix;
    COMPLEX x[4][Z_MAX_DEGREE + 1];
    COMPLEX x_full[4][LOCAL_MAX_DEGREE + 1];
    COMPLEX lambda2_a[LOCAL_MAX_DEGREE + 1];
    COMPLEX lambda2_b[LOCAL_MAX_DEGREE + 1];
    COMPLEX lambda2[LOCAL_MAX_DEGREE + 1];
    COMPLEX lambda_power[PADE_MAX_DEGREE + 1][LOCAL_MAX_DEGREE + 1];
    COMPLEX f0_polynomial[LOCAL_MAX_DEGREE + 1];
    COMPLEX temp_polynomial[LOCAL_MAX_DEGREE + 1];
    COMPLEX product[LOCAL_MAX_DEGREE + 1];
    REAL f0[PADE_MAX_DEGREE + 1];
    REAL temp[PADE_MAX_DEGREE];
    REAL denominator[PADE_MAX_DEGREE + 1];
    UINT i, j;

    build_z_polynomial(D, index, q, r, eps_t, method_order, &z_matrix);
    for (i = 0; i < 4; i++) {
        transform_z_to_w(&z_matrix.entry[i], z_degree, h, x[i]);
        polynomial_zero(x_full[i]);
        for (j = 0; j <= z_degree; j++)
            x_full[i][j] = x[i][j];
    }

    polynomial_multiply(x_full[0], z_degree, x_full[0], z_degree,
            lambda2_a);
    polynomial_multiply(x_full[1], z_degree, x_full[2], z_degree,
            lambda2_b);
    for (i = 0; i <= LOCAL_MAX_DEGREE; i++)
        lambda2[i] = lambda2_a[i] + lambda2_b[i];

    polynomial_zero(lambda_power[0]);
    lambda_power[0][0] = 1.0;
    for (i = 1; i <= pade_degree; i++)
        polynomial_multiply(lambda_power[i - 1], 2*z_degree*(i - 1),
                lambda2, 2*z_degree, lambda_power[i]);

    pade_coefficients(pade_degree, f0, temp, denominator);
    polynomial_zero(f0_polynomial);
    polynomial_zero(temp_polynomial);
    polynomial_zero(denominator_polynomial);
    for (i = 0; i <= pade_degree; i++) {
        polynomial_add_padded_power(f0_polynomial, lambda_power[i],
                2*z_degree*i, common_degree - 2*z_degree*i, f0[i]);
        polynomial_add_padded_power(denominator_polynomial, lambda_power[i],
                2*z_degree*i, common_degree - 2*z_degree*i, denominator[i]);
        if (i < pade_degree) {
            polynomial_add_padded_power(temp_polynomial, lambda_power[i],
                    2*z_degree*i,
                    common_degree - z_degree - 2*z_degree*i, temp[i]);
        }
    }

    for (i = 0; i < 4; i++) {
        polynomial_zero(numerator[i]);
        polynomial_multiply(temp_polynomial, common_degree - z_degree,
                x_full[i], z_degree, product);
        for (j = 0; j <= common_degree; j++)
            numerator[i][j] = product[j];
    }
    for (j = 0; j <= common_degree; j++) {
        numerator[0][j] += f0_polynomial[j];
        numerator[3][j] += f0_polynomial[j];
    }
}

static UINT local_degree(const UINT method_order, const UINT pade_degree)
{
    if (method_order == 4 && pade_degree >= 2 && pade_degree <= 7)
        return 2*pade_degree;
    if (method_order == 6 && pade_degree >= 3 && pade_degree <= 7)
        return 6*pade_degree;
    if (method_order == 8 && pade_degree >= 3 && pade_degree <= 7)
        return 10*pade_degree;
    return 0;
}

static void reverse_polynomial(COMPLEX * const p, const UINT degree)
{
    UINT i;

    for (i = 0; i < (degree + 1)/2; i++) {
        const COMPLEX value = p[i];
        p[i] = p[degree - i];
        p[degree - i] = value;
    }
}

static UINT checked_numel(const UINT D, const UINT degree,
        const UINT matrix_size)
{
    UINT padded_D;
    const UINT coefficient_count = matrix_size*(degree + 1);

    if (D == 0 || D > INT_MAX)
        return 0;
    padded_D = misc_nextpowerof2(D);
    if (padded_D == 0 || coefficient_count == 0
            || padded_D > UINT_MAX/coefficient_count)
        return 0;
    return coefficient_count*padded_D;
}

UINT akns_fscatter_pade_numel(const UINT D, const UINT method_order,
        const UINT pade_degree)
{
    const UINT degree = local_degree(method_order, pade_degree);

    if (D == 0 || degree == 0)
        return 0;
    return checked_numel(D, degree, 4);
}

UINT akns_fscatter_pade_den_numel(const UINT D, const UINT method_order,
        const UINT pade_degree)
{
    const UINT degree = local_degree(method_order, pade_degree);

    if (D == 0 || degree == 0)
        return 0;
    return checked_numel(D, degree, 1);
}

INT akns_fscatter_pade(const UINT D, COMPLEX const * const q,
        COMPLEX const * const r, const REAL eps_t, const UINT method_order,
        const UINT pade_degree, const REAL h, COMPLEX * const numerator,
        UINT * const numerator_degree, INT * const numerator_exponent,
        COMPLEX * const denominator, UINT * const denominator_degree,
        INT * const denominator_exponent)
{
    const UINT degree = local_degree(method_order, pade_degree);
    UINT matrix_numel, denominator_numel, i, j, component;
    COMPLEX *local_matrices = NULL;
    COMPLEX local_numerator[4][LOCAL_MAX_DEGREE + 1];
    COMPLEX local_denominator[LOCAL_MAX_DEGREE + 1];
    INT ret_code;

    if (D < (method_order == 8 ? 7U : 5U) || D > INT_MAX)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (!(eps_t > 0.0))
        return E_INVALID_ARGUMENT(eps_t);
    if (degree == 0)
        return E_INVALID_ARGUMENT(pade_degree);
    if (!(h > 0.0) || h == INFINITY)
        return E_INVALID_ARGUMENT(h);
    if (numerator == NULL || numerator_degree == NULL
            || denominator == NULL || denominator_degree == NULL)
        return E_INVALID_ARGUMENT(numerator);

    matrix_numel = akns_fscatter_pade_numel(D, method_order, pade_degree);
    denominator_numel = akns_fscatter_pade_den_numel(D, method_order,
            pade_degree);
    if (matrix_numel == 0 || denominator_numel == 0)
        return E_INVALID_ARGUMENT(D);
    local_matrices = malloc(matrix_numel*sizeof(COMPLEX));
    if (local_matrices == NULL)
        return E_NOMEM;

    for (i = 0; i < matrix_numel; i++)
        local_matrices[i] = 0.0;
    for (i = 0; i < denominator_numel; i++)
        denominator[i] = 0.0;

    for (i = 0; i < D; i++) {
        build_local_transition(D, i, q, r, eps_t, method_order, pade_degree,
                h, local_numerator, local_denominator);
        for (component = 0; component < 4; component++) {
            COMPLEX * const destination = local_matrices
                    + component*D*(degree + 1)
                    + (D - 1 - i)*(degree + 1);
            for (j = 0; j <= degree; j++)
                destination[j] = local_numerator[component][j];
        }
        for (j = 0; j <= degree; j++)
            denominator[i*(degree + 1) + j] = local_denominator[j];
    }

    *numerator_degree = degree;
    ret_code = poly_fmult2x2(numerator_degree, D, local_matrices, numerator,
            numerator_exponent);
    if (ret_code != SUCCESS)
        goto release_mem;
    for (component = 0; component < 4; component++)
        reverse_polynomial(numerator + component*(*numerator_degree + 1),
                *numerator_degree);

    *denominator_degree = degree;
    ret_code = poly_fmult(denominator_degree, D, denominator,
            denominator_exponent);
    if (ret_code == SUCCESS)
        reverse_polynomial(denominator, *denominator_degree);

release_mem:
    free(local_matrices);
    return ret_code;
}
