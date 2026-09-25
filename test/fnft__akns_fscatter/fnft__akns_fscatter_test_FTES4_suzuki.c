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

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft__akns_discretization.h"
#include "fnft__akns_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

static void matrix_mult(COMPLEX const A[4], COMPLEX const B[4], COMPLEX C[4])
{
    COMPLEX const tmp[4] = {
        A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3],
        A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3]
    };
    UINT i;

    for (i = 0; i < 4; i++)
        C[i] = tmp[i];
}

static void zero_freq_exponential(COMPLEX M[4], const REAL h,
                                  const COMPLEX q, const COMPLEX r)
{
    const COMPLEX delta = h*CSQRT(-q*r);
    const COMPLEX offdiag = h*misc_CSINC(delta);

    M[0] = CCOS(delta);
    M[1] = q*offdiag;
    M[2] = r*offdiag;
    M[3] = M[0];
}

static void matrix_mult_assign(COMPLEX M[4], COMPLEX const R[4])
{
    COMPLEX tmp[4];
    UINT i;

    matrix_mult(M, R, tmp);
    for (i = 0; i < 4; i++)
        M[i] = tmp[i];
}

static void direct_step(const UINT D, const UINT i,
                        COMPLEX const * const q, COMPLEX const * const r,
                        const REAL eps_t, const COMPLEX z, COMPLEX U[4])
{
    const UINT ip = (i + 1) % D;
    const UINT im = (i + D - 1) % D;
    const COMPLEX dq = (q[ip] - q[im])/24.0;
    const COMPLEX d2q = (q[im] - 2.0*q[i] + q[ip])/48.0;
    const COMPLEX dr = (r[ip] - r[im])/24.0;
    const COMPLEX d2r = (r[im] - 2.0*r[i] + r[ip])/48.0;
    const COMPLEX A_plus[4] = {1.0, 0.0, 0.0, z};
    const COMPLEX A_minus[4] = {z, 0.0, 0.0, 1.0};
    const COMPLEX A_full[4] = {1.0, 0.0, 0.0, z*z*z};
    COMPLEX e_plus[4], e_minus[4], e_7_48B[4], e_3_8B[4], e_m1_48B[4];
    COMPLEX tmp[4];

    zero_freq_exponential(e_plus, eps_t, d2q + dq, d2r + dr);
    zero_freq_exponential(e_minus, eps_t, d2q - dq, d2r - dr);
    zero_freq_exponential(e_7_48B, 7.0*eps_t/48.0, q[i], r[i]);
    zero_freq_exponential(e_3_8B, 3.0*eps_t/8.0, q[i], r[i]);
    zero_freq_exponential(e_m1_48B, -eps_t/48.0, q[i], r[i]);

    U[0] = 1.0;
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = 1.0;
    matrix_mult_assign(U, e_7_48B);
    matrix_mult_assign(U, A_plus);
    matrix_mult_assign(U, e_3_8B);
    matrix_mult_assign(U, A_minus);
    matrix_mult_assign(U, e_m1_48B);
    matrix_mult_assign(U, A_full);
    matrix_mult_assign(U, e_m1_48B);
    matrix_mult_assign(U, A_minus);
    matrix_mult_assign(U, e_3_8B);
    matrix_mult_assign(U, A_plus);
    matrix_mult_assign(U, e_7_48B);

    matrix_mult(e_plus, U, tmp);
    matrix_mult(tmp, e_minus, U);
}

static void direct_scatter(const UINT D, COMPLEX const * const q,
                           COMPLEX const * const r, const REAL eps_t,
                           const COMPLEX z, COMPLEX S[4])
{
    COMPLEX U[4], tmp[4];
    UINT i, j;

    S[0] = 1.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 1.0;
    for (i = 0; i < D; i++) {
        direct_step(D, i, q, r, eps_t, z, U);
        matrix_mult(U, S, tmp);
        for (j = 0; j < 4; j++)
            S[j] = tmp[j];
    }
}

static INT check_invariant(const INT kappa, COMPLEX const S[4])
{
    const COMPLEX v[2] = {0.73 + 0.17*I, -0.29 + 0.41*I};
    const COMPLEX Sv[2] = {
        S[0]*v[0] + S[1]*v[1], S[2]*v[0] + S[3]*v[1]
    };
    const REAL before = CREAL(v[0]*CONJ(v[0])) +
            kappa*CREAL(v[1]*CONJ(v[1]));
    const REAL after = CREAL(Sv[0]*CONJ(Sv[0])) +
            kappa*CREAL(Sv[1]*CONJ(Sv[1]));

    return FABS(after - before) <= 20000*EPSILON ? SUCCESS : E_TEST_FAILED;
}

static INT run_case(const INT kappa, const UINT normalize)
{
    const UINT D = 7, nz = 7;
    const REAL eps_t = 0.09;
    const REAL err_bound = 20000*EPSILON;
    const COMPLEX z[7] = {1.0, CEXP(I*PI/7), CEXP(-I*PI/5),
                          CEXP(I*8*PI/11), CEXP(-I*13*PI/17),
                          0.83 + 0.17*I, 1.09 - 0.11*I};
    COMPLEX q[7], r[7], reference[28], result[28], S[4];
    COMPLEX *transfer_matrix = NULL;
    UINT i, j, deg, numel;
    INT W = 0;
    INT ret_code = SUCCESS;

    for (i = 0; i < D; i++) {
        q[i] = 0.31*COS(0.7*(i + 1)) + 0.23*I*SIN(0.4*(i + 1));
        r[i] = -kappa*CONJ(q[i]);
    }

    numel = akns_fscatter_numel(D, akns_discretization_FTES4_suzuki);
    if (numel == 0)
        return E_INVALID_ARGUMENT(akns_discretization_FTES4_suzuki);
    transfer_matrix = malloc(numel*sizeof(COMPLEX));
    if (transfer_matrix == NULL)
        return E_NOMEM;

    ret_code = akns_fscatter(D, q, r, eps_t, transfer_matrix, &deg,
                             normalize ? &W : NULL,
                             akns_discretization_FTES4_suzuki);
    CHECK_RETCODE(ret_code, leave_fun);
    if (deg != 7*D) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (normalize) {
        const REAL scale = POW(2.0, W);
        for (i = 0; i < 4*(deg + 1); i++)
            transfer_matrix[i] *= scale;
    }

    for (i = 0; i < 4; i++) {
        for (j = 0; j < nz; j++)
            result[i*nz + j] = z[j];
        ret_code = poly_eval(deg, transfer_matrix + i*(deg + 1), nz,
                             result + i*nz);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    for (j = 0; j < nz; j++) {
        direct_scatter(D, q, r, eps_t, z[j], S);
        for (i = 0; i < 4; i++)
            reference[i*nz + j] = S[i];
    }

    if (misc_rel_err(4*nz, result, reference) > err_bound) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    for (i = 0; i < 4; i++)
        S[i] = result[i*nz + 1];
    ret_code = check_invariant(kappa, S);

leave_fun:
    free(transfer_matrix);
    return ret_code;
}

static INT check_mapping(void)
{
    const REAL eps_t = 0.09;
    COMPLEX vals[3] = {-1.3 + 0.2*I, 0.37 - 0.11*I, 1.71 + 0.05*I};
    COMPLEX expected_z[3];
    const COMPLEX reference[3] = {
        -1.3 + 0.2*I, 0.37 - 0.11*I, 1.71 + 0.05*I
    };
    UINT i;
    INT ret_code;

    for (i = 0; i < 3; i++)
        expected_z[i] = CEXP(2.0*I*reference[i]*eps_t/3.0);
    ret_code = akns_discretization_lambda_to_z(3, eps_t, vals,
            akns_discretization_FTES4_suzuki);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(3, vals, expected_z) > 100*EPSILON) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    ret_code = akns_discretization_z_to_lambda(3, eps_t, vals,
            akns_discretization_FTES4_suzuki);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(3, vals, reference) > 100*EPSILON)
        ret_code = E_TEST_FAILED;

leave_fun:
    return ret_code;
}

static INT run_nsev_full_filtering(void)
{
    const UINT D = 128;
    const REAL T[2] = {-25.0, 25.0};
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    COMPLEX *q = NULL, *bound_states = NULL;
    UINT i, K, K_basic;
    REAL const eps_t = (T[1] - T[0])/(D - 1);
    REAL const re_limit = 0.9*PI/((2.0/3.0)*eps_t);
    INT ret_code = SUCCESS;

    opts.discretization = nse_discretization_FTES4_suzuki;
    opts.bound_state_localization = nsev_bsloc_FAST_EIGENVALUE;
    q = malloc(D*sizeof(COMPLEX));
    K = fnft_nsev_max_K(D, &opts);
    bound_states = malloc(K*sizeof(COMPLEX));
    if (q == NULL || bound_states == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < D; i++)
        q[i] = 3.2*I*misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));

    opts.bound_state_filtering = nsev_bsfilt_BASIC;
    ret_code = fnft_nsev(D, q, T, 0, NULL, NULL, &K, bound_states, NULL,
            +1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    K_basic = K;

    K = fnft_nsev_max_K(D, &opts);
    opts.bound_state_filtering = nsev_bsfilt_FULL;
    ret_code = fnft_nsev(D, q, T, 0, NULL, NULL, &K, bound_states, NULL,
            +1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    if (K >= K_basic) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    for (i = 0; i < K; i++) {
        if (FABS(CREAL(bound_states[i])) > re_limit*(1.0 + 10*EPSILON)) {
            ret_code = E_TEST_FAILED;
            goto leave_fun;
        }
    }

leave_fun:
    free(q);
    free(bound_states);
    return ret_code;
}

static INT run_nsev_convergence(void)
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 512;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
    const nsev_testcases_t tc_shifted = nsev_testcases_SECH_FOCUSING2;
    REAL error_bounds[6] = {
        3.0e-6, // reflection coefficient
        4.0e-6, // a
        3.0e-6, // b
        6.0e-6, // bound states
        5.0e-15,// norming constants
        1.3e-5  // residues
    };
    REAL shifted_error_bounds[6] = {
        8.0e-3, // reflection coefficient
        4.0e-3, // a
        2.0e-3, // b
        2.0e-3, // bound states
        3.0e-14,// norming constants
        3.0e-3  // residues
    };
    REAL richardson_bounds[6] = {
        8.0e-8, // reflection coefficient
        4.0e-8, // a
        8.0e-8, // b
        1.2e-7, // bound states
        5.0e-15,// norming constants
        2.2e-7  // residues
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_FTES4_suzuki;

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    D *= 2;
    for (i = 0; i < 6; i++)
        error_bounds[i] /= 16.0;
    error_bounds[4] *= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    D = 1024;
    ret_code = nsev_testcases_test_fnft(tc_shifted, D, shifted_error_bounds,
            &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    D = 512;
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, richardson_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}

INT main()
{
    INT kappa;

    if (akns_discretization_degree(akns_discretization_FTES4_suzuki) != 7 ||
        akns_discretization_method_order(akns_discretization_FTES4_suzuki) != 4)
        return EXIT_FAILURE;
    if (check_mapping() != SUCCESS)
        return EXIT_FAILURE;
    for (kappa = -1; kappa <= 1; kappa += 2) {
        if (run_case(kappa, 0) != SUCCESS || run_case(kappa, 1) != SUCCESS)
            return EXIT_FAILURE;
    }
    if (run_nsev_full_filtering() != SUCCESS)
        return EXIT_FAILURE;
    if (run_nsev_convergence() != SUCCESS)
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
