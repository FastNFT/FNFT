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

#include <stdio.h>
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
    const COMPLEX Z[4] = {1.0, 0.0, 0.0, z};
    const COMPLEX Z2[4] = {1.0, 0.0, 0.0, z*z};
    COMPLEX e_plus[4], e_minus[4], e_half[4], e_full[4];
    COMPLEX tmp1[4], tmp2[4], term1[4], term2[4], base[4];
    UINT j;

    zero_freq_exponential(e_plus, eps_t, d2q + dq, d2r + dr);
    zero_freq_exponential(e_minus, eps_t, d2q - dq, d2r - dr);
    zero_freq_exponential(e_half, eps_t/2.0, q[i], r[i]);
    zero_freq_exponential(e_full, eps_t, q[i], r[i]);

    matrix_mult(Z, e_half, tmp1);
    matrix_mult(tmp1, Z2, tmp2);
    matrix_mult(tmp2, e_half, tmp1);
    matrix_mult(tmp1, Z, term1);

    matrix_mult(Z2, e_full, tmp1);
    matrix_mult(tmp1, Z2, term2);

    for (j = 0; j < 4; j++)
        base[j] = (4.0*term1[j] - term2[j])/3.0;

    matrix_mult(e_plus, base, tmp1);
    matrix_mult(tmp1, e_minus, U);
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

static INT run_case(const INT kappa, const UINT normalize)
{
    const UINT D = 7, nz = 5;
    const REAL eps_t = 0.09;
    const REAL err_bound = 5000*EPSILON;
    const COMPLEX z[5] = {1.0, CEXP(I*PI/7), CEXP(-I*PI/5),
                          CEXP(I*8*PI/11), CEXP(-I*13*PI/17)};
    COMPLEX q[7], r[7], reference[20], result[20], S[4];
    COMPLEX *transfer_matrix = NULL;
    UINT i, j, deg, numel;
    INT W = 0;
    INT ret_code = SUCCESS;

    for (i = 0; i < D; i++) {
        q[i] = 0.31*COS(0.7*(i + 1)) + 0.23*I*SIN(0.4*(i + 1));
        r[i] = -kappa*CONJ(q[i]);
    }

    numel = akns_fscatter_numel(D, akns_discretization_FTES4_4A);
    if (numel == 0)
        return E_INVALID_ARGUMENT(akns_discretization_FTES4_4A);
    transfer_matrix = malloc(numel*sizeof(COMPLEX));
    if (transfer_matrix == NULL)
        return E_NOMEM;

    ret_code = akns_fscatter(D, q, r, eps_t, transfer_matrix, &deg,
                             normalize ? &W : NULL,
                             akns_discretization_FTES4_4A);
    CHECK_RETCODE(ret_code, leave_fun);
    if (deg != 4*D) {
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

#ifdef DEBUG
    printf("kappa=%d normalize=%u error=%g bound=%g\n", kappa,
           (unsigned)normalize,
           misc_rel_err(4*nz, result, reference), err_bound);
#endif
    if (misc_rel_err(4*nz, result, reference) > err_bound)
        ret_code = E_TEST_FAILED;

leave_fun:
    free(transfer_matrix);
    return ret_code;
}

static INT run_nsev_convergence(void)
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 512;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
    REAL error_bounds[6] = {
        3.0e-5, // reflection coefficient
        8.0e-5, // a
        3.0e-5, // b
        6.0e-6, // bound states
        5e-15,  // norming constants
        1.3e-5  // residues
    };
    REAL richardson_bounds[6] = {
        5.5e-8, // reflection coefficient
        7.0e-7, // a
        1.2e-7, // b
        1.8e-9, // bound states
        5e-15,  // norming constants
        3.2e-9  // residues
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_FTES4_4A;

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    D *= 2;
    for (i = 0; i < 6; i++)
        error_bounds[i] /= 16.0;
    error_bounds[4] *= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, richardson_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}

INT main()
{
    INT kappa;

    if (akns_discretization_degree(akns_discretization_FTES4_4A) != 4 ||
        akns_discretization_method_order(akns_discretization_FTES4_4A) != 4)
        return EXIT_FAILURE;
    for (kappa = -1; kappa <= 1; kappa += 2) {
        if (run_case(kappa, 0) != SUCCESS || run_case(kappa, 1) != SUCCESS)
            return EXIT_FAILURE;
    }
    if (run_nsev_convergence() != SUCCESS)
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
