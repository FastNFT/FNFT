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

#include "fnft_nsev_inverse.h"
#include "fnft_nsev.h"
#include "fnft__nse_fscatter.h"
#include "fnft__poly_chirpz.h"
#include "fnft__errwarn.h"
#include "fnft__nse_discretization.h"
#include "fnft__misc.h"
#include <stdio.h>

static INT fnft_nsev_inverse_test(const UINT M, const REAL error_bound,
    const INT kappa, const nse_discretization_t discretization)
{
    const UINT D = 8;
    COMPLEX q_exact[8] = { 0.1, 0.1*I, 0.2, -0.2, 0.0, 0.05+0.05*I, -0.03*I,
        0.06 };
    COMPLEX q[8];
    REAL T[2] = { 0.0, D-1.0 };
    REAL XI[2];

    if (M < D)
        return E_INVALID_ARGUMENT(M);
    if (!(error_bound > 0))
        return E_INVALID_ARGUMENT(error_bound);

    INT ret_code = SUCCESS;
    UINT i;

    const UINT numel = nse_fscatter_numel(D, discretization);
    COMPLEX * const transfer_matrix = malloc(numel * sizeof(COMPLEX));
    COMPLEX * const a_vals = malloc(M * sizeof(COMPLEX));
    COMPLEX * const b_vals = malloc(M * sizeof(COMPLEX));
    COMPLEX * const contspec = malloc(M * sizeof(COMPLEX));
    if (transfer_matrix == NULL || a_vals == NULL || b_vals == NULL
    || contspec == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    fnft_nsev_opts_t opts_fwd = fnft_nsev_default_opts();
    opts_fwd.discretization = discretization;
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft_nsev(D, q_exact, T, M, contspec, XI, NULL, NULL, NULL, kappa,
        &opts_fwd);

    fnft_nsev_inverse_opts_t opts_inv = fnft_nsev_inverse_default_opts();
    opts_inv.discretization = discretization;
    ret_code = fnft_nsev_inverse(M, contspec, XI, 0, NULL, NULL, D, q, T,
        kappa, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    UINT deg;
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    ret_code = nse_fscatter(D, q_exact, eps_t, kappa, transfer_matrix, &deg,
        NULL, discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    const COMPLEX A = 1.0;
    const COMPLEX W = CEXP(2.0*I*FNFT_PI/M);
    ret_code = poly_chirpz(deg, transfer_matrix, A, W, M, a_vals);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = poly_chirpz(deg, transfer_matrix+2*(deg+1), A, W, M, b_vals);
    CHECK_RETCODE(ret_code, leave_fun);

    for (i=0; i<M; i++)
        contspec[i] = b_vals[i]/a_vals[i];

    ret_code = fnft_nsev_inverse(M, contspec, NULL, 0, NULL, NULL, D, q, T,
        kappa, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);

    error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(transfer_matrix);
    free(a_vals);
    free(b_vals);
    free(contspec);
    return ret_code;
}

