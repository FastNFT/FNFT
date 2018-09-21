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

static INT fnft_nsev_inverse_test(const UINT D, const UINT M,
                                  const REAL error_bound,
                                  fnft_nsev_inverse_opts_t * opts_ptr)
{
    REAL XI[2] = { -1, 1 };

    if (!(error_bound > 0))
        return E_INVALID_ARGUMENT(error_bound);
    if (opts_ptr == NULL)
        return E_INVALID_ARGUMENT(opts_ptr);

    INT ret_code = SUCCESS;
    UINT i;

    COMPLEX * const q = malloc(D * sizeof(COMPLEX));
    COMPLEX * const q_exact = malloc(D * sizeof(COMPLEX));
    COMPLEX * const contspec = malloc(M * sizeof(COMPLEX));
    COMPLEX * bound_states = NULL;
    COMPLEX * normconsts = NULL;
    if (q == NULL || q_exact == NULL || contspec == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    REAL T[2] = { -25 , 25 };
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    const REAL A = 3.45; // changing this can change number of bound_states
    const REAL t0 = 0; //Difficult to know normconsts if this is nonzero
    for (i=0; i<D; i++) {
        const REAL t = T[0] + i*eps_t;
        q_exact[i] = I*A*misc_sech(t - t0);
    }
    UINT const K = (UINT) floor(A+0.5);
    bound_states = malloc(K * sizeof(COMPLEX));
    normconsts = malloc(K * sizeof(COMPLEX));
    if (bound_states == NULL || normconsts == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    for (i=0; i<K; i++) {
        bound_states[i] = I*((A+0.5)-(i+1));
        normconsts[i] = -I*POW(-1,i+1);
    }

    if (opts_ptr->contspec_type == fnft_nsev_inverse_cstype_B_OF_TAU) {

        for (i=0; i<D; i++) {
            const REAL t = T[0] + i*eps_t;
            const REAL tau = 2*t;
            const REAL tau0 = 2*t0;
            contspec[i] = I/(2*FNFT_PI)*sin(FNFT_PI*A)*misc_sech((tau-tau0)/2);
        }

    } else if (opts_ptr->contspec_type == fnft_nsev_inverse_cstype_B_OF_XI) {

        ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        const REAL eps_xi = (XI[1] - XI[0])/(M - 1);
        for (i=0; i<M; i++) {
            const REAL xi = XI[0] + i*eps_xi;
            contspec[i] = I*CEXP(-2*I*xi*t0)*SIN(FNFT_PI*A)/COSH(FNFT_PI*xi);
        }

    } else {

        ret_code = E_INVALID_ARGUMENT(opts_ptr->contspec_type);
        goto leave_fun;

    }

    ret_code = fnft_nsev_inverse(M, contspec, XI, K, bound_states, normconsts, D, q, T, +1,
                                 opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        printf("error = %g\n", error);
        goto leave_fun;
    }

leave_fun:
    free(q);
    free(q_exact);
    free(contspec);
    free(bound_states);
    free(normconsts);
    return ret_code;
}
