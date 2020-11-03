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
 * Shrinivas Chimmalgi (TU Delft) 2018.
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
   const REAL error_bound, fnft_nsev_inverse_opts_t * opts_ptr)
{
    if (M < D)
        return E_INVALID_ARGUMENT(M);
    if (!(error_bound > 0))
        return E_INVALID_ARGUMENT(error_bound);

    COMPLEX * q_exact = NULL;
    COMPLEX * q = NULL;
    COMPLEX * contspec = NULL;
    UINT * K_ptr;
    UINT K = 1;
    K_ptr = &K;
    COMPLEX * bound_states = NULL;
    COMPLEX * normconsts_or_residues = NULL;
    REAL T[2] = {-2.0, 2.0};
    REAL XI[2], t = 0;
    UINT i;
    INT const kappa = 1;
    INT ret_code = SUCCESS;

    q_exact = malloc(D * sizeof(COMPLEX));
    q = malloc(D * sizeof(COMPLEX));
    contspec = malloc(M * sizeof(COMPLEX));
    bound_states = malloc(*K_ptr * sizeof(COMPLEX));
    normconsts_or_residues = malloc(*K_ptr * sizeof(COMPLEX));
    if (q_exact == NULL || q == NULL || contspec == NULL
            || bound_states == NULL || normconsts_or_residues == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    REAL al = 2;
    REAL be = 0.55;
    REAL gam = (REAL) CSQRT(FABS(al)*FABS(al) + be*be);

    //Obtain XI grid
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    
    //Compute contspec
    const REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    for (i = 0; i < M; i++)
        contspec[i] = al/((XI[0] + eps_xi*i) - be*I);

    //Discrete spectrum
    bound_states[0] = be*I;
    normconsts_or_residues[0] = -I*al/(gam+be);

    //Compute q_exact
    for (i = 0; i < D; i++){
        t = T[0] + i*(T[1] - T[0])/(D - 1);
        if (t <= 0)
           q_exact[i] = -(2.0*I*gam*al/FABS(al))*misc_sech(2*gam*t + ATANH(be/gam));
        else
           q_exact[i] = 0.0;
    }


    //Perform full inverse NFT
    ret_code = fnft_nsev_inverse(M, contspec, XI, *K_ptr, bound_states, normconsts_or_residues, D, q, T,
            kappa, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    leave_fun:
        free(q_exact);
        free(q);
        free(contspec);
        free(bound_states);
        free(normconsts_or_residues);
    return ret_code;
}
