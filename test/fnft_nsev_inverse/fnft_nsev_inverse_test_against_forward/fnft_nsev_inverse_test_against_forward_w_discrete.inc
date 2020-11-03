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
   const REAL error_bound, fnft_nsev_inverse_opts_t * opts_inv_ptr)
{
    if (M < D)
        return E_INVALID_ARGUMENT(M);
    if (!(error_bound > 0))
        return E_INVALID_ARGUMENT(error_bound);

    COMPLEX * q_exact = NULL;
    COMPLEX * q = NULL;
    COMPLEX * contspec = NULL;
    UINT * K_ptr;
    UINT K = 10;
    K_ptr = &K;
    COMPLEX * bound_states = NULL;
    COMPLEX * normconsts_or_residues = NULL;
    REAL T[2] = {-32.0, 32.0};
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
    //Compute q_exact
    for (i = 0; i < D; i++){
        t = T[0] + i*(T[1] - T[0])/(D - 1);
        q_exact[i] = 3.4*misc_sech(t)*CEXP(-4*I*t);
    }

    //Setting up forward transform options
    fnft_nsev_opts_t opts_fwd = fnft_nsev_default_opts();
    opts_fwd.discretization = nse_discretization_2SPLIT4B;
    if (opts_inv_ptr->discspec_type == fnft_nsev_inverse_dstype_RESIDUES)
       opts_fwd.discspec_type = nsev_dstype_RESIDUES;
    else if (opts_inv_ptr->discspec_type == fnft_nsev_inverse_dstype_NORMING_CONSTANTS)
       opts_fwd.discspec_type = nsev_dstype_NORMING_CONSTANTS;

    //Obtain XI grid
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts_inv_ptr->discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    //Perform full forward NFT
    ret_code = fnft_nsev(D, q_exact, T, M, contspec, XI, K_ptr, bound_states, normconsts_or_residues,
            kappa, &opts_fwd);
    CHECK_RETCODE(ret_code, leave_fun);

    //Perform full inverse NFT
    ret_code = fnft_nsev_inverse(M, contspec, XI, *K_ptr, bound_states, normconsts_or_residues, D, q, T,
            kappa, opts_inv_ptr);
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
