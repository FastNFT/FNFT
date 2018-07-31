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

int main()
{
    const UINT D = 256;
    const UINT M = 4*D;
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
        q_exact[i] = 3.4*misc_sech(t)*CEXP(-6*I*t);
    }
    
    //Setting up inverse transform options
    fnft_nsev_inverse_opts_t opts_inv = fnft_nsev_inverse_default_opts();
    opts_inv.discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANT;
    opts_inv.discspec_inversion_method
            = fnft_nsev_inverse_dsmethod_DEFAULT;
    opts_inv.discretization = nse_discretization_2SPLIT2A;
    opts_inv.contspec_inversion_method
            = fnft_nsev_inverse_csmethod_DEFAULT;
    
    
    //Setting up forward transform options
    fnft_nsev_opts_t opts_fwd = fnft_nsev_default_opts();
//     opts_fwd.discretization = nse_discretization_2SPLIT2A;
    
    //Obtain XI grid
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts_inv.discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    
    //Perform full forward NFT
    ret_code = fnft_nsev(D, q_exact, T, M, contspec, XI, K_ptr, bound_states, normconsts_or_residues,
            kappa, &opts_fwd);
    CHECK_RETCODE(ret_code, leave_fun);
    misc_print_buf(*K_ptr, bound_states, "EV");
    misc_print_buf(*K_ptr, normconsts_or_residues, "b");
    
    bound_states[0] = 3+.9*I;
    bound_states[1] = 3+1.9*I;
    bound_states[2] = 3+2.9*I;
    normconsts_or_residues[0] = -1;
    normconsts_or_residues[1] = 1;
    normconsts_or_residues[2] = -1;
    //Perform full inverse NFT
    ret_code = fnft_nsev_inverse(M, contspec, XI, *K_ptr, bound_states, normconsts_or_residues, D, q, T,
            kappa, &opts_inv);
//         ret_code = fnft_nsev_inverse(M, contspec, XI, 0, NULL, NULL, D, q, T,
//             kappa, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);
    misc_print_buf(D, q, "q");

    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < 1e-2)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    
    leave_fun:
        free(q_exact);
        free(q);
        free(contspec);
        free(bound_states);
        free(normconsts_or_residues);
        if (ret_code == SUCCESS)
            return EXIT_SUCCESS;
        else
            return EXIT_FAILURE;
}
