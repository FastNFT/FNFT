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
    const UINT D = 16384;
    const UINT M = D+2;
    const UINT K = 5;
    COMPLEX * q_exact = NULL;
    COMPLEX * q = NULL;
    COMPLEX bound_states[5] = {0.5*I, 1.5*I, 2.5*I, 3.5*I, 4.5*I};
    COMPLEX normconsts_or_residues[5] = {-1.0, 1.0, -1.0, 1.0, -1.0};
    REAL T[2] = {-15.0, 15.0};
    REAL XI[2];
    UINT i, j;
    INT ret_code;
    q_exact = malloc(D * sizeof(COMPLEX));
    q = malloc(D * sizeof(COMPLEX));
    
    if (q_exact == NULL || q == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < D; i++)
        q_exact[i] = 5.0*misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
    
    //Testing with norming constants
    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();
    opts.discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANT;
    opts.discspec_inversion_method
            = fnft_nsev_inverse_dsmethod_CDT;
    
    ret_code = fnft_nsev_inverse(M, NULL, XI, K, bound_states, normconsts_or_residues, D, q, T,
            1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < 100*EPSILON)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    
    //Testing with residues
    //Converting norming constants to residues
    COMPLEX tmp = 0.0;
    for (i = 0; i < K; i++){
        tmp = 1.0;
        for (j = 0; j < K; j++){
            if (j != i){
                tmp = tmp*(bound_states[i]-bound_states[j])/(bound_states[i]-CONJ(bound_states[j]));
            }
        }
        normconsts_or_residues[i] = (normconsts_or_residues[i]*(2*I*CIMAG(bound_states[i])))/tmp;
    }

    opts.discspec_type = fnft_nsev_inverse_dstype_RESIDUE;
    opts.discspec_inversion_method
            = fnft_nsev_inverse_dsmethod_CDT;
    
    ret_code = fnft_nsev_inverse(M, NULL, XI, 5, bound_states, normconsts_or_residues, D, q, T,
            1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < 100*EPSILON)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    
    leave_fun:
        free(q_exact);
        free(q);
        if (ret_code == SUCCESS)
            return EXIT_SUCCESS;
        else
            return EXIT_FAILURE;
}
