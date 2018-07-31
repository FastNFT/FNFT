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

static INT nsev_inverse_test_addsoliton_cdt(UINT const D, REAL const error_bound)
{
    const UINT M = D+2;
    const UINT K = 3;
    COMPLEX * q_exact = NULL;
    COMPLEX * q = NULL;
    COMPLEX bound_states[3] = {2.5 + 0.9*I, 2.5 + 1.9*I, 2.5 + 2.9*I};
    COMPLEX normconsts_or_residues[3] = {-1.0, 1.0, -1.0};
    REAL T[2] = {-20.0, 20.0};
    REAL XI[2];
    REAL t;
    UINT i, j;
    INT ret_code;
    q_exact = malloc(D * sizeof(COMPLEX));
    q = malloc(D * sizeof(COMPLEX));
    
    if (q_exact == NULL || q == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < D; i++){
        t = T[0] + i*(T[1] - T[0])/(D - 1);
        q_exact[i] = 3.4*misc_sech(t)*CEXP(-5*I*t);
        q[i] = 0.4*misc_sech(t)*CEXP(-5*I*t);
    }
    
    //Testing with norming constants
    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();
    opts.discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANT;
    opts.discspec_inversion_method
            = fnft_nsev_inverse_dsmethod_ADDSOLITON_CDT;
    
    ret_code = fnft_nsev_inverse(M, NULL, XI, K, bound_states, normconsts_or_residues, D, q, T,
            1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    REAL error = misc_rel_err(D, q, q_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
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
    
    for (i = 0; i < D; i++){
        t = T[0] + i*(T[1] - T[0])/(D - 1);
        q[i] = 0.4*misc_sech(t)*CEXP(-5*I*t);
    }
    opts.discspec_type = fnft_nsev_inverse_dstype_RESIDUE;
    opts.discspec_inversion_method
            = fnft_nsev_inverse_dsmethod_ADDSOLITON_CDT;
    
    ret_code = fnft_nsev_inverse(M, NULL, XI, K, bound_states, normconsts_or_residues, D, q, T,
            1, &opts);
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
        free(q_exact);
        free(q);
        return ret_code;
}

INT main()
{

    UINT D = 512;
    REAL error_bound = 0.0029;
    if (nsev_inverse_test_addsoliton_cdt(D, error_bound) != SUCCESS)
        return EXIT_FAILURE;
   //Checking for quadratic error decay
    if (nsev_inverse_test_addsoliton_cdt(D*2, error_bound/4) != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
