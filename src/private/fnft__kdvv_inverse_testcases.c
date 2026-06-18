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
* Fabian Fischer (Hiwi KIT) 2026
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft__kdvv_inverse_testcases.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

void kdvv_print_spectrum(   COMPLEX const * const bound_states,
                            COMPLEX const * const normconsts,
                            COMPLEX const * const contspec,
                            REAL const * const XI,
                            const UINT M,
                            const UINT D,
                            const UINT K)
{
    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    FNFT_REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    printf("Continuous spectrum:\n");
    for (FNFT_UINT i=0; i<M; i++) {
        FNFT_REAL xi = XI[0] + i*eps_xi;
        printf("  continuous_spectrum(xi=%f) \t= %g + %gI\n",
            (double)xi,
            (double)FNFT_CREAL(contspec[i]),
            (double)FNFT_CIMAG(contspec[i])
        );
    }

    printf("Discrete spectrum:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  bound state at %g + %gI with norming constant %g + %gI\n",
            (double)FNFT_CREAL(bound_states[i]),
            (double)FNFT_CIMAG(bound_states[i]),
            (double)FNFT_CREAL(normconsts[i]),
            (double)FNFT_CIMAG(normconsts[i])
        );
    }
}



INT kdvv_testcases_get_spectrum_of_inverse(const fnft_kdvv_params params_i,
                    const REAL err_bnd_bound_states,
                    const REAL err_bnd_spurious_bound_states,
                    const REAL err_bnd_normconst,
                    const REAL err_bnd_contspec)
{
    INT ret_code = SUCCESS;
    
    COMPLEX * q = NULL;
    COMPLEX * contspec_r = NULL;
    COMPLEX * bound_states_r = NULL;
    COMPLEX * normconsts_r = NULL;
    

    q = malloc(params_i.D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);

    ret_code = fnft_kdvv_inverse(   params_i.M, 
                                    params_i.contspec, 
                                    params_i.XI, 
                                    params_i.K, 
                                    params_i.bound_states, 
                                    params_i.normconsts, 
                                    params_i.D, 
                                    q, 
                                    params_i.T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Prepare forward fnft_kdvv
    contspec_r = malloc(params_i.M * sizeof(COMPLEX));
    CHECK_NOMEM(contspec_r, ret_code, leave_fun);

    UINT K_r = params_i.D;

    bound_states_r = malloc(K_r * sizeof(COMPLEX));
    CHECK_NOMEM(bound_states_r, ret_code, leave_fun);

    normconsts_r = malloc(K_r * sizeof(COMPLEX));
    CHECK_NOMEM(normconsts_r, ret_code, leave_fun);

    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();

    ret_code = fnft_kdvv(params_i.D, q, params_i.T, params_i.M, contspec_r, params_i.XI, &K_r, bound_states_r, normconsts_r, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    #ifdef DEBUG
        kdvv_print_spectrum(bound_states_r, normconsts_r, contspec_r, params_i.XI, params_i.M, params_i.D, K_r);
    #endif

    // -- Check results --

    // Simple general tests
    for (UINT i=0; i<params_i.K; i++){
        COMPLEX bsr = bound_states_r[i];
        COMPLEX ncr = normconsts_r[i];

        UINT is_bsr_pure_imaginary = FABS(CREAL(bsr)) < 1e-8;
        UINT is_bsr_positive_imaginary = (CIMAG(bsr) > 0) && is_bsr_pure_imaginary;

        UINT is_ncr_real = FABS(CIMAG(ncr)) < 1e-8;

        if (is_bsr_pure_imaginary &&
            is_bsr_positive_imaginary &&
            is_ncr_real) 
        {
            ret_code = SUCCESS;
        } 
        else {
            ret_code = FNFT_EC_TEST_FAILED;
            break;
        }
    }
    CHECK_RETCODE(ret_code, leave_fun);

    // Check matching of the computed  bound states
    // bound_states_r is sorted in ascending order, K_r is the number of found bound states by fnft_kdvv
    // take only the last params_i.K bound states, because they are the bigger ones and more likely the ones
    // which corresponds to the initial given bound states
    COMPLEX * const candidate_bound_states_r_ptr = &bound_states_r[K_r-params_i.K];
    REAL hausdorff_dist_bound_states = misc_hausdorff_dist_normed(  params_i.K, params_i.bound_states, params_i.K, 
                                                                    candidate_bound_states_r_ptr);
    UINT is_bsr_in_tolerance = hausdorff_dist_bound_states < err_bnd_bound_states;                                                                

    COMPLEX * const candidate_normconsts_r_ptr = &normconsts_r[K_r-params_i.K];
    REAL hausdorff_dist_normconsts = misc_hausdorff_dist_normed(params_i.K, params_i.normconsts, params_i.K, 
                                                                candidate_normconsts_r_ptr);
    UINT is_ncr_in_tolerance = hausdorff_dist_normconsts < err_bnd_normconst;                                                          

    #ifdef DEBUG
        printf("Number of bound_states:\n  K_r = %u\n", (unsigned int)K_r);
        printf("Hausdorff dist bound states:\n  dist = %f\n", hausdorff_dist_bound_states);
        printf("Hausdorff dist normconsts:\n  dist = %f\n", hausdorff_dist_normconsts);
    #endif

    UINT is_contspec_small = 1;                                                                     
    for (UINT i=0; i<params_i.M; i++){
        if (CABS(contspec_r[i]) > err_bnd_contspec){
            is_contspec_small = 0;
        }
    }

    UINT are_spurious_bound_states_small = 1;
    for (UINT i=0; i<(K_r-params_i.K); i++){
        if (CABS(bound_states_r[i]) > err_bnd_spurious_bound_states){
            are_spurious_bound_states_small = 0;
        }
    }

    if (is_bsr_in_tolerance &&
        is_ncr_in_tolerance &&
        is_contspec_small &&
        are_spurious_bound_states_small) 
    {
        ret_code = SUCCESS;
    } 
    else {
        ret_code = FNFT_EC_TEST_FAILED;
    }
    

leave_fun:
    free(q);
    free(contspec_r);
    free(bound_states_r);
    free(normconsts_r);

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}