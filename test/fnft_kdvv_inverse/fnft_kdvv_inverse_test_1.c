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
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft__kdvv_inverse_testcases.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

static INT print_test_results(  COMPLEX * bound_states_r,
                                COMPLEX * normconsts_r,
                                COMPLEX * contspec_r,
                                REAL * XI,
                                const UINT M,
                                const UINT D,
                                const UINT K_r)
{
    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    FNFT_REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    printf("Continuous spectrum:\n");
    for (FNFT_UINT i=0; i<M; i++) {
        FNFT_REAL xi = XI[0] + i*eps_xi;
        printf("  continuous_spectrum(xi=%f) \t= %g + %gI\n",
            (double)xi,
            (double)FNFT_CREAL(contspec_r[i]),
            (double)FNFT_CIMAG(contspec_r[i])
        );
    }

    printf("Discrete spectrum:\n");
    for (FNFT_UINT i=0; i<K_r; i++) {
        printf("  bound state at %g + %gI with norming constant %g + %gI\n",
            (double)FNFT_CREAL(bound_states_r[i]),
            (double)FNFT_CIMAG(bound_states_r[i]),
            (double)FNFT_CREAL(normconsts_r[i]),
            (double)FNFT_CIMAG(normconsts_r[i])
        );
    }
}



static INT run_test(const UINT D,
                    REAL err_bnd_bound_states,
                    REAL err_bnd_spurious_bound_states,
                    REAL err_bnd_normconst,
                    REAL err_bnd_contspec)

{
INT ret_code = SUCCESS;
    
    // General parameters
    UINT M = 10;
    REAL XI[2] = {-2.0, 2.0};
    REAL T[2] = {-20.0, 20.0};

    COMPLEX * q = NULL;
    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);

    // Initialize variables for a specific test case without continuous spectrum
    COMPLEX * contspec_i = NULL;
    contspec_i = malloc(M * sizeof(COMPLEX));
    const UINT K_i = 5;
    COMPLEX bound_states_i[5] = {I*SQRT(1.0/2.0), I*SQRT(2.0/2.0), I*SQRT(3.0/2.0), I*SQRT(4.0/2.0), I*SQRT(5.0/2.0)};
    COMPLEX normconsts_i[5] = {1*10, -1*0.1, 1*1, -1*1e-5, 1*1e7};

    // const UINT K_i = 1;
    // COMPLEX bound_states_i[1] = {SQRT(1.0/2)*I};
    // COMPLEX normconsts_i[1] = {1*10};

    ret_code = fnft_kdvv_inverse(M, contspec_i, XI, K_i, bound_states_i, normconsts_i, D, q, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Prepare forward fnft_kdvv
    COMPLEX * contspec_r = NULL;
    contspec_r = malloc(M * sizeof(COMPLEX));
    CHECK_NOMEM(contspec_r, ret_code, leave_fun);

    UINT K_r = D;

    COMPLEX * bound_states_r = NULL;
    bound_states_r = malloc(K_r * sizeof(COMPLEX));
    CHECK_NOMEM(bound_states_r, ret_code, leave_fun);

    COMPLEX * normconsts_r = NULL;
    normconsts_r = malloc(K_r * sizeof(COMPLEX));
    CHECK_NOMEM(normconsts_r, ret_code, leave_fun);

    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();

    ret_code = fnft_kdvv(D, q, T, M, contspec_r, XI, &K_r, bound_states_r, normconsts_r, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // -- Check results --

    // Simple general tests
    for (UINT i=0; i<K_i; i++){
        COMPLEX bsi = bound_states_i[i];   
        COMPLEX nci = normconsts_i[i];
        COMPLEX bsr = bound_states_r[i];
        COMPLEX ncr = normconsts_r[i];

        UINT is_bsr_pure_imaginary = CABS(CREAL(bsr)) < 1e-8;
        UINT is_bsr_positive_imaginary = (CIMAG(bsr) > 0) && is_bsr_pure_imaginary;

        UINT is_ncr_real = CABS(CIMAG(ncr)) < 1e-8;

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
    // take only the last K_i bound states, because they are the bigger ones and more likely the ones
    // which corresponds to the initial given bound states
    COMPLEX * const candidate_bound_states_r_ptr = &bound_states_r[K_r-K_i];
    REAL hausdorff_dist_bound_states = misc_hausdorff_dist_normed(  K_i, bound_states_i, K_i, 
                                                                    candidate_bound_states_r_ptr);
    UINT is_bsr_in_tolerance = hausdorff_dist_bound_states < err_bnd_bound_states;                                                                

    COMPLEX * const candidate_normconsts_r_ptr = &normconsts_r[K_r-K_i];
    REAL hausdorff_dist_normconsts = misc_hausdorff_dist_normed(K_i, normconsts_i, K_i, 
                                                                candidate_normconsts_r_ptr);
    UINT is_ncr_in_tolerance = hausdorff_dist_normconsts < err_bnd_normconst;                                                         

    UINT is_contspec_small = 1;                                                                     
    for (UINT i=0; i<M; i++){
        if (CABS(contspec_r[i]) > err_bnd_contspec){
            is_contspec_small = 0;
        }
    }

    UINT are_spurious_bound_states_small = 1;
    for (UINT i=0; i<(K_r-K_i); i++){
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

#ifdef DEBUG
    printf("Number of bound_states:\n  K_r = %u\n", (unsigned int)K_r);
    printf("Hausdorff dist bound states:\n  dist = %f\n", hausdorff_dist_bound_states);
    printf("Hausdorff dist normconsts:\n  dist = %f\n", hausdorff_dist_normconsts);
    print_test_results(bound_states_r, normconsts_r, contspec_r, XI, M, D, K_r);
#endif    
    

leave_fun:
    free(q);
    free(contspec_i);
    free(contspec_r);
    free(bound_states_r);
    free(normconsts_r);

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}


INT main()
{
    INT ret_code = SUCCESS;
    
    UINT D = 256;
    REAL err_bnd_bound_states = 1.4e-3;
    REAL err_bnd_spurious_bound_states = 1e-2;
    REAL err_bnd_normconst = 4.5e-2;
    REAL err_bnd_contspec = 1e-1;

    
    ret_code = run_test(D, err_bnd_bound_states, err_bnd_spurious_bound_states, 
                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);

    
    // Check quadratic convergence
    D *= 2;
    err_bnd_bound_states /= 4;
    err_bnd_spurious_bound_states /= 4;
    err_bnd_normconst = 1.1e-2;
    err_bnd_contspec = 5e-3;

    ret_code = run_test(D, err_bnd_bound_states, err_bnd_spurious_bound_states, 
                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);

    D *= 2;
    err_bnd_bound_states /= 4;
    err_bnd_spurious_bound_states /= 4;
    err_bnd_normconst /= 4;
    err_bnd_contspec /= 4;

    ret_code = run_test(D, err_bnd_bound_states, err_bnd_spurious_bound_states, 
                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);
    

leave_fun:
    

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}




