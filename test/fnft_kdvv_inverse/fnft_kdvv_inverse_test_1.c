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


static INT run_test(const UINT D,
                    REAL err_bnd_bound_states,
                    REAL err_bnd_spurious_bound_states,
                    REAL err_bnd_normconst,
                    REAL err_bnd_contspec)

{
INT ret_code = SUCCESS;
    
    // General parameters
    UINT M = 0;
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
    COMPLEX normconsts_i[5] = {1*10, -1*0.1, 1*1, -1*0.00001, 1*10000000};

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

    // TODO: sorting Eigenvalues from bound_states_r?

    // Check result
    // bound_states_r is in ascending order -> start from the end
//     for (UINT i=0; i<K_i; i++){
//         COMPLEX bsi = bound_states_i[i];   
//         COMPLEX nci = normconsts_i[i];
//         COMPLEX bsr = bound_states_r[i];
//         COMPLEX ncr = normconsts_r[i];

//         UINT is_bsr_pure_imaginary = CABS(CREAL(bsr)) < 1e-9;
//         UINT is_bsr_positive_imaginary = (CIMAG(bsr) > 0) && is_bsr_pure_imaginary;
//         UINT is_bsr_in_tolerance = CABS(bsr - bsi)/CABS(bsr) < err_bnd_bound_states;

//         UINT are_spurious_bound_states_small = err_bnd_spurious_bound_states;

//         UINT is_ncr_real = CABS(CIMAG(ncr)) < 1e-9;
//         UINT is_ncr_in_tolerance = CABS(ncr - nci)/CABS(ncr) < err_bnd_normconst;

// // misc_hausdorff_dist

//         UINT is_contspec_small = err_bnd_contspec;


//         if (is_bsr_pure_imaginary &&
//             is_bsr_positive_imaginary &&
//             is_bsr_in_tolerance &&
//             is_ncr_real &&
//             is_ncr_in_tolerance) 
//         {
//             ret_code = SUCCESS;
//         } 
//         else {
//             ret_code = FNFT_EC_TEST_FAILED;
//             break;
//         }
//     }

    for (UINT i=0; i<K_i; i++){
        COMPLEX bsi = bound_states_i[i];   
        COMPLEX nci = normconsts_i[i];
        COMPLEX bsr = bound_states_r[i];
        COMPLEX ncr = normconsts_r[i];

        UINT is_bsr_pure_imaginary = CABS(CREAL(bsr)) < 1e-8;
        UINT is_bsr_positive_imaginary = (CIMAG(bsr) > 0) && is_bsr_pure_imaginary;
        UINT is_bsr_in_tolerance = CABS(bsr - bsi)/CABS(bsr) < 1e-3;

        UINT is_ncr_real = CABS(CIMAG(ncr)) < 1e-8;
        UINT is_ncr_in_tolerance = CABS(ncr - nci)/CABS(ncr) < 1e-2;

        if (is_bsr_pure_imaginary &&
            is_bsr_positive_imaginary &&
            is_bsr_in_tolerance &&
            is_ncr_real &&
            is_ncr_in_tolerance) 
        {
            ret_code = SUCCESS;
        } 
        else {
            ret_code = FNFT_EC_TEST_FAILED;
            break;
        }
    }
    

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


INT main()
{
    INT ret_code = SUCCESS;

    const UINT D = 1001;
    
    run_test(D, 1e-3, 1e-3, 1e-3, 1e-3);
    CHECK_RETCODE(ret_code, leave_fun);
    

leave_fun:
    

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}




