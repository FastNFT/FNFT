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
* Fabian Fischer (Hiwi KIT) 2026.
*/

#define FNFT_ENABLE_SHORT_NAMES
#define DEBUG

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

INT inverse_kdvv_testcases(inverse_kdvv_testcases_t tc, 
                    const UINT D,
                    REAL * const T,
                    UINT * const M_ptr, 
                    REAL * const XI, 
                    UINT * const K_ptr,
                    COMPLEX ** const bound_states_ptr,
                    COMPLEX ** const normconsts_ptr)
{
    INT ret_code = SUCCESS;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (bound_states_ptr == NULL)
        return E_INVALID_ARGUMENT(bound_state_ptr);
    if (normconsts_ptr == NULL)
        return E_INVALID_ARGUMENT(normconst_ptr);

    // Set the number of points in the continuous spectrum *M_ptr and the
    // number of bound states *K_ptr (needed for proper allocation)
    switch (tc) {

        case inverse_kdvv_testcases_5_bound_states:
            *M_ptr = 16;
            *K_ptr = 5;
            break;

        case inverse_kdvv_testcases_19_bound_states:
            *M_ptr = 10;
            *K_ptr = 19;
            break;

        case inverse_kdvv_testcases_8_bound_states_asym:
            *M_ptr = 10;
            *K_ptr = 8;
            break;

        default:
            return E_INVALID_ARGUMENT(tc);
    }

    // Allocate memory for results
    if (*K_ptr > 0) {
        *bound_states_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
        CHECK_NOMEM(*bound_states_ptr, ret_code, release_mem_1);
        *normconsts_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
        CHECK_NOMEM(*normconsts_ptr, ret_code, release_mem_2);
    } else {
        *bound_states_ptr = NULL;
        *normconsts_ptr = NULL;
    }

    // generate test case
    switch (tc) {

        case inverse_kdvv_testcases_5_bound_states:

            (*bound_states_ptr)[0] = I*SQRT(5.0/2.0);
            (*bound_states_ptr)[1] = I*SQRT(4.0/2.0);
            (*bound_states_ptr)[2] = I*SQRT(3.0/2.0);
            (*bound_states_ptr)[3] = I*SQRT(2.0/2.0);
            (*bound_states_ptr)[4] = I*SQRT(1.0/2.0);

            (*normconsts_ptr)[0] = 1*1e7;
            (*normconsts_ptr)[1] = -1*1e-5;
            (*normconsts_ptr)[2] = 1*1;
            (*normconsts_ptr)[3] = -1*0.1;
            (*normconsts_ptr)[4] = 1*10;

            XI[0] = 0.5;
            XI[1] = 23.0;

            T[0] = -20.0;
            T[1] = 20.0;

            break;

        case inverse_kdvv_testcases_19_bound_states:

            (*bound_states_ptr)[0] = I*SQRT(40.0/2.0);
            (*bound_states_ptr)[1] = I*SQRT(39.0/2.0);
            (*bound_states_ptr)[2] = I*SQRT(38.0/2.0);
            (*bound_states_ptr)[3] = I*SQRT(37.0/2.0);
            (*bound_states_ptr)[4] = I*SQRT(36.0/2.0);
            (*bound_states_ptr)[5] = I*SQRT(30.0/2.0);
            (*bound_states_ptr)[6] = I*SQRT(29.0/2.0);
            (*bound_states_ptr)[7] = I*SQRT(28.0/2.0);
            (*bound_states_ptr)[8] = I*SQRT(27.0/2.0);
            (*bound_states_ptr)[9] = I*SQRT(26.0/2.0);
            (*bound_states_ptr)[10] = I*SQRT(20.0/2.0);
            (*bound_states_ptr)[11] = I*SQRT(19.0/2.0);
            (*bound_states_ptr)[12] = I*SQRT(18.0/2.0);
            (*bound_states_ptr)[13] = I*SQRT(17.0/2.0);
            (*bound_states_ptr)[14] = I*SQRT(16.0/2.0);
            (*bound_states_ptr)[15] = I*SQRT(10.0/2.0);
            (*bound_states_ptr)[16] = I*SQRT(9.0/2.0);
            (*bound_states_ptr)[17] = I*SQRT(8.0/2.0);
            (*bound_states_ptr)[18] = I*SQRT(7.0/2.0);

            (*normconsts_ptr)[0] = 1e20;
            (*normconsts_ptr)[1] = -1e-7;
            (*normconsts_ptr)[2] = 1e5;
            (*normconsts_ptr)[3] = -1e3;
            (*normconsts_ptr)[4] = 1e1;
            (*normconsts_ptr)[5] = -1e0;
            (*normconsts_ptr)[6] = 1e2;
            (*normconsts_ptr)[7] = -1e4;
            (*normconsts_ptr)[8] = 1e-6;
            (*normconsts_ptr)[9] = -1e8;
            (*normconsts_ptr)[10] = 1e2;
            (*normconsts_ptr)[11] = -1e4;
            (*normconsts_ptr)[12] = 1e6;
            (*normconsts_ptr)[13] = -1e8;
            (*normconsts_ptr)[14] = 1e-10;
            (*normconsts_ptr)[15] = -1e7;
            (*normconsts_ptr)[16] = 1e-6;
            (*normconsts_ptr)[17] = -1e5;
            (*normconsts_ptr)[18] = 1e-9;

            XI[0] = 0.5;
            XI[1] = 23.0;

            T[0] = -15.0;
            T[1] = 15.0;

            break;

        case inverse_kdvv_testcases_8_bound_states_asym:

            (*bound_states_ptr)[0] = I*SQRT(40.0/2.0);
            (*bound_states_ptr)[1] = I*SQRT(35.0/2.0);
            (*bound_states_ptr)[2] = I*SQRT(28.0/2.0);
            (*bound_states_ptr)[3] = I*SQRT(23.0/2.0);
            (*bound_states_ptr)[4] = I*SQRT(19.0/2.0);
            (*bound_states_ptr)[5] = I*SQRT(16.0/2.0);
            (*bound_states_ptr)[6] = I*SQRT(12.0/2.0);
            (*bound_states_ptr)[7] = I*SQRT(4.0/2.0);

            (*normconsts_ptr)[0] = 1e-3;
            (*normconsts_ptr)[1] = -1e-17;
            (*normconsts_ptr)[2] = 1e-15;
            (*normconsts_ptr)[3] = -1e-20;
            (*normconsts_ptr)[4] = 1e-1;
            (*normconsts_ptr)[5] = -1e-5;
            (*normconsts_ptr)[6] = 1e-12;
            (*normconsts_ptr)[7] = -1e-14;

            XI[0] = 0.5;
            XI[1] = 23.0;

            T[0] = -18.0;
            T[1] = 5.0;

            break;

        default: // unknown test case

            ret_code = E_INVALID_ARGUMENT(tc);
            goto release_mem_2;
    }

    return SUCCESS;

    // the code below is only executed if an error occurs

release_mem_2:
    free(*normconsts_ptr);
release_mem_1:
    free(*bound_states_ptr);

    return ret_code;
}

// Compares computed with exact nonlinear Fourier spectrum.
static INT inverse_kdvv_check_nfs(const UINT M, const UINT K_computed, const UINT K_exact,
                            COMPLEX const * const contspec_computed,
                            COMPLEX const * const bound_states_computed,
                            COMPLEX const * const bound_states_exact,
                            COMPLEX const * const norming_constants_computed,
                            COMPLEX const * const norming_constants_exact,
                            REAL dists[4])
{
    // Check last argument
    if (dists == NULL)
        return E_INVALID_ARGUMENT(dists);

    
    // Check matching of the computed bound states
    // bound_states_computed is sorted in ascending order, K_computed is the number of found bound states by fnft_kdvv
    // take only the last K_exact bound states, because they are the bigger ones and more likely the ones
    // which corresponds to the initial given bound states
    COMPLEX const * const candidate_bound_states_r_ptr = &bound_states_computed[K_computed-K_exact];
    dists[0] = misc_hausdorff_dist_normed(K_exact, bound_states_exact, K_exact, 
                                                            candidate_bound_states_r_ptr);

    COMPLEX const * const candidate_normconsts_r_ptr = &norming_constants_computed[K_computed-K_exact];
    dists[1] = misc_hausdorff_dist_normed(K_exact, norming_constants_exact, K_exact, 
                                                            candidate_normconsts_r_ptr);                                                          

    dists[2] = 0;                                                                     
    for (UINT i=0; i<M; i++){
        if (CABS(contspec_computed[i]) > dists[2]){
            dists[2] = CABS(contspec_computed[i]);
        }
    }

    dists[3] = 0;
    for (UINT i=0; i<(K_computed-K_exact); i++){
        if (CABS(bound_states_computed[i]) > dists[3]){
            dists[3] = CABS(bound_states_computed[i]);
        }
    }                       

    return SUCCESS;
}



INT inverse_kdvv_testcases_test_fnft(inverse_kdvv_testcases_t tc, UINT D,
                             const REAL error_bounds[6], void * const opts)
{
    COMPLEX * q = NULL;
    COMPLEX * contspec_computed = NULL;
    COMPLEX * bound_states_computed = NULL;
    COMPLEX * norming_constants_computed = NULL;
    REAL T[2];
    REAL XI[2];
    COMPLEX * bound_states_exact = NULL;
    COMPLEX * norming_constants_exact = NULL;
    UINT K;
    UINT K_exact = 0;
    UINT M;
    REAL errors[4] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN};
    INT ret_code;

    // Check inputs: opts has not yet been implemented (state 08/2026)!
    if (!(opts == NULL))
        return E_INVALID_ARGUMENT(opts);

    // Load test case
    ret_code = inverse_kdvv_testcases(tc, D, T, &M, XI, &K_exact, &bound_states_exact, 
                                        &norming_constants_exact);
    CHECK_RETCODE(ret_code, release_mem);

#ifdef DEBUG
    for (UINT i=0; i<K_exact; i++) {
        printf("bound_state: %12.1e\n", CREAL(bound_states_exact[i])*1e8);
        printf("norming_constant: %12.1e\n", CABS(norming_constants_exact[i]));
    }
#endif

    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, release_mem);

    // Compute the inverse kdvv NFT
    ret_code = fnft_kdvv_inverse(   0, NULL, NULL, K_exact, bound_states_exact, 
                                    norming_constants_exact, D, q, T, NULL);
    CHECK_RETCODE(ret_code, release_mem);

    // Allocate memory for forward NFT
    contspec_computed = malloc(M * sizeof(COMPLEX));
    CHECK_NOMEM(contspec_computed, ret_code, release_mem);
    K = D;
    bound_states_computed = malloc(K * sizeof(COMPLEX));
    CHECK_NOMEM(bound_states_computed, ret_code, release_mem);
    norming_constants_computed = malloc(K * sizeof(COMPLEX));
    CHECK_NOMEM(norming_constants_computed, ret_code, release_mem);

    // Compute NFT out of results of the inverse
    fnft_kdvv_opts_t opts_fnft = fnft_kdvv_default_opts();

    ret_code = fnft_kdvv(D, q, T, M, contspec_computed, XI, &K, bound_states_computed, 
                                    norming_constants_computed, &opts_fnft);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the errors
    ret_code = inverse_kdvv_check_nfs(M, K, K_exact, contspec_computed,
                                bound_states_computed, bound_states_exact, norming_constants_computed, 
                                norming_constants_exact, errors);
    CHECK_RETCODE(ret_code, release_mem);

    // -- Check results --

    // Simple general tests
    for (UINT i=0; i<K_exact; i++){
        COMPLEX bsc = bound_states_computed[i];
        COMPLEX ncc = norming_constants_computed[i];

        UINT is_bsc_pure_imaginary = FABS(CREAL(bsc)) < 1e-8;
        UINT is_bsc_positive_imaginary = (CIMAG(bsc) > 0) && is_bsc_pure_imaginary;

        UINT is_ncc_real = FABS(CIMAG(ncc)) < 1e-8;

        if (is_bsc_pure_imaginary &&
            is_bsc_positive_imaginary &&
            is_ncc_real) 
        {
            ret_code = SUCCESS;
        } 
        else {
            ret_code = FNFT_EC_TEST_FAILED;
            break;
        }
    }
    CHECK_RETCODE(ret_code, release_mem);

    // Check if the errors are below the specified bounds. Organized such that
    // the line number tells us which error was too high. The conditions are
    // written in this way to ensure that they fail if an error is NAN.

    // deviation in bound states
    if (!(errors[0] <= error_bounds[0])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    // deviation in norming constants
    if (!(errors[1] <= error_bounds[1])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    // deviation in contspec: computed contspec should be zero everywhere
    if (!(errors[2] <= error_bounds[2])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    // spurious bound states:
    // bound states which arises by computation and which are not included in the exact bound states
    if (!(errors[3] <= error_bounds[3])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

#ifdef DEBUG
    for (UINT i=0; i<4; i++)
        printf("kdvv_testcases_test_fnft: error_bounds[%i] = %2.1e <= %2.1e\n",
               (int)i, errors[i], error_bounds[i]);
#endif

#ifdef DEBUG
        misc_print_buf(K_exact, norming_constants_exact, "norming_constants_exact");
        misc_print_buf(K_exact, norming_constants_computed, "norming_constants_computed");
#endif
    
    ///// Clean up /////

release_mem:
    free(q);
    free(contspec_computed);
    free(bound_states_computed);
    free(bound_states_exact);
    free(norming_constants_computed);
    free(norming_constants_exact);

    return ret_code;
}   