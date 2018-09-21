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
* Sander Wahls (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__nsep_testcases.h"
#include "fnft__misc.h" // for misc_filter
#include "fnft__nse_discretization.h" // for nse_discretization_degree
#ifdef DEBUG
#include <stdio.h> // for printf
#endif

INT nsep_testcases(nsep_testcases_t tc, const UINT D,
    COMPLEX ** const q_ptr, REAL * const T,
    UINT * const K_ptr, COMPLEX ** const mainspec_ptr,
    UINT * const M_ptr, COMPLEX ** const auxspec_ptr,
    REAL ** const sheet_indices_ptr, INT * const kappa_ptr,
    REAL * const remove_box)
{
    UINT i, j;
    INT ret_code;
    REAL eps_t;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q_ptr == NULL)
        return E_INVALID_ARGUMENT(q_ptr);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);
    if (mainspec_ptr == NULL)
        return E_INVALID_ARGUMENT(mainspec_ptr);
    if (auxspec_ptr == NULL)
        return E_INVALID_ARGUMENT(auxspec_ptr);
    if (sheet_indices_ptr == NULL)
        return E_INVALID_ARGUMENT(sheet_indices_ptr);
    if (kappa_ptr == NULL)
        return E_INVALID_ARGUMENT(kappa_ptr);
    if (remove_box == NULL)
        return E_INVALID_ARGUMENT(remove_box);

    // Set the number of points in the continuous spectrum *M_ptr and the
    // number of bound states *K_ptr (needed for proper allocation)
    switch (tc) {

    case nsep_testcases_PLANE_WAVE_FOCUSING:
        *K_ptr = 100; // compute many -> unneeded ones get filtered later
        *M_ptr = *K_ptr - 2;
        break;

    case nsep_testcases_CONSTANT_DEFOCUSING:
        *K_ptr = 100; // compute many
        *M_ptr = *K_ptr;
        break;

    default:
        return E_INVALID_ARGUMENT(tc);
    }

    // Default remove_box: Nothing is removed.
    remove_box[0] = 0;
    remove_box[1] = 0;
    remove_box[2] = 0;
    remove_box[3] = 0;

    // Allocate memory for results
    *q_ptr = malloc(D * sizeof(COMPLEX));
    if (*q_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_1;
    }
    *mainspec_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
    if (*mainspec_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_2;
    }

    *auxspec_ptr = malloc((*M_ptr) * sizeof(COMPLEX));
    if (*auxspec_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_3;
    }
    *sheet_indices_ptr = malloc((*M_ptr) * sizeof(COMPLEX));
    if (*sheet_indices_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_4;
    }

    // generate test case
    switch (tc) {

    case nsep_testcases_PLANE_WAVE_FOCUSING:

        T[0] = 0.0; // location of first sample and begin of period
        T[1] = 2.0*PI; // end of period
        eps_t = 2.0*PI/D; // location of last sample is T[1]-eps_t

        for (i=0; i<D; i++)
            (*q_ptr)[i] = 2.0*CEXP( 3.0*I*(T[0] + i*eps_t) );

        // Make sure that *K_ptr and *M_ptr are multiples of 2
        if ((int)*K_ptr % 2 != 0) {
            ret_code = E_ASSERTION_FAILED;
            goto release_mem_4;
        }
        if ((int)*M_ptr % 2 != 0) {
            ret_code = E_ASSERTION_FAILED;
            goto release_mem_4;
        }

        // The main spectrum (first *K_ptr points only)
        for (i=0; i<*K_ptr; i+=2) {
            j = i/2;
            (*mainspec_ptr)[i] = -1.5 + I*CSQRT(4.0 - j*j/(REAL)4.0);
            (*mainspec_ptr)[i+1] = -1.5 - I*CSQRT(4.0 - j*j/(REAL)4.0);
        }

        // All main spectrum points are double expect for the ones with the
        // maximal imaginary part. Each of the double main spectrum points
        // traps a point from the auxiliary spectrum.
        for (i=0; i<*M_ptr; i+=2) {
            j = i/2 + 1; // "+1" => skip max imag parts
            (*auxspec_ptr)[i] = -1.5 + I*CSQRT(4.0 - j*j/(REAL)4.0);
            (*auxspec_ptr)[i+1] = -1.5 - I*CSQRT(4.0 - j*j/(REAL)4.0);
        }

        // The values at -1.5 converge very slowly => remove
        remove_box[0] = -1.5 - 0.1;
        remove_box[1] = -1.5 + 0.1,
        remove_box[2] = -0.1;
        remove_box[3] = 0.1,

        *kappa_ptr = +1;

        break;

   case nsep_testcases_CONSTANT_DEFOCUSING:

    /*  The following Matlab code was used to construct this test case. Note
        that Matlab misses two points in the main spectrum ([2] and [3] below).

        syms lam
        q = (1+2j)/sym(5);
        T = expm([-1j*lam q ; conj(q) 1j*lam]);

        disp('Main spectrum');
        del = trace(T);
        ret = solve(del == 2, lam, 'ReturnConditions', true);
        for i=1:length(ret.lam)
            [ret.lam(i) ret.conditions(i)]
        end
        ret = solve(del == -2, lam, 'ReturnConditions', true);
        for i=1:length(ret.lam)
            [ret.lam(i) ret.conditions(i)]
        end

        disp('Auxiliary spectrum');
        b = T(2,1);
        ret = solve(b == 0, lam, 'ReturnConditions', true);
        for i=1:length(ret.lam)
            [ret.lam(i) ret.conditions(i)]
        end
        */

        T[0] = 0; // location of 1st sample
        T[1] = 1.0; // end of period
        eps_t = (T[1] - T[0])/D; // location of last sample is T[1]-eps_t

        for (i=0; i<D; i++)
            (*q_ptr)[i] = (1.0 + 2.0*I)/5.0;

        // Main spectrum
        REAL pi2 = POW(PI, 2.0);
        (*mainspec_ptr)[0] = 1.0/SQRT(5.0);
        (*mainspec_ptr)[1] = -(*mainspec_ptr)[0];
        (*mainspec_ptr)[2] = SQRT(5.0*pi2 + 1.0)/SQRT(5.0);
        (*mainspec_ptr)[3] = -(*mainspec_ptr)[2];
        for (j=1;; j++) {
            i = 3+4*j;
            if (i >= *K_ptr)
                break;
            (*mainspec_ptr)[i-3] = SQRT(20.0*pi2*j*j + 1.0)/SQRT(5.0);
            (*mainspec_ptr)[i-2] = -(*mainspec_ptr)[i-3];
            (*mainspec_ptr)[i-1] = SQRT(20.0*pi2*j*j + 20.0*pi2*j + 5.0*pi2 + 1.0)/SQRT(5.0);
            (*mainspec_ptr)[i] = -(*mainspec_ptr)[i-1];
        }
        *K_ptr = i-4;

        // Auxiliary spectrum
        (*auxspec_ptr)[0] = SQRT(5.0*pi2 + 1.0)/SQRT(5.0);
        (*auxspec_ptr)[1] = -(*auxspec_ptr)[0];
        for (j=1;; j++) {
            i = 1 + 4*j;
            if (i >= *K_ptr)
                break;
            (*auxspec_ptr)[i-3] = SQRT(20.0*pi2*j*j + 1.0)/SQRT(5.0);
            (*auxspec_ptr)[i-2] = -(*auxspec_ptr)[i-3];
            (*auxspec_ptr)[i-1] = SQRT(20.0*pi2*j*j + 20.0*pi2*j + 5.0*pi2 + 1.0)/SQRT(5.0);
            (*auxspec_ptr)[i] = -(*auxspec_ptr)[i-1];
        }
        *M_ptr = i-4;

        *kappa_ptr = -1;

        break;

    default: // unknown test case

        ret_code = E_INVALID_ARGUMENT(tc);
        goto release_mem_4;
    }

    return SUCCESS;

    // the code below is only executed if an error occurs

release_mem_4:
    free(*sheet_indices_ptr);
release_mem_3:
    free(*auxspec_ptr);
release_mem_2:
    free(*mainspec_ptr);
release_mem_1:
    free(*q_ptr);

    return ret_code;
}

// Compares computed with true nonlinear Fourier spectrum.
static INT nsep_compare_nfs(const UINT K1, const UINT K2,
    COMPLEX const * const mainspec_1,
    COMPLEX const * const mainspec_2,
    const UINT M1, const UINT M2,
    COMPLEX const * const auxspec_1,
    COMPLEX const * const auxspec_2,
    REAL const * const sheet_indices_1,
    REAL const * const sheet_indices_2,
    REAL *dists)
{
    // Check last argument
    if (dists == NULL)
        return E_INVALID_ARGUMENT(dists);

    if (sheet_indices_1 != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices_1, "Pass NULL");
    if (sheet_indices_2 != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices_2, "Pass NULL");

    // Compare main spectra
    if (K1 == 0 && K2 == 0) {
        dists[0] = 0.0;
    } else if (mainspec_1 == NULL || mainspec_2 == NULL || K1 == 0
    || K2 == 0) {
        dists[0] = NAN;
    } else {
        dists[0] = misc_hausdorff_dist(K1, mainspec_1, K2, mainspec_2);
    }

    // Compare auxiliary spectra
    if (M1 == 0 && M2 == 0) {
        dists[1] = 0.0;
    } else if (auxspec_1 == NULL || auxspec_2 == NULL || M1 == 0
    || M2 == 0) {
        dists[1] = NAN;
    } else {
        dists[1] = misc_hausdorff_dist(M1, auxspec_1, M2, auxspec_2);
    }

    // Not yet implemented
    dists[2] = 0;

    return SUCCESS;
}

INT nsep_testcases_test_fnft(nsep_testcases_t tc, UINT D, REAL error_bounds[3],
fnft_nsep_opts_t * opts_ptr) {
    COMPLEX * q = NULL;
    COMPLEX * mainspec = NULL;
    COMPLEX * auxspec = NULL;
    REAL * sheet_indices = NULL;
    REAL T[2];
    COMPLEX * mainspec_exact = NULL;
    COMPLEX * auxspec_exact = NULL;
    REAL * sheet_indices_exact = NULL;
    REAL remove_box[4];
    UINT K, K_exact, M, M_exact;
    INT kappa = 0;
    fnft_nsep_opts_t default_opts;
    REAL errs[3];
    INT ret_code;

    // Load test case
    ret_code = nsep_testcases(tc, D, &q, T, &K_exact, &mainspec_exact,
        &M_exact, &auxspec_exact, &sheet_indices_exact, &kappa, remove_box);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);

    // Allocate memory
    if (opts_ptr == NULL) {
        default_opts = fnft_nsep_default_opts();
        opts_ptr = &default_opts;
    }
    K = 2*nse_discretization_degree(opts_ptr->discretization)*D + 1;
    M = K;
    mainspec = malloc(K * sizeof(COMPLEX));
    auxspec = malloc(M * sizeof(COMPLEX));
    if ( q == NULL || mainspec == NULL || auxspec == NULL ) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Compute the NFT (sheet_indices stay NULL because they are not yet
    // implemented)
    ret_code = fnft_nsep(D, q, T, &K, mainspec, &M, auxspec, sheet_indices,
        kappa, opts_ptr);
    CHECK_RETCODE(ret_code, release_mem);

    // Filter exact spectra as well
    ret_code = misc_filter(&K_exact, mainspec_exact, NULL,
        opts_ptr->bounding_box);
    CHECK_RETCODE(ret_code, release_mem);
    ret_code = misc_filter(&M_exact, auxspec_exact, NULL,
        opts_ptr->bounding_box);
    CHECK_RETCODE(ret_code, release_mem);

    // Apply the remove box to all spectra
    ret_code = misc_filter_inv(&K_exact, mainspec_exact, NULL, remove_box);
    CHECK_RETCODE(ret_code, release_mem);
    ret_code = misc_filter_inv(&K, mainspec, NULL, remove_box);
    CHECK_RETCODE(ret_code, release_mem);
    ret_code = misc_filter_inv(&M_exact, auxspec_exact, NULL, remove_box);
    CHECK_RETCODE(ret_code, release_mem);
    ret_code = misc_filter_inv(&M, auxspec, NULL, remove_box);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the errors
    ret_code = nsep_compare_nfs(K, K_exact, mainspec, mainspec_exact,
        M, M_exact, auxspec, auxspec_exact, NULL /* not yet implemented */,
        NULL /* not yet implemented */,  errs);
    CHECK_RETCODE(ret_code, release_mem);

#ifdef DEBUG
    printf(
        "nsep_testcases_test_fnft: %2.1e<=%2.1e %2.1e<=%2.1e %2.1e<=%2.1e\n",
        errs[0], error_bounds[0],
        errs[1], error_bounds[1],
        errs[2], error_bounds[2]
    );
    //print_buf2(K_exact, mainspec_exact, "mainspec_exact");
    //print_buf2(K, mainspec, "mainspec");
    //print_buf2(M_exact, auxspec_exact, "auxspec_exact");
    //print_buf2(M, auxspec, "auxspec");
#endif

    // Check if the errors are below the specified bounds. Organized such that
    // the line number tells us which error was too high. The conditions are
    // written in this way to ensure that they fail if an error is NAN.
    if (!(errs[0] <= error_bounds[0])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[1] <= error_bounds[1])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[2] <= error_bounds[2])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    ///// Clean up /////

release_mem:
    free(q);
    free(mainspec);
    free(auxspec);
    free(mainspec_exact);
    free(auxspec_exact);
    free(sheet_indices_exact);

    return ret_code;
}
