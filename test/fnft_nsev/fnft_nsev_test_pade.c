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
 * Igor Chekhovskoy 2026.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"

#include <stdio.h>
#include <stdlib.h>

static INT test_zero_signal(const nse_discretization_t discretization,
        const UINT pade_degree, const INT kappa, const INT normalization_flag)
{
    const UINT D = 16, M = 5;
    const REAL T[2] = {-1.0, 1.0};
    const REAL XI[2] = {-2.0, 2.0};
    COMPLEX q[16] = {0.0};
    COMPLEX spectrum[15];
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    UINT i;
    INT ret_code;

    opts.discretization = discretization;
    opts.pade_degree = pade_degree;
    opts.contspec_type = nsev_cstype_BOTH;
    opts.normalization_flag = normalization_flag;
    ret_code = fnft_nsev(D, q, T, M, spectrum, XI, NULL, NULL, NULL,
            kappa, &opts);
    if (ret_code != SUCCESS)
        return ret_code;
    for (i = 0; i < M; i++) {
        if (FNFT_CABS(spectrum[i]) > 2e-11
                || FNFT_FABS(FNFT_CABS(spectrum[M + i]) - 1.0) > 2e-11
                || FNFT_CABS(spectrum[2*M + i]) > 2e-11) {
            fprintf(stderr, "Padé zero-signal failure: degree=%lu kappa=%d norm=%d i=%lu rho=%.3e a=(%.6e,%.6e) b=%.3e\n",
                    (unsigned long)pade_degree, (int)kappa,
                    (int)normalization_flag, (unsigned long)i,
                    FNFT_CABS(spectrum[i]), FNFT_CREAL(spectrum[M + i]),
                    FNFT_CIMAG(spectrum[M + i]),
                    FNFT_CABS(spectrum[2*M + i]));
            return E_TEST_FAILED;
        }
    }
    return SUCCESS;
}

static REAL continuous_spectrum_error(const nsev_testcases_t testcase,
        const UINT D, const nse_discretization_t discretization,
        const UINT pade_degree, INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL, *computed = NULL;
    REAL T[2], XI[2], error = 0.0, norm = 0.0;
    UINT M, K, i;
    INT kappa;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();

    *ret_code = nsev_testcases(testcase, D, &q, T, &M, &exact, &ab, XI,
            &K, &bound_states, &normconsts, &residues, &kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    computed = malloc(M*sizeof(COMPLEX));
    if (computed == NULL) {
        *ret_code = E_NOMEM;
        goto release_mem;
    }
    opts.discretization = discretization;
    opts.pade_degree = pade_degree;
    opts.contspec_type = nsev_cstype_REFLECTION_COEFFICIENT;
    *ret_code = fnft_nsev(D, q, T, M, computed, XI, NULL, NULL, NULL,
            kappa, &opts);
    if (*ret_code != SUCCESS)
        goto release_mem;
    for (i = 0; i < M; i++) {
        error += FNFT_CABS(computed[i] - exact[i]);
        norm += FNFT_CABS(exact[i]);
    }
    error /= 1.0 + norm;

release_mem:
    free(q);
    free(exact);
    free(ab);
    free(bound_states);
    free(normconsts);
    free(residues);
    free(computed);
    return error;
}

static INT test_convergence(const nsev_testcases_t testcase,
        const nse_discretization_t discretization, const UINT pade_degree,
        const UINT coarse_D, const REAL minimum_ratio)
{
    INT ret_code;
    const REAL coarse = continuous_spectrum_error(testcase, coarse_D,
            discretization, pade_degree, &ret_code);
    REAL fine;

    if (ret_code != SUCCESS)
        return ret_code;
    fine = continuous_spectrum_error(testcase, 2*coarse_D, discretization,
            pade_degree, &ret_code);
    if (ret_code != SUCCESS)
        return ret_code;
    if (!(fine > 0.0 && coarse/fine >= minimum_ratio)) {
        fprintf(stderr, "Padé convergence failure: degree=%lu coarse=%.3e fine=%.3e ratio=%.2f\n",
                (unsigned long)pade_degree, coarse, fine, coarse/fine);
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

static REAL bound_state_error(const UINT D, INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL, *computed = NULL;
    REAL T[2], XI[2], error = 0.0;
    UINT M, K, K_exact, i, j;
    INT kappa;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();

    *ret_code = nsev_testcases(nsev_testcases_SECH_FOCUSING, D, &q, T,
            &M, &exact, &ab, XI, &K_exact, &bound_states, &normconsts,
            &residues, &kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    opts.discretization = nse_discretization_FES4_PADE;
    opts.pade_degree = 2;
    opts.bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE;
    K = fnft_nsev_max_K(D, &opts);
    computed = malloc(K*sizeof(COMPLEX));
    if (computed == NULL) {
        *ret_code = E_NOMEM;
        goto release_mem;
    }
    *ret_code = fnft_nsev(D, q, T, 0, NULL, XI, &K, computed, NULL,
            kappa, &opts);
    if (*ret_code != SUCCESS)
        goto release_mem;
    if (K != K_exact) {
        fprintf(stderr, "Padé bound-state count failure: D=%lu computed=%lu exact=%lu\n",
                (unsigned long)D, (unsigned long)K, (unsigned long)K_exact);
        *ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    for (i = 0; i < K; i++) {
        REAL minimum = INFINITY;
        for (j = 0; j < K_exact; j++) {
            const REAL distance = FNFT_CABS(computed[i] - bound_states[j]);
            if (distance < minimum)
                minimum = distance;
        }
        if (minimum > error)
            error = minimum;
    }

release_mem:
    free(q);
    free(exact);
    free(ab);
    free(bound_states);
    free(normconsts);
    free(residues);
    free(computed);
    return error;
}

static INT test_bound_state_localization(void)
{
    INT ret_code;
    const REAL coarse = bound_state_error(256, &ret_code);
    REAL fine;

    if (ret_code != SUCCESS)
        return ret_code;
    fine = bound_state_error(512, &ret_code);
    if (ret_code != SUCCESS)
        return ret_code;
    if (!(fine > 0.0 && fine < 1e-3 && coarse/fine >= 3.0)) {
        fprintf(stderr, "Padé bound-state localization failure: coarse=%.3e fine=%.3e ratio=%.2f\n",
                coarse, fine, coarse/fine);
        return E_TEST_FAILED;
    }
    return SUCCESS;
}

INT main(void)
{
    INT ret_code, kappa, normalization_flag;

    for (normalization_flag = 0; normalization_flag <= 1;
            normalization_flag++) {
        for (kappa = -1; kappa <= 1; kappa += 2) {
            ret_code = test_zero_signal(nse_discretization_FES4_PADE, 2,
                    kappa, normalization_flag);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            ret_code = test_zero_signal(nse_discretization_FES6_PADE, 4,
                    kappa, normalization_flag);
            if (ret_code != SUCCESS)
                return EXIT_FAILURE;
            if (normalization_flag == 0 && kappa == 1) {
                ret_code = test_zero_signal(nse_discretization_FES6_PADE, 7,
                        kappa, normalization_flag);
                if (ret_code != SUCCESS)
                    return EXIT_FAILURE;
            }
        }
    }

    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,
            nse_discretization_FES4_PADE, 2, 512, 8.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_convergence(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            nse_discretization_FES6_PADE, 3, 128, 24.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,
            nse_discretization_FES6_PADE, 4, 512, 24.0);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_bound_state_localization();
    return ret_code == SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
