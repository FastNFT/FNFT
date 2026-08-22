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
#include "fnft__errwarn.h"

#include <stdio.h>
#include <stdlib.h>

static REAL continuous_spectrum_error(const nsev_testcases_t testcase,
        const UINT D, fnft_nsev_opts_t * const opts, INT * const ret_code)
{
    COMPLEX *q = NULL, *exact = NULL, *ab = NULL, *bound_states = NULL;
    COMPLEX *normconsts = NULL, *residues = NULL, *computed = NULL;
    REAL T[2], XI[2], error = 0.0, norm = 0.0;
    UINT M, K, i;
    INT kappa;

    *ret_code = nsev_testcases(testcase,D,&q,T,&M,&exact,&ab,XI,&K,
            &bound_states,&normconsts,&residues,&kappa);
    if (*ret_code != SUCCESS)
        goto release_mem;
    computed = malloc(M*sizeof(COMPLEX));
    if (computed == NULL) {
        *ret_code = E_NOMEM;
        goto release_mem;
    }
    opts->contspec_type = nsev_cstype_REFLECTION_COEFFICIENT;
    *ret_code = fnft_nsev(D,q,T,M,computed,XI,NULL,NULL,NULL,kappa,opts);
    if (*ret_code != SUCCESS)
        goto release_mem;
    for (i=0; i<M; i++) {
        error += CABS(computed[i]-exact[i]);
        norm += CABS(exact[i]);
    }
    error /= 1.0+norm;

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
        const UINT coarse_D, const UINT richardson_flag,
        const REAL minimum_ratio)
{
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    INT ret_code;
    REAL coarse, fine;

    opts.discretization = nse_discretization_ES6;
    opts.richardson_extrapolation_flag = richardson_flag;
    coarse = continuous_spectrum_error(testcase,coarse_D,&opts,&ret_code);
    CHECK_RETCODE(ret_code, leave_fun);
    fine = continuous_spectrum_error(testcase,2*coarse_D,&opts,&ret_code);
    CHECK_RETCODE(ret_code, leave_fun);
    if (!(fine > 0.0 && coarse/fine >= minimum_ratio)) {
        fprintf(stderr,"ES6 convergence failure: testcase=%d RE=%lu coarse=%.3e fine=%.3e ratio=%.2f\n",
                (int)testcase,(unsigned long)richardson_flag,coarse,fine,
                coarse/fine);
        return E_TEST_FAILED;
    }
leave_fun:
    return ret_code;
}

static INT test_full_focusing(void)
{
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    REAL error_bounds[6] = {
        2.0e-3, 2.0e-3, 2.0e-3, 2.0e-3, 2.0e-10, 2.0e-3
    };
    INT ret_code;

    opts.discretization = nse_discretization_ES6;
    opts.bound_state_localization = nsev_bsloc_NEWTON;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_FOCUSING2,
            1024,error_bounds,&opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = 0;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_FOCUSING2,
            1024,error_bounds,&opts);
leave_fun:
    return ret_code;
}

static INT test_pade_consistency(void)
{
    fnft_nsev_opts_t es6 = fnft_nsev_default_opts();
    fnft_nsev_opts_t pade3 = fnft_nsev_default_opts();
    fnft_nsev_opts_t pade4 = fnft_nsev_default_opts();
    INT ret_code;
    REAL es6_error, pade3_error, pade4_error;

    es6.discretization = nse_discretization_ES6;
    pade3.discretization = nse_discretization_FES6_PADE;
    pade3.pade_degree = 3;
    pade4.discretization = nse_discretization_FES6_PADE;
    pade4.pade_degree = 4;
    es6_error = continuous_spectrum_error(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            256,&es6,&ret_code);
    CHECK_RETCODE(ret_code, leave_fun);
    pade3_error = continuous_spectrum_error(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            256,&pade3,&ret_code);
    CHECK_RETCODE(ret_code, leave_fun);
    pade4_error = continuous_spectrum_error(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            256,&pade4,&ret_code);
    CHECK_RETCODE(ret_code, leave_fun);
    if (!(es6_error < 1e-4 && pade3_error < 1e-4 && pade4_error < 1e-4)) {
        fprintf(stderr,"ES6 consistency failure: exact=%.3e pade3=%.3e pade4=%.3e\n",
                es6_error,pade3_error,pade4_error);
        return E_TEST_FAILED;
    }
leave_fun:
    return ret_code;
}

INT main(void)
{
    INT ret_code;

    ret_code = test_convergence(nsev_testcases_SECH_FOCUSING_CONTSPEC,
            512,0,40.0);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,512,0,40.0);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_convergence(nsev_testcases_SECH_DEFOCUSING,256,1,128.0);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_full_focusing();
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_pade_consistency();
    CHECK_RETCODE(ret_code, failure);
    return EXIT_SUCCESS;
failure:
    return EXIT_FAILURE;
}
