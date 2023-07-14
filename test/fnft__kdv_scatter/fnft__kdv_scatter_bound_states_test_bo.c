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
* Sander Wahls (TU Delft) 2023.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <assert.h>
#include "fnft__kdv_scatter.h"

static INT test_aux(const REAL ebounds[3], kdv_discretization_t discr,
                    const UINT normalization_flag)
{
    const REAL T[2] = {-6, 6};
    COMPLEX q[1024];
    COMPLEX * q_preprocessed = NULL;
    COMPLEX * r_preprocessed = NULL;
    COMPLEX lambda[1] = {I};
    COMPLEX a_aprime_b[3] = {0};
    UINT i;
    UINT vanilla_flag = 0;
    INT W = 0;
    INT * const W_ptr = (normalization_flag) ? &W : NULL;
    INT ret_code = SUCCESS;

    ret_code = kdv_discretization_vanilla_flag(&vanilla_flag, discr);
    CHECK_RETCODE(ret_code, leave_fun);

    const UINT D = sizeof(q)/sizeof(q[0]);
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    for (i=0; i<D; i++) {
        const REAL t = T[0] + i*eps_t;
        const COMPLEX s = misc_sech(t);
        const COMPLEX val = 2*s*s;
        q[i] = val;
    }

    const UINT upsampling_factor = kdv_discretization_upsampling_factor(discr);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(discr);
        goto leave_fun;
    }

    UINT Dsub = D;
    UINT first_last_index[2] = {0};
    ret_code = kdv_discretization_preprocess_signal(D, q, eps_t, 0, &Dsub, &q_preprocessed, &r_preprocessed, first_last_index, discr);
    CHECK_RETCODE(ret_code, leave_fun);

    const UINT D_effective = D * upsampling_factor;
    ret_code = fnft__kdv_scatter_bound_states(D_effective, q_preprocessed, r_preprocessed, T, 1, lambda, a_aprime_b, a_aprime_b+1, a_aprime_b+2, W_ptr, discr, 0);
    CHECK_RETCODE(ret_code, leave_fun);

    if (W_ptr != NULL) {
        const REAL scl = POW(2, W);
        a_aprime_b[0] *= scl;
        a_aprime_b[1] *= scl;
        // b is not scaled
    }

    const COMPLEX a_exact = 0;
    const COMPLEX aprime_exact = -I/2;
    const COMPLEX b_exact = 1;
    REAL errs[3] = { CABS(a_exact - a_aprime_b[0]),
                     CABS(aprime_exact - a_aprime_b[1]),
                     CABS(b_exact - a_aprime_b[2]) };

#ifdef DEBUG
    misc_print_buf(3, a_aprime_b , "a_aprime_b");
    printf("errs = [%g, %g, %g];\n", errs[0], errs[1], errs[2]);
#endif

    if (!(errs[0] <= ebounds[0])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (!(errs[1] <= ebounds[1])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (!(errs[2] <= ebounds[2])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    // repeat with skip_b_flag turned on
    a_aprime_b[0] = 0;
    a_aprime_b[1] = 0;
    a_aprime_b[2] = 0;
    ret_code = fnft__kdv_scatter_bound_states(D_effective, q_preprocessed, r_preprocessed, T, 1, lambda, a_aprime_b, a_aprime_b+1, a_aprime_b+2, W_ptr, discr, 1);
    CHECK_RETCODE(ret_code, leave_fun);

    if (W_ptr != NULL) {
        const REAL scl = POW(2, W);
        a_aprime_b[0] *= scl;
        a_aprime_b[1] *= scl;
    }

    errs[0] = CABS(a_exact - a_aprime_b[0]);
    errs[1] = CABS(aprime_exact - a_aprime_b[1]);

#ifdef DEBUG
    misc_print_buf(2, a_aprime_b , "a_aprime");
    printf("errs = [%g, %g];\n", errs[0], errs[1]);
#endif

    if (!(errs[0] <= ebounds[0])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (!(errs[1] <= ebounds[1])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(q_preprocessed);
    free(r_preprocessed);
    return ret_code;
}

static INT test(const REAL ebounds[3], kdv_discretization_t discr)
{
#ifdef DEBUG
    printf("# discr = %i, without normalization\n", discr);
#endif

    INT ret_code = test_aux(ebounds, discr, 0);
    CHECK_RETCODE(ret_code, leave_fun);

#ifdef DEBUG
    printf("\n# discr = %i, with normalization\n", discr);
#endif

    ret_code = test_aux(ebounds, discr, 1);
    CHECK_RETCODE(ret_code, leave_fun);

#ifdef DEBUG
    printf("\n", discr);
#endif


leave_fun:
    return ret_code;
}

int main()
{
    REAL ebounds[3] = {2e-6, 1e-5, 2e-5};

    INT ret_code = test(ebounds, kdv_discretization_BO_VANILLA);
    CHECK_RETCODE(ret_code, on_error);

    ret_code = test(ebounds, kdv_discretization_BO);
    CHECK_RETCODE(ret_code, on_error);

    ebounds[0] = 8e-11;
    ebounds[1] = 2e-5;
    ebounds[2] = 7e-10;
    ret_code = test(ebounds, kdv_discretization_CF4_2);
    CHECK_RETCODE(ret_code, on_error);

    ret_code = test(ebounds, kdv_discretization_CF5_3);
    CHECK_RETCODE(ret_code, on_error);

    ret_code = test(ebounds, kdv_discretization_CF6_4);
    CHECK_RETCODE(ret_code, on_error);

    ebounds[0] = 3e-10;
    ebounds[1] = 2e-5;
    ebounds[2] = 3e-9;
    ret_code = test(ebounds, kdv_discretization_ES4);
    CHECK_RETCODE(ret_code, on_error);

    ebounds[0] = 2e-10;
    ebounds[1] = 2e-5;
    ebounds[2] = 2e-9;
    ret_code = test(ebounds, kdv_discretization_TES4);
    CHECK_RETCODE(ret_code, on_error);

    return EXIT_SUCCESS;

on_error:
    return EXIT_FAILURE;
}

