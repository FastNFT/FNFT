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
 * Sander Wahls (TU Delft) 2018, 2020.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__akns_fscatter.h"
#include "fnft__kdv_finvscatter.h"
#include "fnft__misc.h"

static INT kdv_finvscatter_test(const UINT D,
                                const REAL error_bound,
                                const kdv_discretization_t kdv_discretization)
{
    INT ret_code = SUCCESS;

    COMPLEX *r = NULL, *r_exact = NULL, *q_exact = NULL, *result = NULL;
    REAL eps_t = 0.12;
    UINT deg, i;

    akns_discretization_t akns_discretization;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization,
                                                         &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    const UINT len = akns_fscatter_numel(D, akns_discretization);

    r = malloc(D * sizeof(COMPLEX));
    r_exact = malloc(D * sizeof(COMPLEX));
    q_exact = malloc(D * sizeof(COMPLEX));
    result = malloc(len * sizeof(COMPLEX));
    if (r == NULL || r_exact == NULL || q_exact == NULL || result == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    for (i=0; i<D; i++) {
        r_exact[i] = ((REAL)(i+1)/(D+1)/D)*CEXP(I*(REAL)i/D);
        q_exact[i] = -1.0;
        r[i] = 0.0;
    }

    ret_code = akns_fscatter(D, q_exact, r_exact, eps_t, result, &deg, NULL,
        akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = kdv_finvscatter(deg, result, r, eps_t, kdv_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    misc_print_buf(D, r, "r");
    misc_print_buf(D, r_exact, "r_exact");

    REAL error = misc_rel_err(D, r, r_exact);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(r);
    free(r_exact);
    free(q_exact);
    free(result);
    return ret_code;
}
