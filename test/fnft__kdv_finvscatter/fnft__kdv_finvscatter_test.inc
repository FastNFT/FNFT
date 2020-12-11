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
 * Peter J Prins (TU Delft) 2020.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__akns_fscatter.h"
#include "fnft__kdv_fscatter.h"
#include "fnft__kdv_finvscatter.h"
#include "fnft__misc.h"

static INT kdv_finvscatter_test(const UINT D,
                                const REAL error_bound,
                                const kdv_discretization_t kdv_discretization)
{
    INT ret_code = SUCCESS;

    COMPLEX *q_exact = NULL, *r_exact = NULL, *r = NULL, *q = NULL, *transfer_matrix = NULL;
    REAL eps_t = 0.12;
    UINT deg, i, vanilla_flag;
    REAL error;

    CHECK_RETCODE(ret_code, leave_fun);
    const UINT len = kdv_fscatter_numel(D, kdv_discretization);

    q_exact = malloc(D * sizeof(COMPLEX));
    r_exact = malloc(D * sizeof(COMPLEX));
    transfer_matrix = malloc(len * sizeof(COMPLEX));
    if (q_exact == NULL || r_exact == NULL || transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    INT kappa = 1;
    ret_code = kdv_discretization_vanilla_flag(&vanilla_flag, kdv_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    if (vanilla_flag) {

        for (i=0; i<D; i++) {
            q_exact[i] = ((REAL)(i+1)/(D+1)/D)*CEXP(I*(REAL)i/D);
            r_exact[i] = -1;
        }
        ret_code = kdv_fscatter(D, q_exact, r_exact, eps_t, kappa, transfer_matrix, &deg, NULL,
                                kdv_discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        q = malloc(D * sizeof(COMPLEX));
        if (q == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        ret_code = kdv_finvscatter(deg, transfer_matrix, q, eps_t,
                                   kdv_discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        error = misc_rel_err(D, q, q_exact);

    } else {

        for (i=0; i<D; i++) {
            q_exact[i] = -1;
            r_exact[i] = ((REAL)(i+1)/(D+1)/D)*CEXP(I*(REAL)i/D);
        }
        ret_code = kdv_fscatter(D, q_exact, r_exact, eps_t, kappa, transfer_matrix, &deg, NULL,
                                kdv_discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        r = malloc(D * sizeof(COMPLEX));
        if (r == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        ret_code = kdv_finvscatter(deg, transfer_matrix, r, eps_t,
                                   kdv_discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        error = misc_rel_err(D, r, r_exact);

    }
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(q_exact);
    free(r_exact);
    free(r);
    free(q);
    free(transfer_matrix);
    return ret_code;
}
