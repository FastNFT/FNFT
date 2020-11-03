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
 * Sander Wahls (TU Delft) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_nsev_inverse.h"
#include "fnft_nsev.h"
#include "fnft__errwarn.h"
#include "fnft__nse_discretization.h"
#include "fnft__misc.h"
#include <stdio.h>

static INT fnft_nsev_inverse_test(const REAL error_bound,
       fnft_nsev_inverse_opts_t * opts_ptr)
{
#include "fnft_nsev_inverse_test_sech_defocusing_data_2048.inc"
#include "fnft_nsev_inverse_test_sech_defocusing_data_4096.inc"

    if (!(error_bound > 0))
        return E_INVALID_ARGUMENT(error_bound);
    if (opts_ptr == NULL)
        return E_INVALID_ARGUMENT(opts_ptr);

    const INT kappa = -1;
    INT ret_code = SUCCESS;

    COMPLEX * const q = malloc(4096 * sizeof(COMPLEX)) ;
    if (q == NULL) {
       ret_code = E_NOMEM;
       goto leave_fun;
    }

    ret_code = fnft_nsev_inverse(M_2048, contspec_2048, XI_2048, 0, NULL, NULL,
        2048, q, T_2048, kappa, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL error = misc_rel_err(2048, q, q_exact_2048);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    ret_code = fnft_nsev_inverse(M_4096, contspec_4096, XI_4096, 0, NULL, NULL,
                                 4096, q, T_4096, kappa, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);

    error = misc_rel_err(4096, q, q_exact_4096);
#ifdef DEBUG
    printf("error = %g\n", error);
#endif
    if (!(error < error_bound/4)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(q);
    return ret_code;
}
