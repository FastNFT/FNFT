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
 * Sander Wahls (TU Delft) 2022.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft_nsev.h"
#include "fnft__nsev_testcases.h"

static INT run_test(UINT D, fnft_nsev_opts_t * opts_ptr, REAL err_bnd)
{
    const REAL error_bounds[6] = {
        INFINITY,
        INFINITY,
        INFINITY,
        err_bnd, // bound states
        INFINITY,
        INFINITY};
    INT ret_code = fnft__nsev_testcases_test_fnft(
            nsev_testcases_SECH_FOCUSING2, D, error_bounds, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);
leave_fun:
    return ret_code;
}

int main()
{
    // We make sure that gprf achieves the same precision
    // as the default method (with appropriately chosen tol).

    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    INT ret_code = run_test(1024, &opts, 3.5e-2);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.bound_state_localization = nsev_bsloc_GRPF;
    opts.tol = 1e-3;
    ret_code = run_test(1024, &opts, 3.5e-2);
    CHECK_RETCODE(ret_code, leave_fun);

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_4SPLIT4B;
    ret_code = run_test(1024, &opts, 8e-5);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.discretization = nse_discretization_CF4_2;
    opts.bound_state_localization = nsev_bsloc_GRPF;
    opts.tol = 1e-5;
    ret_code = run_test(1024, &opts, 8e-5);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
