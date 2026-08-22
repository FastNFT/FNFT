/*
 * This file is part of FNFT.
 *
 * FNFT is free software; you can redistribute it and/or
 * modify it under the terms of the version 2 of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * FNFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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

INT main()
{
    INT ret_code;
    UINT D = 512;
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    REAL error_bounds[6] = {
        6.6e-4, INFINITY, INFINITY, 0.0, 0.0, 0.0
    };
    opts.discretization = nse_discretization_CT4;

    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds,&opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = 0;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds,&opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = 1;

    D *= 2;
    for (UINT i=0; i<6; i++)
        error_bounds[i] /= 16.0;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds,&opts);
    CHECK_RETCODE(ret_code, leave_fun);

    D = 512;
    REAL error_bounds_RE[6] = {
        8.2e-5, INFINITY, INFINITY, 0.0, 0.0, 0.0
    };
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds_RE,&opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = 0;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds_RE,&opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = 1;

    D *= 2;
    for (UINT i=0; i<6; i++)
        error_bounds_RE[i] /= 32.0;
    ret_code = nsev_testcases_test_fnft(nsev_testcases_SECH_DEFOCUSING,
            D,error_bounds_RE,&opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code == SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
