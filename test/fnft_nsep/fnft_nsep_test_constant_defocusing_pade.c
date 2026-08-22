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

#include "fnft__nsep_testcases.h"
#include "fnft__errwarn.h"

static INT run_case(const nse_discretization_t discretization,
        const UINT pade_degree, const INT normalization_flag)
{
    const fnft__nsep_testcases_t testcase =
            nsep_testcases_CONSTANT_DEFOCUSING;
    REAL error_bounds[3] = {
        2.0e-3,
        2.0e-3,
        0.0
    };
    fnft_nsep_opts_t opts = fnft_nsep_default_opts();
    const UINT D = discretization == nse_discretization_FES6_PADE ? 16 : 64;

    opts.discretization = discretization;
    opts.pade_degree = pade_degree;
    opts.localization = fnft_nsep_loc_GRIDSEARCH;
    opts.filtering = fnft_nsep_filt_MANUAL;
    opts.bounding_box[0] = -10.0;
    opts.bounding_box[1] = 10.0;
    opts.bounding_box[2] = -10.0;
    opts.bounding_box[3] = 10.0;
    opts.normalization_flag = normalization_flag;
    return nsep_testcases_test_fnft(testcase, D, error_bounds, &opts);
}

INT main(void)
{
    INT ret_code, normalization_flag;

    for (normalization_flag = 0; normalization_flag <= 1;
            normalization_flag++) {
        ret_code = run_case(nse_discretization_FES4_PADE, 2,
                normalization_flag);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = run_case(nse_discretization_FES6_PADE, 3,
                normalization_flag);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = run_case(nse_discretization_FES6_PADE, 4,
                normalization_flag);
        CHECK_RETCODE(ret_code, leave_fun);
        ret_code = run_case(nse_discretization_FES8_PADE, 3,
                normalization_flag);
        if (ret_code == SUCCESS) {
            ret_code = E_TEST_FAILED;
            goto leave_fun;
        }
        ret_code = SUCCESS;
    }

leave_fun:
    return ret_code == SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
