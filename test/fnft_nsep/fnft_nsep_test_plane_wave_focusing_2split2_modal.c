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

#include "fnft__nsep_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    const fnft__nsep_testcases_t tc = nsep_testcases_PLANE_WAVE_FOCUSING;
    UINT D = 1024;
    REAL error_bounds[3] = {
        3.5e-4, // main spectrum
        3.1e-4, // aux spectrum
        0.0     // sheet indices (zero since not yet implemented)
    };
    fnft_nsep_opts_t opts;

    opts = fnft_nsep_default_opts();
    opts.discretization = nse_discretization_2SPLIT2_MODAL;
    opts.localization = fnft_nsep_loc_MIXED;
    opts.filtering = fnft_nsep_filt_MANUAL;
    opts.bounding_box[0] = -10;
    opts.bounding_box[1] = 10;
    opts.bounding_box[2] = -10;
    opts.bounding_box[3] = 10;

    ret_code = nsep_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check for error decay. The poly_roots_fftgridsearch routine used to
    // find the spectrum on real line guarantees only linear error decay
    D *= 2;
    error_bounds[0] /= 2.0;
    error_bounds[1] /= 4.0;
    error_bounds[2] /= 4.0;
    ret_code = nsep_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Repeat tests without real spectrum
    opts.bounding_box[2] = 0.1;
    D = 1024;
    error_bounds[0] = 4.4e-5;
    error_bounds[1] = 4.4e-5;
    error_bounds[2] = 0.0;
    ret_code = nsep_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Error decay should now be quadratic.
    D *= 2;
    for (i=0; i<3; i++)
        error_bounds[i] /= 4.0;
    ret_code = nsep_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
