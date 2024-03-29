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
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J Prins (TU Delft) 2020.
 * Sander Wahls (KIT) 2023.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 400;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING2;
    REAL error_bounds[6] = {
        1.7e-4,     // reflection coefficient
        7.9e-5,     // a
        2.4e-5,     // b
        4.3e-5,     // bound states
        5.0e-14,    // norming constants
        5.3e-5      // residues
    };
    
    opts = fnft_nsev_default_opts();
    opts.bound_state_localization = nsev_bsloc_NEWTON;
    opts.discretization = nse_discretization_CF6_4;
    
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = nsev_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nsev_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Check for 6th-order error decay (error_bounds[4] corresponding
    // to the norming constants stays as it is already close to machine precision.
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds[i] /= 64.0;
    error_bounds[4] *= 64.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D = 512;
    REAL error_bounds_RE[6] = {
        7e-7,     // reflection coefficient
        4e-7,     // a
        9.9e-8,     // b
        4.6e-7,     // bound states
        5e-14,      // norming constants
        6.1e-7      // residues
    };
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Repeat without normalization
    opts.normalization_flag = !opts.normalization_flag;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = !opts.normalization_flag;

    // Check for at least 7th-order error decay (error_bounds_RE[4] corresponding
    // to the norming constants stays as it is already close to machine precision)
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= 128.0;
    error_bounds_RE[4] *= 128.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun)
    
    leave_fun:
        if (ret_code != SUCCESS)
            return EXIT_FAILURE;
        else
            return EXIT_SUCCESS;
}

