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
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 1024;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING2;
    REAL error_bounds[6] = {
        7.7e-3,     // reflection coefficient
        3.8e-3,     // a
        1.2e-3,     // b
        1.2e-3,     // bound states
        2.2e-14,    // norming constants
        3.3e-2      // residues
    };
    opts = fnft_nsev_default_opts();
    opts.bound_state_localization = nsev_bsloc_NEWTON;
    opts.discretization = nse_discretization_TES4;
    
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = nsev_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nsev_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    // Check for 4th-order error decay (error_bounds[4] corresponding
    // to the norming constants stays as it is already close to machine precision,
    // error_bounds[5] corresponding to the residues only decays with 2nd-order)
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds[i] /= 16.0;
    error_bounds[4] *= 16.0;
    error_bounds[5] *= 4.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    
    
    D = 1024;
    REAL error_bounds_RE[6] = {
        1.1e-3,     // reflection coefficient
        1.6e-4,     // a
        5.6e-5,     // b
        5.5e-5,     // bound states
        2.2e-14,    // norming constants
        3.9e-2      // residues
    };
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    // Check for at least 5th-order error decay (error_bounds_RE[4] corresponding
    // to the norming constants stays as it is already close to machine precision, 
    // error_bounds_RE[5] corresponding to the residues only decays with 2nd-order)
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= 32.0;
    error_bounds_RE[4] *= 32.0;
    error_bounds_RE[5] *= 8.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    leave_fun:
        if (ret_code != SUCCESS)
            return EXIT_FAILURE;
        else
            return EXIT_SUCCESS;
}

