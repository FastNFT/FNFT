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
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
* Sander Wahls (KIT) 2023.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    const nsev_testcases_t tc = nsev_testcases_SECH_DEFOCUSING;
    UINT D = 256;
    REAL error_bounds[6] = { 
        3.0e-5,     // reflection coefficient
        INFINITY,   // a
        INFINITY,   // b
        0.0,        // bound states
        0.0,        // norming constants
        0.0         // residues 
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_CF4_3;

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = nsev_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nsev_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check for 4th-order error decay
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds[i] /= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D = 512;
    REAL error_bounds_RE[6] = { 
        2.1e-7,     // reflection coefficient
        INFINITY,   // a
        INFINITY,   // b
        0.0,        // bound states
        0.0,        // norming constants
        0.0         // residues 
    };
    opts.richardson_extrapolation_flag = 1; 
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);    

    // Repeat without normalization
    opts.normalization_flag = !opts.normalization_flag;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    opts.normalization_flag = !opts.normalization_flag;

    // Check for at least 5th-order error decay
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= 32.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);    

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	return EXIT_SUCCESS;
}

