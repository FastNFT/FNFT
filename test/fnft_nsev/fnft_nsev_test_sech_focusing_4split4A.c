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
* Shrinivas Chimmalgi (TU Delft) 2017-2019.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 512;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
    REAL error_bounds[6] = { 
        1.6e-6,//10.0e-8,     // reflection coefficient
        4.3e-6,//2.7e-7,     // a
        1.5e-6,//8.3e-8,     // b
        2.2e-7,//1.4e-8,     // bound states
        5e-15,      // norming constants
        1.1e-6//6.4e-8      // residues 
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_4SPLIT4A;

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = nsev_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nsev_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);


    // Check for fourth-order error decay (error_bounds[4] stays as it is because it is
    // already close to machine precision)
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds[i] /= 16.0;
    error_bounds[4] *= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
        // Check for Richardson
    D /= 2;
    REAL error_bounds_RE[6] = {
        4.4e-8,//8.9e-10,     // reflection coefficient
        5.6e-7,//1.9e-8,     // a
        9.3e-8,//3.6e-9,     // b
        3.1e-9,//4.8e-11,     // bound states
        5e-14,      // norming constants
        3.4e-9//5.3e-11      // residues
    };
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= 16.0;
    error_bounds_RE[4] *= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

