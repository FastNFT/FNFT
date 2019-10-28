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
* Shrinivas Chimmalgi (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_slow_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 4096;
    const nsev_slow_testcases_t tc = nsev_slow_testcases_SECH_FOCUSING;
    REAL error_bounds[6] = { 
        3.9e-6,     // reflection coefficient
        6.3e-6,     // a
        2.0e-6,     // b
        1.6e-5,     // bound states
        INFINITY,//5e-14,      // norming constants
        INFINITY//2.1e-6      // residues
    };

    opts = fnft_nsev_default_opts();
    opts.bound_state_localization = nsev_bsloc_NEWTON;
    opts.discretization = nse_discretization_BO;

    ret_code = nsev_slow_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = nsev_slow_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nsev_slow_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
 
    // Check for quadratic error decay (error_bounds[4] stays as it is
    // already close to machine precision)
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds[i] /= 4.0;
    error_bounds[4] *= 4.0;
    ret_code = nsev_slow_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

