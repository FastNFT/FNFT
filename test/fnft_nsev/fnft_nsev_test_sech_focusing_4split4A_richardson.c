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
* Sander Wahls (TU Delft) 2017-2018, 2021.
* Shrinivas Chimmalgi (TU Delft) 2017-2020.
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

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_4SPLIT4A;
   
    // Check for at least 5th-order error decay on resulting from application
    // of Richardson extrapolation to 4th-order method.(error_bounds[4] stays 
    // as it is because it is already close to machine precision)
    D = 512;
    REAL error_bounds_RE[6] = {
        4.4e-8, // reflection coefficient
        5.6e-7, // a
        1.1e-7, // b,
        3.1e-9, // bound states
        5e-14,  // norming constants
        3.4e-9 // residues
    };
    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    UINT DN = 800;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= POW((REAL)DN/(REAL)D,5);
    error_bounds_RE[4] *= POW((REAL)DN/(REAL)D,5);
    ret_code = nsev_testcases_test_fnft(tc, DN, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

