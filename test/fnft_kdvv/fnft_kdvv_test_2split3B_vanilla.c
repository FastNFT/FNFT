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
* Peter J Prins (TU Delft) 2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdvv_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code;
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    opts.grid_spacing = 0.01;
    opts.discretization = kdv_discretization_2SPLIT3B_VANILLA;

    // Test staircase potential
    kdvv_testcases_t tc = kdvv_testcases_NEGATIVE_RECT;
    UINT D = 3;
    REAL eb_rect[6] = {  // error bounds
        1.2e-1,     // continuous spectrum
        3.5e-2,     // a(xi)
        4.0e-2,     // b(xi)
        0,          // bound states
        0,          // norming constants
        0           // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb_rect, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Test smooth potential
    tc = kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH;
    D = 1024;
    REAL eb[6] = {  // error bounds
        1.4e-3,     // continuous spectrum
        1.7e-4,     // a(xi)
        3.5e-4,     // b(xi)
        7.8e-5,     // bound states
        1.3e-3,     // norming constants
        1.3e-3      // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

