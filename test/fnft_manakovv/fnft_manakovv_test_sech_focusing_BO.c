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
* Contributor:
* Lianne de Vries (TU Delft) 2021.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakovv_testcases.h"
#include "fnft__errwarn.h"

// TODO: this one does not work now, just used for debugging

INT main()
{
    INT ret_code, i;
    fnft_manakovv_opts_t opts;
    UINT D = 512;
    const manakovv_testcases_t tc = manakovv_testcases_SECH_FOCUSING;
    REAL error_bounds[5] = { 
        2.3e-4,     // reflection coefficient 1
        2.3e-4,     // reflection coefficient 2
        5.7e-4,     // a
        5.4e-5,     // b1
        5.4e-5     // b2
    };

    opts = fnft_manakovv_default_opts();
    opts.discretization = manakov_discretization_BO;


    ret_code = manakovv_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check the case where D is not a power of two. The error bounds have to
    // be tight but not too tight for this to make sense!
    ret_code = manakovv_testcases_test_fnft(tc, D+1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = manakovv_testcases_test_fnft(tc, D-1, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
