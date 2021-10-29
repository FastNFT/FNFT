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
* Lianne de Vries (TU Delft student) 2021
* Sander Wahls (TU Delft) 2021
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakovv_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code;
    fnft_manakovv_opts_t opts;
    UINT D = 512;
    const manakovv_testcases_t tc = manakovv_testcases_SECH_FOCUSING;
    REAL error_bounds[5] = { 
        2.2e-4,     // reflection coefficient 1
        2.2e-4,     // reflection coefficient 2
        3.2e-4,     // a
        1.7e-4,     // b1
        1.7e-4     // b2
    };

    opts = fnft_manakovv_default_opts();
    opts.discretization = manakov_discretization_2SPLIT4A;

    ret_code = manakovv_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

