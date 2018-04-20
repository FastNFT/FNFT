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

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 126;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
    // The error bounds do not matter. This code used to segfault due to a bug
    // in line 426 of fnft_nsev.c ("<=" should have been ">=").
    REAL error_bounds[6] = { 
        INFINITY,     // reflection coefficient
        INFINITY,     // a
        INFINITY,     // b
        INFINITY,     // bound states
        INFINITY,      // norming constants
        INFINITY      // residues
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_2SPLIT4B;

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

