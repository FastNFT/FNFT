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
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code;
    fnft_nsev_opts_t opts;
    UINT D = 2048;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
    REAL error_bounds[6] = {
        FNFT_INF,     // reflection coefficient
        FNFT_INF,     // a
        FNFT_INF,     // b
        4.01e-2,     // bound states
        FNFT_INF,      // norming constants
        FNFT_INF      // residues 
    };

    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_2SPLIT2A;
    opts.Dsub = D/4;
    opts.niter = 0; // no refinement

    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check for quadratic error decay
    opts.Dsub = D/2;
    error_bounds[3] /= 4.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	return EXIT_SUCCESS;
}
