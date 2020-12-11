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
* Peter J Prins (TU Delft) 2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdvv_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code;
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    const kdvv_testcases_t tc = kdvv_testcases_SECH_SQUARED;
    opts.discretization = kdv_discretization_CF4_3_VANILLA;
    UINT D = 128;
    REAL eb[6] = {  // error bounds
        5.9e-4,     // continuous spectrum
        3.5e-7,     // a(xi)
        4.0e-4,     // b(xi)
        9.7e-7,     // bound states
        1.6e-5,     // norming constants
        1.7e-5      // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D+1, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D-1, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // check for 4th order error decay
    for (UINT n=0; n<5; n++){
        D *= 2;
        for (UINT i=0; i<6; i++)
            eb[i] /= 16.0;
        ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

