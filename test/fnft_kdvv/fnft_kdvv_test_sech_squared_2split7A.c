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
    const kdvv_testcases_t tc = kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH;
    opts.discretization = kdv_discretization_2SPLIT7A_VANILLA;
    UINT D = 256;
    REAL eb[6] = {  // error bounds
        5.0e-4,     // continuous spectrum
        2.1e-4,     // a(xi)
        3.1e-6,     // b(xi)
        FNFT_INF,   // bound states
        FNFT_INF,   // norming constants
        FNFT_INF    // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Limited test, because it takes long
//    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
//    CHECK_RETCODE(ret_code, leave_fun);
//
//    ret_code = kdvv_testcases_test_fnft(tc, D+1, eb, &opts);
//    CHECK_RETCODE(ret_code, leave_fun);
//
//    ret_code = kdvv_testcases_test_fnft(tc, D-1, eb, &opts);
//    CHECK_RETCODE(ret_code, leave_fun);
//
//    // check for quadratic error decay
//    for (UINT n=0; n<1; n++){
//    D *= 2;
//        for (UINT i=0; i<6; i++)
//             eb[i] /= 4.0;
//        ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
//        CHECK_RETCODE(ret_code, leave_fun);
//    }

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

