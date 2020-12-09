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
    opts.discretization = kdv_discretization_CF6_4;
    UINT D = 256;
    REAL eb[6] = {  // error bounds
        9.7e-8,     // continuous spectrum
        2.6e-9,     // a(xi)
        1.1e-8,     // b(xi)
        1.5e-9,     // bound states
        9.4e-9,     // norming constants
        2.4e-9      // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D+1, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D-1, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // check for 6th order error decay
    for (UINT n=0; n<2; n++){
        D *= 2;
        for (UINT i=0; i<5; i++)
            eb[i] /= 64.0;
        eb[5] /= 8.0; // Only 3rd order decay for residues
        ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Test Richardson extrapolation
    opts.richardson_extrapolation_flag = 1;
    D = 256;
    REAL eb_RE[6] = {  // error bounds
        2.4e-10,    // continuous spectrum
        8.1e-11,    // a(xi)
        1.9e-10,    // b(xi)
        1.4e-11,    // bound states
        9.4e-9,     // norming constants
        4.1e-9      // residues
    };

    ret_code = kdvv_testcases_test_fnft(tc, D, eb_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D+1, eb_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdvv_testcases_test_fnft(tc, D-1, eb_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // check for 7th order error decay
    for (UINT n=0; n<1; n++){
        D *= 2;
        for (UINT i=0; i<4; i++)
            eb_RE[i] /= 128.0;
        eb_RE[4] /= 64.0; // Only 6th order error decay for norming constants
        eb_RE[5] /= 16.0; // Only 4th order error decay for residues
        ret_code = kdvv_testcases_test_fnft(tc, D, eb_RE, &opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}

