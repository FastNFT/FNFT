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
* Fabian Fischer (Hiwi KIT) 2026.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft__kdvv_inverse_testcases.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

/* This is a testcase with an assymetric window */

INT main()
{
    INT ret_code = SUCCESS;

    inverse_kdvv_testcases_t testcase = inverse_kdvv_testcases_8_bound_states_asym;
    REAL error_bounds[4] = {
        3.8e-3,         // bound states
        2e-1,           // norming constants
        2.5e-1,         // continuous spectrum
        0.12,           // spurious bound states
    };

    UINT D = 256;

    ret_code = inverse_kdvv_testcases_test_fnft(testcase, D, error_bounds, NULL);
    CHECK_RETCODE(ret_code, leave_fun);


    for (UINT n=0; n<3; n++){
        D *= 2;
        for (UINT i=0; i<4; i++)
            error_bounds[i] /= 4.0;
        ret_code = inverse_kdvv_testcases_test_fnft(testcase, D, error_bounds, NULL);
        CHECK_RETCODE(ret_code, leave_fun);
    }    

leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
