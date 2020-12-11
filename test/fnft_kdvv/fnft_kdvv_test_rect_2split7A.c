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
 * Peter J Prins (TU Delft) 2017-2018.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdvv_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code;
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    const kdvv_testcases_t tc = kdvv_testcases_RECT;
    opts.discretization = kdv_discretization_2SPLIT7A;
    UINT D = 4;
    REAL eb[6] = {  // error bounds
        2.64e-05,   // continuous spectrum
        FNFT_INF,   // a(xi)
        FNFT_INF,   // b(xi)
        FNFT_INF,   // bound states
        FNFT_INF,   // norming constants
        FNFT_INF    // residues
    };
    
    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    //CHECK_RETCODE(ret_code, leave_fun); // Uncomment if more tests are added
    
//leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

