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
* Sander Wahls (TU Delft) 2017-2018, 2023.
* Peter J Prins (TU Delft) 2020.
* Sander Wahls (KIT) 2023.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdvv_testcases.h"
#include "fnft__errwarn.h"

INT test(const UINT D, const REAL eb[6], fnft_kdv_discretization_t discretization, const REAL grid_spacing)
{
    INT ret_code;
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    const kdvv_testcases_t tc = kdvv_testcases_SECH_SQUARED;
    opts.discretization = discretization;
    opts.grid_spacing = grid_spacing;

    opts.niter = 0; // no Newton refinement!
    ret_code = kdvv_testcases_test_fnft(tc, D, eb, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}

int main() {

    // The grid search code does not depend on the discretization, so
    // we just test two choices. The error bounds are from the
    // corresponding tests with Newton refinement. This is just to
    // verify that grid search alone can archive the same errors
    // (with much smaller grid spacing, of course).

    UINT D = 256;
    REAL eb_4split4B[6] = {  // error bounds
        3.0e-6,     // continuous spectrum
        7.8e-4,     // a(xi)
        3.6e-6,     // b(xi)
        1.7e-6,     // bound states
        6.5e-6,     // norming constants
        5.9e-6      // residues
    };

    REAL grid_spacing = 0.00001;
    INT ret_code = test(D, eb_4split4B, kdv_discretization_4SPLIT4B, grid_spacing);
    CHECK_RETCODE(ret_code, failure);

    REAL eb_BO[6] = {       // error bounds
        8.3e-4,     // continuous spectrum
        1.6e-4,     // a(xi)
        2.7e-6,     // b(xi)
        1.3e-3,     // bound states
        2.1e-2,     // norming constants
        2.0e-2      // residues
    };
    grid_spacing = 0.00025;
    ret_code = test(D, eb_BO, kdv_discretization_BO, grid_spacing);
    CHECK_RETCODE(ret_code, failure);

    return EXIT_SUCCESS;
failure:
    return EXIT_FAILURE;
}

