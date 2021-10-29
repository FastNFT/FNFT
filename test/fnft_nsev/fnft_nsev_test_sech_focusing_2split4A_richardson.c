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
 * Sander Wahls (TU Delft) 2017-2018, 2021.
 * Shrinivas Chimmalgi (TU Delft) 2017.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nsev_testcases.h"
#include "fnft__errwarn.h"

INT main()
{
    INT ret_code, i;
    fnft_nsev_opts_t opts;
    UINT D = 4096;
    const nsev_testcases_t tc = nsev_testcases_SECH_FOCUSING;
 
    opts = fnft_nsev_default_opts();
    opts.discretization = nse_discretization_2SPLIT4A;
    
    // Check for Richardson
    REAL error_bounds_RE[6] = {
        2.7e-8,     // reflection coefficient
        6.7e-8,     // a
        2.4e-8,     // b
        1.2e-9,     // bound states
        5e-14,      // norming constants
        4.5e-9      // residues
    };

    opts.richardson_extrapolation_flag = 1;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D *= 2;
    for (i=0; i<6; i++)
        error_bounds_RE[i] /= 16.0;
    error_bounds_RE[4] *= 16.0;
    ret_code = nsev_testcases_test_fnft(tc, D, error_bounds_RE, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
leave_fun:
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

