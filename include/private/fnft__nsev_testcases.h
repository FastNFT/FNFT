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

/**
 * @file fnft__nsev_testcases.h
 * @brief Provides test cases for the tests of \link fnft_nsev \endlink.
 * @ingroup nse
 */

#ifndef FNFT__NSEV_TESTCASES_H
#define FNFT__NSEV_TESTCASES_H

#include "fnft_nsev.h" // for fnft_nsev_opts_t
#include "fnft__nse_discretization.h"

/**
 * List of currently implemented test cases for the NSE with vanishing
 * boundary conditions.\n \n
 * @ingroup nse
 * fnft__nsev_testcases_SECH_FOCUSING : Test for focusing NSE (kappa = +1) with a sech potential.\n \n
 * fnft__nsev_testcases_SECH_DEFOCUSING : Test for defocusing NSE (kappa = -1) with a sech potential.\n \n
 * fnft__nsev_testcases_TRUNCATED_SOLITON : Test for focusing NSE (kappa = +1) with a truncated solitonic potential.
 */
typedef enum {
    fnft__nsev_testcases_SECH_FOCUSING, 
    fnft__nsev_testcases_SECH_DEFOCUSING,
    fnft__nsev_testcases_TRUNCATED_SOLITON
} fnft__nsev_testcases_t;

/**
 * @brief Routine to run tests for \link fnft_nsev \endlink.
 * 
 * @ingroup nse
 * It runs the specified test case tc with the specificed number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in eb. 
 * @param[in] tc Type of test case.
 * @param[in] D Number of samples.
 * @param[in] eb Real valued array with 6 elements corresponding to various error bounds.
 * @param[in] opts \link fnft_nsev_opts_t \endlink options for the tests. 
 * @return If all errors stay below bounds the routine \link FNFT_SUCCESS \endlink. Otherwise,
 * it returns an error code (normally, \link FNFT_EC_TEST_FAILED \endlink).
 */
FNFT_INT fnft__nsev_testcases_test_fnft(fnft__nsev_testcases_t tc, FNFT_UINT D,
	const FNFT_REAL eb[6], fnft_nsev_opts_t * const opts);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nsev_testcases_t fnft__nsev_testcases_t
#define nsev_testcases_SECH_FOCUSING fnft__nsev_testcases_SECH_FOCUSING
#define nsev_testcases_SECH_DEFOCUSING fnft__nsev_testcases_SECH_DEFOCUSING
#define nsev_testcases_TRUNCATED_SOLITON fnft__nsev_testcases_TRUNCATED_SOLITON
#define nsev_testcases(...) fnft__nsev_testcases(__VA_ARGS__)
#define nsev_testcases_test_fnft(...) fnft__nsev_testcases_test_fnft(__VA_ARGS__)
#endif

#endif
