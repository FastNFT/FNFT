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
* Sander Wahls (TU Delft) 2017.
*/

#ifndef FNFT__NSEP_TESTCASES_H
#define FNFT__NSEP_TESTCASES_H

#include "fnft_nsep.h"
#include "fnft__nse_discretization.h"

/**
 * List of currently implemented test cases for the NSE with periodic
 * boundary conditions.
 * @ingroup nse
 * fnft__nsep_testcases_PLANE_WAVE_FOCUSING - Test for focusing NSE (kappa = +1) with a plane wave signal.\n
 * fnft__nsep_testcases_CONSTANT_DEFOCUSING - Test for defocusing NSE (kappa = -1) with a constant signal.\n
 */
typedef enum {
    fnft__nsep_testcases_PLANE_WAVE_FOCUSING,
    fnft__nsep_testcases_CONSTANT_DEFOCUSING
} fnft__nsep_testcases_t;

/**
 * @brief This routine is used by the tests for \link fnft_nsep \endlink. \n
 * @ingroup nse
 * It runs the specified test case tc with the specificed number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in error_bounds. 
 * @param[in] tc Type of test case.
 * @param[in] D Number of samples.
 * @param[in] error_bounds Real valued array with 3 elements corresponding to various error bounds.
 * @param[in] opts \link fnft_nsep_opts_t \endlink options for the tests. 
 * @return If all errors stay below bounds the routine \link FNFT_SUCCESS \endlink. Otherwise,
 * it returns an error code (normally, \link FNFT_EC_TEST_FAILED \endlink).
 */
FNFT_INT fnft__nsep_testcases_test_fnft(fnft__nsep_testcases_t tc, FNFT_UINT D,
    FNFT_REAL error_bounds[3], fnft_nsep_opts_t * const opts);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nsep_testcases_PLANE_WAVE_FOCUSING fnft__nsep_testcases_PLANE_WAVE_FOCUSING
#define nsep_testcases_CONSTANT_DEFOCUSING fnft__nsep_testcases_CONSTANT_DEFOCUSING
#define nsep_testcases_t fnft__nsep_testcases_t
#define nsep_testcases(...) fnft__nsep_testcases(__VA_ARGS__)
#define nsep_testcases_test_fnft(...) fnft__nsep_testcases_test_fnft(__VA_ARGS__)
#endif

#endif
