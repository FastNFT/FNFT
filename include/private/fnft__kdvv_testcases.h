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
* Sander Wahls (TU Delft) 2018.
*/

/**
 * @file fnft__kdvv_testcases.h
 * @brief Provides test cases for the tests of \link fnft_kdvv \endlink.
 *
 * @ingroup kdv
 */

#ifndef FNFT__KDVV_TESTCASES_H
#define FNFT__KDVV_TESTCASES_H

#include "fnft_kdvv.h" // for fnft_kdvv_opts_t

/**
 * List of currently implemented test cases for the KdV with vanishing
 * boundary conditions.
 *
 * fnft__kdvv_testcases_SECH - A squared sech potential.\n
 * fnft__kdvv_testcases_RECT - A rectangular potential, amplitude=1.\n
 * fnft__kdvv_testcases_NEGATIVE_RECT - A rectangular potential, amplitude=-1.\n
 *
 * @ingroup kdv
 */
typedef enum {
    fnft__kdvv_testcases_SECH_SQUARED,
    fnft__kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH,
    fnft__kdvv_testcases_NEGATIVE_RECT,
} fnft__kdvv_testcases_t;

/**
 * @brief Routine to run tests for \link fnft_kdvv \endlink.
 *
 * This routine is used by the tests for \link fnft_kdvv \endlink. It runs the
 * specified test case tc with the specified number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in eb.
 *
 * @param[in] tc Type of test case.
 * @param[in] D Number of samples.
 * @param[in] eb Real valued array with 6 elements corresponding to various
 * error bounds.
 * @param[in] opts \link fnft_kdvv_opts_t \endlink options for the tests.
 * @return If all errors stay below bounds the routine
 * \link FNFT_SUCCESS \endlink. Otherwise, it returns an error code
 * (normally, \link FNFT_EC_TEST_FAILED \endlink).
 *
 * @ingroup kdv
 */
FNFT_INT fnft__kdvv_testcases_test_fnft(fnft__kdvv_testcases_t tc, FNFT_UINT D,
    const FNFT_REAL eb[6], fnft_kdvv_opts_t * const opts); 

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdvv_testcases_SECH_SQUARED fnft__kdvv_testcases_SECH_SQUARED
#define kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH fnft__kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH
#define kdvv_testcases_NEGATIVE_RECT fnft__kdvv_testcases_NEGATIVE_RECT
#define kdvv_testcases_t fnft__kdvv_testcases_t
#define kdvv_testcases(...) fnft__kdvv_testcases(__VA_ARGS__)
#define kdvv_testcases_test_fnft(...) fnft__kdvv_testcases_test_fnft(__VA_ARGS__)
#endif

#endif
