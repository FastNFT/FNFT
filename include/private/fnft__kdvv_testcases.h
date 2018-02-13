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
 * @ingroup kdv
 */

#ifndef FNFT__KDVV_TESTCASES_H
#define FNFT__KDVV_TESTCASES_H

#include "fnft_kdvv.h" // for fnft_kdvv_opts_t

/**
 * List of currently implemented test cases for the KdV with vanishing
 * boundary conditions.
 * @ingroup kdv
 */
typedef enum {
    fnft__kdvv_testcases_SECH,
    fnft__kdvv_testcases_RECT
} fnft__kdvv_testcases_t;

/**
 * This routine is used by the tests for \link fnft_kdvv \endlink. It runs the
 * specified test case tc with the specificed number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in eb. If yes, the routine \link FNFT_SUCCESS \endlink. Otherwise,
 * it returns an error code (normally, \link FNFT_EC_TEST_FAILED \endlink).
 * @ingroup nse
 */
FNFT_INT fnft__kdvv_testcases_test_fnft(fnft__kdvv_testcases_t tc, FNFT_UINT D,
    const FNFT_REAL eb[6], fnft_kdvv_opts_t * const opts); 

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdvv_testcases_SECH fnft__kdvv_testcases_SECH
#define kdvv_testcases_RECT fnft__kdvv_testcases_RECT
#define kdvv_testcases_t fnft__kdvv_testcases_t
#define kdvv_testcases(...) fnft__kdvv_testcases(__VA_ARGS__)
#define kdvv_testcases_test_fnft(...) fnft__kdvv_testcases_test_fnft(__VA_ARGS__)
#endif

#endif
