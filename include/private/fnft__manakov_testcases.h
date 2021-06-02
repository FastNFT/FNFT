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
* Lianne de Vries (TU Delft) 2021.
*/

/**
 * @file fnft__manakov_testcases.h
 * @brief Provides test cases for the tests of \link fnft_manakov \endlink.
 *
 * @ingroup manakov
 */

#ifndef FNFT__MANAKOV_TESTCASES_H
#define FNFT__MANAKOV_TESTCASES_H

#include "fnft_manakov.h" // for fnft_manakov_opts_t

/**
 * List of currently implemented test cases for the Manakov equation\n \n
 *
 * fnft__manakov_testcases_SECH_FOCUSING : Test for focusing Manakov (kappa = +1) with a sech potential.\n \n
 * fnft__manakov_testcases_SECH_DEFOCUSING : Test for defocusing Manakov (kappa = -1) with a sech potential.\n \n
 * fnft__manakov_testcases_RECTANGLE_FOCUSING : Test for focusing Manakov (kappa = 1) with a rectangle potential.\n
 *
 * @ingroup manakov
 */
typedef enum {
    fnft__manakov_testcases_SECH_FOCUSING, 
    fnft__manakov_testcases_SECH_DEFOCUSING,
    fnft__manakov_testcases_RECTANGLE_FOCUSING
} fnft__manakov_testcases_t;

/**
 * @brief Routine to run tests for \link fnft_manakov \endlink.\n
 * @ingroup manakov
 *
 * This routine is used by the tests for \link fnft_manakov \endlink.
 * It runs the specified test case tc with the specificed number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in eb. 
 * @param[in] tc Type of test case.
 * @param[in] D Number of samples.
 * @param[in] eb Real valued array with 5 elements corresponding to various error bounds.
 * @param[in] opts \link fnft_manakov_opts_t \endlink options for the tests. 
 * @return If all errors stay below bounds the routine \link FNFT_SUCCESS \endlink. Otherwise,
 * it returns an error code (normally, \link FNFT_EC_TEST_FAILED \endlink).
 */
FNFT_INT fnft__manakov_testcases_test_fnft(fnft__manakov_testcases_t tc, FNFT_UINT D,
	const FNFT_REAL eb[5], fnft_manakov_opts_t * const opts);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_testcases_t fnft__manakov_testcases_t
#define manakov_testcases_SECH_FOCUSING fnft__manakov_testcases_SECH_FOCUSING
#define manakov_testcases_SECH_DEFOCUSING fnft__manakov_testcases_SECH_DEFOCUSING
#define manakov_testcases_RECTANGLE_FOCUSING fnft__manakov_testcases_RECTANGLE_FOCUSING
#define manakov_testcases(...) fnft__manakov_testcases(__VA_ARGS__)
#define manakov_testcases_test_fnft(...) fnft__manakov_testcases_test_fnft(__VA_ARGS__)
#endif

#endif
