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
 * @ingroup nse
 */

#ifndef FNFT__MANAKOV_TESTCASES_H
#define FNFT__MANAKOV_TESTCASES_H

#include "fnft_manakov.h" // for fnft_manakov_opts_t

/**
 * List of currently implemented test cases for the NSE with vanishing
 * boundary conditions.\n \n
 *
 * fnft__nsev_testcases_SECH_FOCUSING : Test for focusing NSE (kappa = +1) with a sech potential.\n \n
 * fnft__nsev_testcases_SECH_DEFOCUSING : Test for defocusing NSE (kappa = -1) with a sech potential.\n \n
 * fnft__nsev_testcases_TRUNCATED_SOLITON : Test for focusing NSE (kappa = +1) with a truncated solitonic potential.\n \n
 * fnft__nsev_testcases_SECH_FOCUSING2 : Test for focusing NSE (kappa = +1) with a sech potential with five bound states.\n \n
 * fnft__nsev_testcases_SECH_FOCUSING_CONTSPEC : Test for focusing NSE (kappa = +1) with a sech potential. Only checking continuous spectrum.\n \n
 *
 * @ingroup nse
 */
typedef enum {
    fnft__manakov_testcases_SECH_FOCUSING, 
    fnft__manakov_testcases_SECH_DEFOCUSING,
    fnft__manakov_testcases_TRUNCATED_SOLITON,
    fnft__manakov_testcases_SECH_FOCUSING2,
    fnft__manakov_testcases_SECH_FOCUSING_CONTSPEC
} fnft__manakov_testcases_t;

/**
 * @brief Routine to run tests for \link fnft_nsev \endlink.\n
 * @ingroup nse
 *
 * This routine is used by the tests for \link fnft_nsev \endlink.
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
FNFT_INT fnft__manakov_testcases_test_fnft(fnft__manakov_testcases_t tc, FNFT_UINT D,
	const FNFT_REAL eb[7], fnft_manakov_opts_t * const opts);

INT fnft__manakov_testcases(fnft__manakov_testcases_t tc, const FNFT_UINT D,
	FNFT_COMPLEX** const q1_ptr, FNFT_COMPLEX** const q2_ptr,
	FNFT_REAL* const T,
	FNFT_UINT* const M_ptr, FNFT_COMPLEX** const contspec_ptr,
	FNFT_COMPLEX** const ab_ptr,
	FNFT_REAL* const XI, FNFT_UINT* const K_ptr,
	FNFT_COMPLEX** const bound_states_ptr,
	FNFT_COMPLEX** const normconsts_ptr,
	FNFT_COMPLEX** residues_ptr, FNFT_INT* const kappa_ptr);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_testcases_t fnft__manakov_testcases_t
#define manakov_testcases_SECH_FOCUSING fnft__manakov_testcases_SECH_FOCUSING
#define manakov_testcases_SECH_DEFOCUSING fnft__manakov_testcases_SECH_DEFOCUSING
#define manakov_testcases_TRUNCATED_SOLITON fnft__manakov_testcases_TRUNCATED_SOLITON
#define manakov_testcases_SECH_FOCUSING2 fnft__manakov_testcases_SECH_FOCUSING2
#define manakov_testcases_SECH_FOCUSING_CONTSPEC fnft__manakov_testcases_SECH_FOCUSING_CONTSPEC
#define manakov_testcases(...) fnft__manakov_testcases(__VA_ARGS__)
#define manakov_testcases_test_fnft(...) fnft__manakov_testcases_test_fnft(__VA_ARGS__)

#endif

#endif
