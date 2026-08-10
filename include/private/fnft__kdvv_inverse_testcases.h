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

#ifndef FNFT__KDVV_INVERSE_TESTCASES_H
#define FNFT__KDVV_INVERSE_TESTCASES_H

#include "fnft_kdvv_inverse.h"
#include "fnft_kdvv.h"


/**
 * List of currently implemented test cases for the inverse KdV.
 *
 *  fnft__inverse_kdvv_testcases_5_bound_states: potential with 5 solitions
 *  fnft__inverse_kdvv_testcases_19_bound_states: potential with 19 solitions
 *  fnft__inverse_kdvv_testcases_8_bound_states_asym: potential with 8 solitions
 *      and assymetric window
 *
 * @ingroup kdv
 */
typedef enum {
    fnft__inverse_kdvv_testcases_5_bound_states,
    fnft__inverse_kdvv_testcases_19_bound_states,
    fnft__inverse_kdvv_testcases_8_bound_states_asym,
} fnft__inverse_kdvv_testcases_t;

 
    
/**
 * @brief Routine to print the continuous spectrum and discrete spectrum, consisting of 
 * bound states and norming constants, as a result of \link fnft_kdvv \endlink.
 *
 * This routine is used by the tests for \link fnft_kdvv_inverse \endlink.
 *
 * @param[in] bound_states
 * @param[in] normconsts 
 * @param[in] contspec
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum.
 * @param[in] M Number of points at which the continuous spectrum is computed.
 * @param[in] D Number of samples of the potential.
 * @param[in] K Number of bound states (same than number of norming constants)
 *
 * @ingroup kdv
 */      
void fnft__kdvv_print_spectrum(  FNFT_COMPLEX const * const bound_states,
                            FNFT_COMPLEX const * const normconsts,
                            FNFT_COMPLEX const * const contspec,
                            FNFT_REAL const * const XI,
                            const UINT M,
                            const UINT D,
                            const UINT K);

/**
 * @brief Routine to run tests for \link fnft_kdvv_inverse \endlink.
 *
 * This routine is used by the tests for \link fnft_kdvv_inverse \endlink. It runs 
 * the specified test case tc with the specified number of samples D and the
 * options opts, and tests if several errors stay below the provided error
 * bounds in error_bounds.
 *
 * @param[in] tc Type of test case.
 * @param[in] D Number of samples.
 * @param[in] error_bounds Real valued array with 4 elements corresponding to various
 * error bounds.
 * @param[in] opts options for the tests.
 *  Note: Has not yet been implemented! If not NULL pointer is handed over, error is
 *      returned!
 * @return If all errors stay below bounds the routine
 * \link FNFT_SUCCESS \endlink. Otherwise, it returns an error code
 * (normally, \link FNFT_EC_TEST_FAILED \endlink).
 *
 * @ingroup kdv
 */
FNFT_INT fnft__inverse_kdvv_testcases_test_fnft( fnft__inverse_kdvv_testcases_t tc, 
                                            UINT D,
                                            const FNFT_REAL error_bounds[4], 
                                            void * const opts);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define inverse_kdvv_testcases_5_bound_states fnft__inverse_kdvv_testcases_5_bound_states
#define inverse_kdvv_testcases_19_bound_states fnft__inverse_kdvv_testcases_19_bound_states
#define inverse_kdvv_testcases_8_bound_states_asym fnft__inverse_kdvv_testcases_8_bound_states_asym
#define inverse_kdvv_testcases_t fnft__inverse_kdvv_testcases_t
#define inverse_kdvv_testcases(...) fnft__inverse_kdvv_testcases(__VA_ARGS__)
#define inverse_kdvv_testcases_test_fnft(...) fnft__inverse_kdvv_testcases_test_fnft(__VA_ARGS__)
#define kdvv_testcases_get_spectrum_of_inverse(...) fnft__kdvv_inverse_testcases_get_spectrum_of_inverse(__VA_ARGS__)
#define kdvv_print_spectrum(...) fnft__kdvv_print_spectrum(__VA_ARGS__)
#endif

#endif
