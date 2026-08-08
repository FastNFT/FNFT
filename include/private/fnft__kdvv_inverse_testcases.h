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
 * @struct fnft_kdvv_params
 * @brief Stores necessary parameter to use \link fnft_kdvv_inverse \endlink.
 * @ingroup fnft
 * @ingroup data_types
 * 
 * @var D
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var T
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var K
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var bound_states
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var normconsts
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var M
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var XI
 *  see \link fnft_kdvv_inverse \endlink
 * 
 * @var contspec
 *  see \link fnft_kdvv_inverse \endlink
 */
typedef struct {
    UINT D;
    REAL T[2];
    UINT K;
    COMPLEX * bound_states;
    COMPLEX * normconsts;
    UINT M;
    REAL * XI;
    COMPLEX * contspec;
} fnft_kdvv_params;


typedef enum {
    fnft__inverse_kdvv_testcases_5_bound_states,
    fnft__inverse_kdvv_testcases_19_bound_states,
    fnft__inverse_kdvv_testcases_8_bound_states_asym,
} fnft__inverse_kdvv_testcases_t;

FNFT_INT fnft__inverse_kdvv_testcases_test_fnft( fnft__inverse_kdvv_testcases_t tc, 
                                            UINT D,
                                            const FNFT_REAL * error_bounds, 
                                            void * const opts);

/**
 * @brief Routine to run tests for \link fnft_kdvv_inverse \endlink.
 *
 * This routine is used by the tests for \link fnft_kdvv_inverse \endlink.
 *
 * @param[in] params_i \link fnft_kdvv_params \endlink
 * @param[in] err_bnd_bound_states 
 * @param[in] err_bnd_spurious_bound_states 
 * @param[in] err_bnd_normconst
 * @param[in] err_bnd_contspec
 * @return If all errors stay below bounds the routine
 * \link FNFT_SUCCESS \endlink. Otherwise, it returns an error code
 * (normally, \link FNFT_EC_TEST_FAILED \endlink).
 *
 * @ingroup kdv
 */
FNFT_INT fnft__kdvv_inverse_testcases_get_spectrum_of_inverse(
    const fnft_kdvv_params params_i,
    const FNFT_REAL err_bnd_bound_states,
    const FNFT_REAL err_bnd_spurious_bound_states,
    const FNFT_REAL err_bnd_normconst,
    const FNFT_REAL err_bnd_contspec);
 
    
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
 * @return void
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
