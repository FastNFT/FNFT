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
* Lianne de Vries (TU Delft) 2021
*/
#ifndef FNFT__MANAKOV_DISCRETIZATION_H
#define FNFT__MANAKOV_DISCRETIZATION_H

#include "fnft_manakov_discretization_t.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"


/**
 * @brief This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 * 
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @returns polynomial degree, or 0 for unknown discretizations
 * 
 * @ingroup manakov
 */
FNFT_UINT fnft__manakov_discretization_degree(fnft_manakov_discretization_t
        discretization);

/**
 * @brief This routine returns the scaling for effective number of samples based on the discretization.
 * Returns 0 for unknown discretization
 * 
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @returns upsampling factor. Factor with which to multiply the number of samples to get the effective number of samples for the method
 * 
 * @ingroup manakov
 */
FNFT_UINT fnft__manakov_discretization_upsampling_factor(fnft_manakov_discretization_t discretization);

/**
 * @brief Returns the order of the chosen discretization method. 0 for unknown discretizations
 * 
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @returns order of the chosen discretization method
 * 
 * @ingroup manakov
 */
FNFT_UINT fnft__manakov_discretization_method_order(fnft_manakov_discretization_t discretization);

/**
 * @brief This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization.
 * 
 * @param[in] n Length of array vals (number of \f$\lambda\f$ 's)
 * @param[in] eps_t timestep size
 * @param[in, out] vals Pointer to first element of array of size n. Upon entry, this array contains the continuous-time
 * values \f$\lambda\f$. Upon exit, the array contains the corresponding discrete-time values \f$z\f$. 
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * 
 * @ingroup manakov
 */
FNFT_INT fnft__manakov_discretization_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t,
	FNFT_COMPLEX* const vals, fnft_manakov_discretization_t discretization);

/**
 * @brief This routine returns the phase factor for reflection coefficient (rho).
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 * 
 * @param[in] eps_t timestep size
 * @param[in] T1 time at the last timestep
 * @param[in, out] phase_factor_rho pointer to real-valued variable where the computed phase factor is will be stored
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @return \link FNFT_SUCCESS \endlink or the FNFT__E_INVALID_ARGUMENT(discretization) code
 * 
 * @ingroup manakov
 */
FNFT_INT fnft__manakov_discretization_phase_factor_rho(const FNFT_REAL eps_t, const FNFT_REAL T1,
	FNFT_REAL* const phase_factor_rho, fnft_manakov_discretization_t discretization);

/**
 * @brief This routine returns the phase factor for the NFT coefficient a.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 * 
 * @param[in] eps_t timestep size
 * @param[in] D effective number of samples
 * @param[in] T real array containing the time t of the first and last sample T=[T[0] T[1]]
 * @param[in, out] phase_factor_a pointer to real-valued variable where the computed phase factor is will be stored
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @return \link FNFT_SUCCESS \endlink or the FNFT__E_INVALID_ARGUMENT(discretization) code
 * 
 * @ingroup manakov
 */
FNFT_INT fnft__manakov_discretization_phase_factor_a(const FNFT_REAL eps_t, const FNFT_UINT D, FNFT_REAL const* const T,
	FNFT_REAL* const phase_factor_a, fnft_manakov_discretization_t discretization);

/**
 * @brief This routine returns the phase factor for the NFT coefficient b1, b2.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 * 
 * @param[in] eps_t timestep size
 * @param[in] D effective number of samples
 * @param[in] T real array containing the time t of the first and last sample T=[T[0] T[1]]
 * @param[in, out] phase_factor_b pointer to real-valued variable where the computed phase factor is will be stored
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @return \link FNFT_SUCCESS \endlink or the FNFT__E_INVALID_ARGUMENT(discretization) code
 * 
 * @ingroup manakov
 */
FNFT_INT fnft__manakov_discretization_phase_factor_b(const FNFT_REAL eps_t, const FNFT_UINT D, FNFT_REAL const* const T,
	FNFT_REAL* const phase_factor_b, fnft_manakov_discretization_t discretization);

/**
 * @brief This routine computes the downsampled and interpolated version of the samples if required
 * 
 * For richardson extrapolation downsampling is needed. This simply means only using a certain number from the original samples.
 * Not all choices for the number od samples are valid as at this stage the samples should still be equidistant and no interpolation
 * is carried out. If the chosen number of samples for downsampling is not valid, the routine chooses a number as close as possible to the desired amount of samples.
 * Certain discretization types from \link fnft_manakov_discretization_t \endlink require
 * non-equidistant samples. In that case, band-limited Fourier interpolation is used after the optional downsampling to get the required samples.
 * Both the downsampling and interpolation step are combined in this function, resulting in a preprocessed array
 * of samples to be used by the other functions in \link fnft_manakov \endlink
 * 
 * @param[in] D the original number of samples
 * @param[in] q1 Pointer to complex array of size D containg the samples of the first element of the potential q1
 * @param[in] q2 Pointer to complex array of size D containg the samples of the second element of the potential q2
 * @param[in] eps_t timestep size
 * @param[in] kappa dispersion constant. kappa=-1 for normal dispersion and kappa=+1 for anomalous dispersion
 * @param[in, out] Dsub_ptr Pointer to the desired number of samples after subsampling. Upon exit, this contains the
 * actual number of samples after subsampling which ensures that the samples stay equidistant
 * @param[out] q1_preprocessed_ptr Pointer to a complex array of size Dsub*upsampling_factor containing the optionally downsampled and
 * interpolated (if the discretization required interpolation) samples
 * @param[out] first_last_index Pointer to an integer array containing the index of the first and last index from the original samples
 * used in the preprocessed arrays
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * 
 * @ingroup manakov
 */
FNFT_UINT fnft__manakov_discretization_preprocess_signal(const FNFT_UINT D, FNFT_COMPLEX const* const q1,
		FNFT_COMPLEX const* const q2, FNFT_REAL const eps_t, const FNFT_INT kappa,
	FNFT_UINT* const Dsub_ptr, FNFT_COMPLEX** q1_preprocessed_ptr, FNFT_COMPLEX** q2_preprocessed_ptr,
	FNFT_UINT* const first_last_index, manakov_discretization_t discretization);

/**
 * @brief This routine returns the boundary coefficient based on the discretization.
 * 
 * @param[in] discretization discretization Should be of type \link fnft_manakov_discretization_t \endlink
 * @returns the boundary coefficient
 * 
 * @ingroup manakov
 */
FNFT_REAL fnft__manakov_discretization_boundary_coeff(fnft_manakov_discretization_t discretization);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_discretization_degree(...) fnft__manakov_discretization_degree(__VA_ARGS__)
#define manakov_discretization_upsampling_factor(...) fnft__manakov_discretization_upsampling_factor(__VA_ARGS__)
#define manakov_discretization_method_order(...) fnft__manakov_discretization_method_order(__VA_ARGS__)
#define manakov_discretization_lambda_to_z(...) fnft__manakov_discretization_lambda_to_z(__VA_ARGS__)
#define manakov_discretization_phase_factor_rho(...) fnft__manakov_discretization_phase_factor_rho(__VA_ARGS__)
#define manakov_discretization_phase_factor_a(...) fnft__manakov_discretization_phase_factor_a(__VA_ARGS__)
#define manakov_discretization_phase_factor_b(...) fnft__manakov_discretization_phase_factor_b(__VA_ARGS__)
#define manakov_discretization_boundary_coeff(...) fnft__manakov_discretization_boundary_coeff(__VA_ARGS__)
#define manakov_discretization_preprocess_signal(...) fnft__manakov_discretization_preprocess_signal(__VA_ARGS__)
#endif

#endif