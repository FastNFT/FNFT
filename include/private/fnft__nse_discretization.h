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
 * Shrinivas Chimmalgi (TU Delft) 2017.
 */

/**
 * @brief Properties of the discretizations for the nonlinear Schroedinger
 * equation.
 *
 * @file fnft__nse_discretization.h
 * @ingroup nse
 */

#ifndef FNFT__NSE_DISCRETIZATION_H
#define FNFT__NSE_DISCRETIZATION_H

#include "fnft_nse_discretization_t.h"
#include "fnft__akns_discretization.h"
#include "fnft__errwarn.h"



/**
 * @brief This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 *
 * @param[in] nse_discretization The type of discretization to be used. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @returns polynomial degree, or 0 for discretizations not supported by \link fnft__nse_fscatter \endlink.
 *
 * @ingroup nse
 */
FNFT_UINT fnft__nse_discretization_degree(fnft_nse_discretization_t
        nse_discretization);

/**
 * @brief This routine returns the boundary coefficient based on the
 * discretization.
 *
 * The boundary coefficient is the fraction of the step size that a discretized
 * potential extends beyond the last sample. This routine returns this value
 * based on the discretization of type \link fnft_nse_discretization_t \endlink.
 * @param[in] nse_discretization The type of discretization to be used. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @returns the boundary coefficient, or NAN for discretizations not supported
 * by \link fnft__nse_fscatter \endlink.
 *
 * @ingroup nse
 */
FNFT_REAL fnft__nse_discretization_boundary_coeff(fnft_nse_discretization_t nse_discretization);

/**
 * @brief This routine returns akns discretization related to the given
 * nse discretization.
 *
 * The function is used by nse specific functions to convert discretization type from 
 * \link fnft_nse_discretization_t \endlink to \link fnft__akns_discretization_t \endlink.
 * @param[in] nse_discretization The type of nse discretization. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @param[out] akns_discretization The pointer to the converted discretization
 * of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup nse
 */
FNFT_INT fnft__nse_discretization_to_akns_discretization(fnft_nse_discretization_t nse_discretization, 
        fnft__akns_discretization_t * const akns_discretization);
        
        
/**
 * @brief This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 * 
 * This routine maps continuous-time domain value lambda to discrete-time domain value
 * z = exp(2i*lambda*eps_t/degree1step), where degree1step is based on the discretization 
 * of type \link fnft_nse_discretization_t \endlink. Changes discretization to 
 * \link fnft__akns_discretization_t \endlink type and calls \link fnft__akns_lambda_to_z \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued continuous-time domain spectral parameter lambda. The values are replaced with
 * discrete-time domain values z.
 * @param[in] discretization Discretization of type \link fnft_nse_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup nse
 */
FNFT_INT fnft__nse_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft_nse_discretization_t discretization);

/**
 * @brief This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 * 
 * This routine maps discrete-time domain value z to continuous-time domain value
 * lambda = degree1step*log(z)/(2i*eps_t), where degree1step is based on the discretization 
 * of type \link fnft_nse_discretization_t \endlink. Changes discretization to 
 * \link fnft__akns_discretization_t \endlink type and calls \link fnft__akns_z_to_lambda \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued discrete-time domain spectral parameter z. The values are replaced with
 * continuous-time domain values lambda.
 * @param[in] discretization Discretization of type \link fnft_nse_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup nse
 */
FNFT_INT fnft__nse_z_to_lambda(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_discretization_degree(...) fnft__nse_discretization_degree(__VA_ARGS__)
#define nse_discretization_boundary_coeff(...) fnft__nse_discretization_boundary_coeff(__VA_ARGS__)
#define nse_discretization_to_akns_discretization(...) fnft__nse_discretization_to_akns_discretization(__VA_ARGS__)
#define nse_lambda_to_z(...) fnft__nse_lambda_to_z(__VA_ARGS__)
#define nse_z_to_lambda(...) fnft__nse_z_to_lambda(__VA_ARGS__)
#endif

#endif
