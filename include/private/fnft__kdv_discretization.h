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
 * Shrinivas Chimmalgi (TU Delft) 2017-2019.
 * Peter J. Prins (TU Delft) 2018.
 */

/**
 * @brief Properties of the discretizations for the Korteweg-de Vries equation.
 *
 * @file fnft__kdv_discretization.h
 * @ingroup kdv
 */
#ifndef FNFT__KDV_DISCRETIZATION_H
#define FNFT__KDV_DISCRETIZATION_H

#include "fnft_kdv_discretization_t.h"
#include "fnft__akns_discretization.h"


/**
 * @brief This routine returns the max degree \f$d\f$ of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 *
 * It defines the step size of the frequency grid
 * \f$z = \text{e}^{2j\xi\epsilon_t d}\f$ based on the discretization type.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_kdv_discretization_t \endlink.
 * @returns polynomial degree, or 0 for discretizations not supported by
 * \link fnft__kdv_fscatter \endlink.
 *
 * @ingroup kdv
 */
FNFT_UINT fnft__kdv_discretization_degree(fnft_kdv_discretization_t
        discretization);

/**
 * @brief This routine returns the boundary coefficient based on the
 * discretization.
 *
 * The boundary coefficient is the fraction of the step size that a discretized
 * potential extends beyond the last sample. This routine returns this value
 * based on the discretization of type \link fnft_kdv_discretization_t \endlink.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_kdv_discretization_t \endlink.
 * @returns the boundary coefficient, or NAN for discretizations not supported
 * by \link fnft__kdv_fscatter \endlink.
 *
 * @ingroup kdv
 */
FNFT_REAL fnft__kdv_discretization_boundary_coeff(fnft_kdv_discretization_t discretization);


/**
 * @brief This routine returns akns discretization related to the given
 * kdv discretization.
 *
 * The function is used by kdv specific functions to convert discretization type from 
 * \link fnft_kdv_discretization_t \endlink to \link fnft__akns_discretization_t \endlink.
 * @param[in] kdv_discretization The type of kdv discretization. Should be
 * of type \link fnft_kdv_discretization_t \endlink.
 * @param[out] akns_discretization The pointer to the converted discretization
 * of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup kdv
 */
FNFT_INT fnft__kdv_discretization_to_akns_discretization(fnft_kdv_discretization_t kdv_discretization, 
        fnft__akns_discretization_t * const akns_discretization);

/**
 * @brief This routine returns the scaling for effective number of samples based on the
 * discretization.
 *
 * Higher order methods use more than one sample per integration step. This routine returns
 * the value upsampling_factor based on the discretization of type \link fnft_kdv_discretization_t \endlink.
 * D_effective = upsampling_factor * D.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_nse_discretization_t \endlink.
 * @returns the upsampling_factor value, or 0 for discretizations not supported
 * by \link fnft__kdv_fscatter \endlink.
 *
 * @ingroup kdv
 */
FNFT_UINT fnft__kdv_discretization_upsampling_factor(fnft_kdv_discretization_t discretization);
        
        
/**
 * @brief This routine maps \f$\lambda\f$ from continuous-time domain to
 * \f$z\f$ in the discrete-time domain based on the discretization. 
 * 
 * This routine maps continuous-time domain value \f$\lambda\f$ to discrete-time domain value
 * \f$z = e^{2j\lambda\epsilon_t degree1step}\f$, where degree1step is based on the discretization 
 * of type \link fnft_kdv_discretization_t \endlink. Changes discretization to 
 * \link fnft__akns_discretization_t \endlink type and calls \link fnft__akns_discretization_lambda_to_z \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued continuous-time domain spectral parameter \f$\lambda\f$. The values are replaced with
 * discrete-time domain values \f$z\f$.
 * @param[in] discretization Discretization of type \link fnft_kdv_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup kdv
 */
FNFT_INT fnft__kdv_discretization_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft_kdv_discretization_t discretization);

/**
 * @brief This routine maps \f$z\f$ from the discrete-time domain to
 * \f$\lambda\f$ in the continuous-time domain based on the discretization. 
 * 
 * This routine maps discrete-time domain value \f$z\f$ to continuous-time domain value
 * \f$\lambda = degree1step\log(z)/(2j\epsilon_t)\f$, where degree1step is based on the discretization 
 * of type \link fnft_kdv_discretization_t \endlink. Changes discretization to 
 * \link fnft__akns_discretization_t \endlink type and calls \link fnft__akns_discretization_z_to_lambda \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued discrete-time domain spectral parameter \f$z\f$. The values are replaced with
 * continuous-time domain values \f$\lambda\f$.
 * @param[in] discretization Discretization of type \link fnft_kdv_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup kdv
 */
FNFT_INT fnft__kdv_discretization_z_to_lambda(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft_kdv_discretization_t discretization);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_discretization_degree(...) fnft__kdv_discretization_degree(__VA_ARGS__)
#define kdv_discretization_boundary_coeff(...) fnft__kdv_discretization_boundary_coeff(__VA_ARGS__)
#define kdv_discretization_to_akns_discretization(...) fnft__kdv_discretization_to_akns_discretization(__VA_ARGS__)
#define kdv_discretization_upsampling_factor(...) fnft__kdv_discretization_upsampling_factor(__VA_ARGS__)
#define kdv_discretization_lambda_to_z(...) fnft__kdv_discretization_lambda_to_z(__VA_ARGS__)
#define kdv_discretization_z_to_lambda(...) fnft__kdv_discretization_z_to_lambda(__VA_ARGS__)
#endif

#endif
