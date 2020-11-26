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
 * Shrinivas Chimmalgi (TU Delft) 2017.
 * Peter J. Prins (TU Delft) 2018.
 */

/**
 * @brief Properties of the discretizations for the Korteweg-de Vries equation.
 *
 * @file fnft__akns_discretization.h
 * @ingroup akns
 */
#ifndef FNFT__AKNS_DISCRETIZATION_H
#define FNFT__AKNS_DISCRETIZATION_H

#include "fnft__akns_discretization_t.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"


/**
 * @brief This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 *
 * It defines the step size of the frequency grid
 * \f$z = \text{e}^{2j\xi\epsilon_t d}\f$ based on the discretization type.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft_kdv_discretization_t \endlink.
 * @returns polynomial degree, or 0 for discretizations not supported by
 * \link fnft__akns_fscatter \endlink.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_degree(fnft__akns_discretization_t
        discretization);

/**
 * @brief This routine returns the boundary coefficient based on the
 * discretization.
 *
 * The boundary coefficient is the fraction of the step size that a discretized
 * potential extends beyond the last sample. This routine returns this value
 * based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the boundary coefficient, or NAN for discretizations not supported
 * by \link fnft__akns_fscatter \endlink.
 *
 * @ingroup akns
 */
FNFT_REAL fnft__akns_discretization_boundary_coeff(fnft__akns_discretization_t discretization);

/**
 * @brief This routine returns the scaling for effective number of samples based on the
 * discretization.
 *
 * Higher order methods use more than one sample per integration step. This routine returns
 * the value upsampling_factor based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * D_effective = upsampling_factor * D.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the upsampling_factor value, or 0 for unknown discretizations.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_upsampling_factor(fnft__akns_discretization_t discretization);

/**
 * @brief This routine returns the order of the method based on the
 * discretization.
 *
 * Different numerical methods have different orders of accuray. This routine returns
 * the order of the order based on the discretization of type \link fnft__akns_discretization_t \endlink.
 * When the step-size of the signal samples is reduced by a factor \f$s\f$, the error in the
 * computed values is expected to decrease by a factor \f$s^{order}\f$.
 * @param[in] discretization The type of discretization to be used. Should be
 * of type \link fnft__akns_discretization_t \endlink.
 * @returns the method_order value, or 0 for unknown discretization.
 *
 * @ingroup akns
 */
FNFT_UINT fnft__akns_discretization_method_order(fnft__akns_discretization_t discretization);

/**
 * @brief This routine maps \f$\lambda\f$ from continuous-time domain to
 * \f$z\f$ in the discrete-time domain based on the discretization. 
 * 
 * This routine maps continuous-time domain value \f$\lambda\f$ to discrete-time domain value
 * \f$z = e^{2j \lambda \epsilon_t degree1step}\f$, where degree1step is based on the discretization 
 * of type \link fnft__akns_discretization_t \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued continuous-time domain spectral parameter \f$\lambda\f$. The values are replaced with
 * discrete-time domain values \f$z\f$.
 * @param[in] discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft__akns_discretization_t discretization);

/**
 * @brief This routine maps \f$z\f$ from the discrete-time domain to
 * \f$\lambda\f$ in the continuous-time domain based on the discretization. 
 * 
 * This routine maps discrete-time domain value \f$z\f$ to continuous-time domain value
 * \f$\lambda = degree1step\log(z)/(2j\epsilon_t)\f$, where degree1step is based on the discretization 
 * of type \link fnft__akns_discretization_t \endlink.
 * @param[in] n Number of values to be mapped.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[in,out] vals Pointer to location of first element of array containing
 * complex-valued discrete-time domain spectral parameter \f$z\f$. The values are replaced with
 * continuous-time domain values \f$\lambda\f$.
 * @param[in] discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_z_to_lambda(const FNFT_UINT n, const FNFT_REAL eps_t, 
        FNFT_COMPLEX * const vals, fnft__akns_discretization_t discretization);

/**
 * @brief This routine computes various weights required by some methods
 * based on the discretization. 
 * 
 * This routing computes the special weights required for the 
 * higher-order methods CF\f$^{[4]}_2\f$, CF\f$^{[4]}_3\f$, CF\f$^{[5]}_3\f$ 
 * and CF\f$^{[6]}_4\f$. The weights are used in \link fnft__nse_discretization_preprocess_signal \endlink,
 * \link fnft__akns_scatter_matrix \endlink and \link fnft__nse_scatter_bound_states\endlink.
 * The weights for CF\f$^{[4]}_3\f$ are taken from Alvermann and Fehske (<a href="https://doi.org/10.1016/j.jcp.2011.04.006">Journal of Computational Phys. 230, 2011</a>)
 * and the weights for the others are from Blanes, Casas and Thalhammer(<a href="https://doi.org/10.1016/j.cpc.2017.07.016">Computer Phys. Comm. 220, 2017</a>).
 * The weights are mentioned as matrices in the references. This routine returns 
 * them in row-major order.
 * @param[in,out] weights_ptr Pointer to the starting location of weights.
 * @param[in] akns_discretization Discretization of type \link fnft__akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT fnft__akns_discretization_method_weights(FNFT_COMPLEX **weights_ptr,
        fnft__akns_discretization_t akns_discretization);

/**
 * @brief  This routine preprocesses the signal by resampling and subsampling based on the discretization.
 * The preprocessing is necessary for higher-order methods.
 *
 * This routine preprocess q to generate q_preprocessed and r_preprocessed
 * based on the discretization. The preprocessing may involve resampling
 * and sub-sampling.
 * The routine is based on the following papers:
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot; Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] r_from_q A function pointer array of length 3. `r_from_q[0](q)` has to return the value of an r-sample given the corresponding q-sample, before preprocessing. `r_from_q[1](q)` and `r_from_q[2](q)` have to return respectively the sample of the first/second derivative of r given the corresponding sample of the first/second derivative of q.
 * @param[in] eps_t Real-valued discretization step-size.
 * @param[out] q_preprocessed_ptr Pointer to the starting location of preprocessed signal q_preprocessed.
 * @param[out] r_preprocessed_ptr Pointer to the starting location of preprocessed signal r_preprocessed.
 * @param[in,out] Dsub_ptr Pointer to number of processed samples. Upon entry, *Dsub_ptr
 *             should contain a desired number of samples. Upon exit, *Dsub_ptr
 *             has been overwritten with the actual number of samples that the
 *             routine has chosen. It is usually close to the desired one.
 * @param[out] first_last_index Vector of length two. Upon exit, it contains
 *             the original index of the first and the last sample used to build
 *             q_preprocessed.
 * @param[in] discretization Discretization of type \link fnft_akns_discretization_t \endlink.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
*/
FNFT_INT fnft__akns_discretization_preprocess_signal(FNFT_UINT const D,
                                                     FNFT_COMPLEX const * const q,
                                                     FNFT_COMPLEX (*r_from_q[3])(FNFT_COMPLEX const),
                                                     FNFT_REAL const eps_t,
                                                     FNFT_UINT * const Dsub_ptr,
                                                     FNFT_COMPLEX **q_preprocessed_ptr,
                                                     FNFT_COMPLEX **r_preprocessed_ptr,
                                                     FNFT_UINT * const first_last_index,
                                                     fnft__akns_discretization_t discretization);
        
#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_discretization_degree(...) fnft__akns_discretization_degree(__VA_ARGS__)
#define akns_discretization_boundary_coeff(...) fnft__akns_discretization_boundary_coeff(__VA_ARGS__)
#define akns_discretization_upsampling_factor(...) fnft__akns_discretization_upsampling_factor(__VA_ARGS__)
#define akns_discretization_method_order(...) fnft__akns_discretization_method_order(__VA_ARGS__)
#define akns_discretization_lambda_to_z(...) fnft__akns_discretization_lambda_to_z(__VA_ARGS__)
#define akns_discretization_z_to_lambda(...) fnft__akns_discretization_z_to_lambda(__VA_ARGS__)
#define akns_discretization_method_weights(...) fnft__akns_discretization_method_weights(__VA_ARGS__)
#define akns_discretization_preprocess_signal(...) fnft__akns_discretization_preprocess_signal(__VA_ARGS__)

#endif

#endif
