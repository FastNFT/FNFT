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
*/

/**
 * @file fnft__poly_eval.h
 * @brief Evaluation of polynomials and their derivatives using Horner's method.
 * @ingroup poly
 */

#ifndef FNFT__POLY_EVAL_H
#define FNFT__POLY_EVAL_H

#include "fnft.h"

/**
 * @brief Evaluation of polynomials.
 * @ingroup poly
 *
 * Evaluates a polynomial
 *
 *   \f[ p(z)=p_0+p_1 z^1+p_2 z^2+...+p_{deg} z^{deg} \f]
 *
 * at \a nz points
 *
 *  \f[ z=z_1,...,z_{nz} \f]
 *
 * in the complex plane using Horner's method.
 *
 * @see http://cnx.org/contents/hnzXCQRy@7/Horners-Method-for-Evaluating-
 *
 * @param[in] deg Degree of the polynomial
 * @param[in] p Array containing the deg+1 coefficients of the polynomial in
 *  descending order (i.e., \f$ p_{deg}, p_{deg-1}, \dots, p_1, p_0 \f$).
 * @param[in] nz Number of points z at which the polynomial should be evaluated
 * @param[in,out] z Array of \a nz points \f$z_1,\dots,z_{nz}\f$ at which the
 *  polynomial should be evaluated. These values will be overwritten with
 *  the corresponding values \f$p(z_1),\dots,p(z_{nz})\f$ of the polynomial.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_eval(const FNFT_UINT deg, FNFT_COMPLEX const * const p, \
    const FNFT_UINT nz, FNFT_COMPLEX * const z);

/**
 * @brief Evaluation of polynomials and their derivatives.
 * @ingroup poly
 *
 * Evaluates a polynomial
 *
 *   \f[ p(z)=p_0+p_1 z^1+p_2 z^2+...+p_{deg} z^{deg} \f]
 *
 * and its derivative
 *
 *   \f[ p'(z)=\frac{dp}{dz}(z)=p_1+2p_2 z^1+...+deg~p_{deg} z^{deg-1} \f]
 *
 * at \a nz points
 *
 *  \f[ z=z_1,...,z_{nz} \f]
 *
 * in the complex plane using Horner's method.
 *
 * @see http://cnx.org/contents/hnzXCQRy@7/Horners-Method-for-Evaluating-
 *
 * @param[in] deg Degree of the polynomial
 * @param[in] p Array containing the deg+1 coefficients of the polynomial in
 *  descending order (i.e., \f$ p_{deg}, p_{deg-1}, \dots, p_1, p_0 \f$).
 * @param[in] nz Number of points z at which the polynomial should be evaluated
 * @param[in,out] z Array of \a nz points \f$z_1,\dots,z_{nz}\f$ at which the
 *  polynomial and its derivative should be evaluated. These values will be
 *  overwritten with the corresponding values \f$p(z_1),\dots,p(z_{nz})\f$ of
 *  the polynomial.
 * @param[out] deriv Array of length \a nz. It will be filled with the values
 *  \f$p'(z_1),\dots,p'(z_{nz})\f$ of the derivative.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_evalderiv(const FNFT_UINT deg, FNFT_COMPLEX const * const p, \
    const FNFT_UINT nz, FNFT_COMPLEX * const z, \
    FNFT_COMPLEX * const deriv);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_eval(...) fnft__poly_eval(__VA_ARGS__)
#define poly_evalderiv(...) fnft__poly_evalderiv(__VA_ARGS__)
#endif

#endif
