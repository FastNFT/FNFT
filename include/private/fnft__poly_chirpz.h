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
*/

/**
 * @file fnft__poly_chirpz.h
 * @brief Fast evaluation of polynomials using the Chirp Z-transform.
 * @ingroup poly
 */

#ifndef FNFT__POLY_CHIRPZ_H
#define FNFT__POLY_CHIRPZ_H

#include "fnft.h"

/**
 * @brief Fast evaluation of a polynomial on a spiral in the complex plane.\n
 * @ingroup poly
 *
 * This routine implements the Chirp Z-transform. Given a polynomial
 *
 *   \f[ p(z)=p_0+p_1 z^1+p_2 z^2+...+p_{deg} z^{deg} \f]
 *
 * and complex numbers \a A and \a W, this routine evaluates the polynomial
 * \f$ p(z) \f$ at the \a M points \f$ z=1/w_k \f$, where
 *
 *  \f[ w_k=AW^{-m}, \quad m=0,1,\dots,M-1, \f]
 *
 * using only \f$ O\{(N+M)\log(N+M)\} \f$ floating point operations.
 *
 * @see https://doi.org/10.1109/TAU.1969.1162034
 *
 * @param[in] deg Degree of the polynomial
 * @param[in] p Array containing the deg+1 coefficients of the polynomial in
 *  descending order (i.e., \f$ p_{deg}, p_{deg-1}, \dots, p_1, p_0 \f$).
 * @param[in] A First constant defining the spiral at which the polynomial will
 *  be evaluated.
 * @param[in] W Second constant defining the spiral at which the polynomial
 *  will be evaluated
 * @param[in] M Number of points at which the polynomial will be evaluated.
 * @param[out] result Array of M points. Will be filled with
 *  \f$ p(1/w_1),...,p(1/w_M) \f$.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_chirpz(const FNFT_UINT deg, FNFT_COMPLEX const * const p, \
    const FNFT_COMPLEX A, const FNFT_COMPLEX W, const FNFT_UINT M, \
    FNFT_COMPLEX * const result);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_chirpz(...) fnft__poly_chirpz(__VA_ARGS__)
#endif

#endif
