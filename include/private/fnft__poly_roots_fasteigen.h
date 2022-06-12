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
 * @file fnft__poly_roots_fasteigen.h
 * @brief Fast root finding of polynomials.
 * @ingroup poly
 */

#ifndef FNFT__POLY_ROOTS_FASTEIGEN_H
#define FNFT__POLY_ROOTS_FASTEIGEN_H

#include "fnft.h"

/**
 * @brief Fast computation of polynomial roots.
 * 
 * @ingroup poly
 * @ingroup roots
 * This routine compute the roots of a polynomial
 *
 *   \f[ p(z)=p_0+p_1 z^1+p_2 z^2+...+p_{deg} z^{deg} \f]
 *
 * using only \f$ O\{ deg^2 \}\f$ floating point operations.
 *
 * @see https://arxiv.org/abs/1611.02435v2
 *
 * @param[in] deg Degree of the polynomial
 * @param[in] p Array containing the deg+1 coefficients of the polynomial in
 *  descending order (i.e., \f$ p_{deg}, p_{deg-1}, \dots, p_1, p_0 \f$).
 * @param[out] roots Array of deg+1 points. Will be filled with the roots of
 *  \f$ p(z) \f$.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_roots_fasteigen(const FNFT_UINT deg,
    FNFT_COMPLEX const * const p, FNFT_COMPLEX * const roots);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_roots_fasteigen(...) fnft__poly_roots_fasteigen(__VA_ARGS__)
#endif

#endif
