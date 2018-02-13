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
 * @file fnft__poly_fmult.h
 * @brief Fast multiplication of n polynomials.
 * @ingroup poly
 */

#ifndef FNFT__POLY_FMULT_H
#define FNFT__POLY_FMULT_H

#include "fnft.h"

/**
 * Fast multiplication of n polynomials of degree d. Their coefficients are
 * stored in the array p and will be overwritten. If W_ptr != NULL, the
 * result has been normalized by a factor 2^W. Upon exit, W has been stored
 * in *W_ptr.
 * @returns \link FNFT_SUCCESS \endlink or an error code.
 */
FNFT_INT fnft__poly_fmult(FNFT_UINT * const d, FNFT_UINT n, FNFT_COMPLEX * const p,
    FNFT_INT * const W_ptr);

/**
 * Fast multiplication of n 2x2 matrix-valued polynomials of degree d. Their
 * coefficients are stored in the array p and will be overwritten. If
 * W_ptr != NULL, the result has been normalized by a factor 2^W. Upon exit,
 * W has been stored in *W_ptr.
 * @returns \link FNFT_SUCCESS \endlink or an error code.
 */
FNFT_INT fnft__poly_fmult2x2(FNFT_UINT *d, FNFT_UINT n, FNFT_COMPLEX * const p, 
    FNFT_COMPLEX * const result, FNFT_INT * const W_ptr);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_fmult(...) fnft__poly_fmult(__VA_ARGS__)
#define poly_fmult2x2(...) fnft__poly_fmult2x2(__VA_ARGS__)
#endif

#endif
