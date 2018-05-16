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
#include "fnft__fft_wrapper_plan_t.h"

FNFT_INT fnft__poly_fmult_two_polys_len(const FNFT_UINT deg);
FNFT_UINT fnft__poly_fmult_two_polys_lenmen(const FNFT_UINT deg);
FNFT_INT fnft__poly_fmult_two_polys(
    const FNFT_UINT deg,
    FNFT_COMPLEX const * const p1, 
    FNFT_COMPLEX const * const p2,
    FNFT_COMPLEX * const result,
    fnft__fft_wrapper_plan_t plan_fwd,
    fnft__fft_wrapper_plan_t plan_inv,
    FNFT_COMPLEX * const buf0,
    FNFT_COMPLEX * const buf1,
    FNFT_COMPLEX * const buf2,
    const FNFT_INT add_flag);
FNFT_INT fnft__poly_fmult_two_polys2x2(const FNFT_UINT deg,
    FNFT_COMPLEX const * const p1_11,
    const FNFT_UINT p1_stride,
    FNFT_COMPLEX const * const p2_11,
    const FNFT_UINT p2_stride,
    FNFT_COMPLEX * const result_11,
    const FNFT_UINT result_stride,
    fnft__fft_wrapper_plan_t plan_fwd,
    fnft__fft_wrapper_plan_t plan_inv,
    FNFT_COMPLEX * const buf0,
    FNFT_COMPLEX * const buf1,
    FNFT_COMPLEX * const buf2);

/**
 * @brief Number of elements that the input p to
 * \link fnft__poly_fmult \endlink should have.
 *
 * @ingroup poly
 * Specifies how much memory (in number of elements) the user needs to allocate
 * for the input p of the routine \link fnft__poly_fmult \endlink.
 * @param [in] deg Degree of the polynomials
 * @param [in] n Number of polynomials
 * @return A number m. The input p to \link fnft__poly_fmult \endlink should be
 * a array with m entries.
 */
FNFT_UINT fnft__poly_fmult_numel(const FNFT_UINT deg, const FNFT_UINT n);

/**
 * @brief Fast multiplication of multiple polynomials of same degree.
 * 
 * @ingroup poly
 * Fast multiplication of n polynomials of degree d. Their coefficients are
 * stored in the array p and will be overwritten. If W_ptr != NULL, the
 * result has been normalized by a factor 2^W. Upon exit, W has been stored
 * in *W_ptr.
 * @param[in, out] d Upon entry, degree of the input polynomials. Upon exit,
 *  degree of their product.
 * @param[in] n Number of polynomials.
 * @param[in,out] p Complex valued array with m entries, where m is determined
 *  using \link fnft__poly_fmult_memneeded \endlink. Upon entry, the first
 *  (*d+1)*n elements of this array contain the coefficients of the
 *  polynomials. Upon exit, the first *d+1 elements contain the result.
 * @param[in] W_ptr Pointer to normalization flag. 
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_fmult(FNFT_UINT * const d, FNFT_UINT n, FNFT_COMPLEX * const p,
    FNFT_INT * const W_ptr);

/**
 * @brief Number of elements that the inputs p and result to
 * \link fnft__poly_fmult2x2 \endlink should have.
 *
 * @ingroup poly
 * Specifies how much memory (in number of elements) the user needs to allocate
 * for the inputs p and result of the routine
 * \link fnft__poly_fmult2x2 \endlink.
 * @param [in] deg Degree of the polynomials
 * @param [in] n Number of polynomials
 * @return A number m. The inputs p and result to
 * \link fnft__poly_fmult2x2 \endlink should be a array with m entries.
 */
FNFT_UINT fnft__poly_fmult2x2_numel(const FNFT_UINT deg, const FNFT_UINT n);

/**
 * @brief Fast multiplication of multiple 2x2 matrix-valued polynomials of same degree.
 * 
 * @ingroup poly
 * Fast multiplication of n 2x2 matrix-valued polynomials of degree d. Their
 * coefficients are stored in the array p and will be overwritten. If
 * W_ptr != NULL, the result has been normalized by a factor 2^W. Upon exit,
 * W has been stored in *W_ptr.
 * @param[in] d Pointer to a \link FNFT_UINT \endlink containing the degree of
 * the polynomials.
 * @param[in] n Number of 2x2 matrix-valued polynomials.
 * @param[in,out] p Complex valued array which holds the coefficients of 
 * the polynomials being multiplied. Should be of length m*(*d+1), where
 * m is obtained using \link fnft__poly_fmult_memneeded \lendlink.
 * WARNING: p is overwritten.
 * @param[out] result Complex valued array that holds the result of the
 * multiplication. Should be of the same size as p.
 * @param[in] W_ptr Pointer to normalization flag. 
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_fmult2x2(FNFT_UINT *d, FNFT_UINT n, FNFT_COMPLEX * const p, 
    FNFT_COMPLEX * const result, FNFT_INT * const W_ptr);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_fmult_two_polys_len(...) fnft__poly_fmult_two_polys_len(__VA_ARGS__)
#define poly_fmult_two_polys_lenmen(...) fnft__poly_fmult_two_polys_lenmen(__VA_ARGS__)
#define poly_fmult_two_polys(...) fnft__poly_fmult_two_polys(__VA_ARGS__)
#define poly_fmult_two_polys2x2(...) fnft__poly_fmult_two_polys2x2(__VA_ARGS__)
#define poly_fmult_numel(...) fnft__poly_fmult_numel(__VA_ARGS__)
#define poly_fmult2x2_numel(...) fnft__poly_fmult2x2_numel(__VA_ARGS__)
#define poly_fmult(...) fnft__poly_fmult(__VA_ARGS__)
#define poly_fmult2x2(...) fnft__poly_fmult2x2(__VA_ARGS__)
#endif

#endif
