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
 * @file fnft__poly_roots_fftgridsearch.h
 * @brief Approximation of polynomial roots on the unit circle using gridding.
 * @ingroup poly
 */

#ifndef FNFT__POLY_ROOTS_FFTGRIDSEARCH_H
#define FNFT__POLY_ROOTS_FFTGRIDSEARCH_H

#include "fnft.h"

/**
 * @brief Unit circle roots of a polynomial via grid search.
 * 
 * @ingroup poly
 * This routine approximates the unit circle roots of a polynomial
 *
 *   \f[ p(z)=p_0 + p_{1}z^1 + \dots p_{deg} + p_N z^{deg}, \f]
 *
 * where the star denote the complex conjugate. (That is, it looks for the
 * solutions of \f$ p(z)=0 \f$ that satisfy \f$ |z|=1 \f$.) The polynomial
 * \f$ p(z) \f$ is evaluated on the grid \f$ z=r_k e^{j\theta_m} \f$,
 * where
 *
 *  \f[ \theta_m=\Phi_0 + m \epsilon, \ \ \ \ m=0,1,\dots,M-1,
 \ \epsilon=\frac{\Phi_1 - \Phi_0}{M-1}\f]
 *
 * and
 *
 * \f[ r_k = 1 - k\epsilon, \ \ \ \ k=-1,0,1. \f]
 *
 * The evaluations are carried out using the Chirp transform. The grid points
 * on the unit circle where the absolute value of the polynomial is lower than
 * or equal to the absolute value at the eight neighboring grid points are
 * candidates for a root. (This is the first stage of the Lindsey-Fox
 * algorithm.) The polynomial is then approximated locally around the
 * center grid poFNFT_INT and the root of this linerization is computed. The
 * root of the linearization is kept if it is not too far from the grid point.
 *
 * @see https://doi.org/10.1109/MSP.2003.1253552
 * @see poly_roots_fftgridsearch_paraherm
 * @see poly_chirpz
 *
 * @param[in] deg The degree of the polynomial.
 * @param[in] p Array containing the deg+1 coefficients of the Laurent
 *  polynomial in descending order (i.e.,
 *  \f$ p_{deg}, p_{deg-1}, \dots, p_{1}, p_{0} \f$).
 * @param[in,out] M_ptr Upon entry, *M_ptr contains the desired value for the
 *  number of points M. Upon return, *M_ptr has been overwritten with the
 *  number of detected roots.
 * @param[in] PHI Array with two entries, \f$ \Phi_0 \f$ and \f$ \Phi_1 \f$.
 *  The first value should be lower than the second one.
 * @param[out] roots Array of M points. Will be filled with the detected roots.
 *  (The number of detected roots is stored in *M_ptr.)
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_roots_fftgridsearch(const FNFT_UINT deg,
    FNFT_COMPLEX const * const p, FNFT_UINT * const M_ptr,
    FNFT_REAL const * const PHI, FNFT_COMPLEX * const roots);

/**
 * @brief Unit circle roots of a parahermitian Laurent polynomial via grid
 *  search.
 * @ingroup poly
 *
 *
 * This routine approximates the unit circle roots of a para-hermitian
 * Laurent polynomial
 *
 *   \f[ p(z)=p_{-N} + p_{-N+1} + \dots p_{N-1} + p_N,
\ \ \ \ p(z)=p^*(1/z^*), \f]
 *
 * where the star denote the complex conjugate. (That is, it looks for the
 * solutions of \f$ p(z)=0 \f$ that satisfy \f$ |z|=1 \f$.) The Laurent
 * polynomial \f$ p(z) \f$ is evaluated on the grid \f$ z=e^{j\theta_m} \f$,
 * where
 *
 *  \f[ \theta_m=\Phi_0 + m \epsilon, \ \ \ \ m=0,1,\dots,M-1,
\ \epsilon=\frac{\Phi_1 - \Phi_0}{M-1}.\f]
 *
 * The evaluations are carried out using the Chirp transform. Since the Laurent
 * polynomial is para-hermitian, it will be
 * real on the unit circle. If the sign of the Laurent polynomial changes
 * between two consequtive points on the grid, a root is detected. The
 * position of the roots is determined using linear interpolation.
 *
 * @see poly_roots_fftgridsearch
 * @see poly_chirpz
 *
 * @param[in] deg This should be the value of \f$ 2N+1 \f$.
 * @param[in] p Array containing the deg+1 coefficients of the Laurent
 *  polynomial in descending order (i.e.,
 *  \f$ p_{N}, p_{N-1}, \dots, p_{-N+1}, p_{-N} \f$).
 * @param[in,out] M_ptr Upon entry, *M_ptr contains the desired value for the
 *  number of points M. Upon return, *M_ptr has been overwritten with the
 *  number of detected roots.
 * @param[in] PHI Array with two entries, \f$ \Phi_0 \f$ and \f$ \Phi_1 \f$.
 *  The first value should be lower than the second one.
 * @param[out] roots Array of M points. Will be filled with the detected roots.
 *  (The number of detected roots is stored in *M_ptr.)
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_roots_fftgridsearch_paraherm(const FNFT_UINT deg,
    FNFT_COMPLEX const * const p, FNFT_UINT * const M_ptr,
    FNFT_REAL const * const PHI, FNFT_COMPLEX * const roots);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_roots_fftgridsearch(...) fnft__poly_roots_fftgridsearch(__VA_ARGS__)
#define poly_roots_fftgridsearch_paraherm(...) fnft__poly_roots_fftgridsearch_paraherm(__VA_ARGS__)
#endif

 
#endif
