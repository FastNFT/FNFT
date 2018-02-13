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
 * @file fnft__misc.h
 * @ingroup misc
 */

#ifndef FNFT__MISC_H
#define FNFT__MISC_H

#include "fnft.h"

/**
 * @brief Helper function for debugging. Prints an array in MATLAB style.
 * @ingroup misc
 * This function prints an array in MATLAB style.
 * @param len Length of the array to be printed.
 * @param buf Array to be printed.
 * @param varname Name of the array being printed.
 */
void fnft__misc_print_buf(FNFT_INT len, FNFT_COMPLEX *buf, char *varname);

/**
 * @brief Relative l1 error between two vectors.
 * @ingroup misc
 * This function computes the relative l1 error between two vectors.\n
 * \f$err = \frac{\sum_{i=0}^{i=len-1} |vec\_numer[i]-vec\_exact[i]|}{\sum_{i=0}^{i=len-1} |vec\_exact[i]|}\f$.
 * @param len Length of the vectors
 * @param vec_numer Complex array of numerically computed result of length len.
 * @param vec_exact Complex array of exact result of length len.
 * @return Returns the real valued relative error err.
 */
FNFT_REAL fnft__misc_rel_err(FNFT_INT len, FNFT_COMPLEX *vec_numer,
    FNFT_COMPLEX *vec_exact);

/**
 * @brief Hausdorff distance between two vectors.
 * @ingroup misc
 * This function computes the Hausdorff distance between two vectors vecA and vecB.
 * @param lenA Length of vector vecA.
 * @param vecA Complex vector of length lenA.
 * @param lenB length of vector vecB.
 * @param vecB Complex vector of length lenB.
 * @return Returns the real valued Hausdorff distance between the vectors vecA and vecB.
 */
FNFT_REAL fnft__misc_hausdorff_dist(const FNFT_UINT lenA,
    FNFT_COMPLEX const * const vecA, const FNFT_UINT lenB,
    FNFT_COMPLEX const * const vecB);

/**
 * Hyperbolic secant.
 * @ingroup misc
 */
FNFT_COMPLEX fnft__misc_sech(FNFT_COMPLEX Z);

/**
 * @brief Squared l2 norm. 
 * @ingroup misc
 * This function computes the quantity\n
 * \f$ val = \frac{b-a}{2N}.(|Z[0]|^2+|Z[N-1]|^2)+\sum_{i=1}^{i=N-2}\frac{b-a}{N}.|Z[i]|^2\f$.
 * @param N Number of elements in the array.
 * @param Z Complex valued array of length N.
 * @param a Real number corresponding to first element of Z.
 * @param b Real number corresponding to last element of Z.
 * @return Returns the quantity val. Returns NAN if N<2 or a>=b.
 */
FNFT_REAL fnft__misc_l2norm2(const FNFT_UINT N, FNFT_COMPLEX const * const Z,
    const FNFT_REAL a, const FNFT_REAL b);

/**
 * This function filters the array vals. Only values that satisfy
 *
 *      bounding_box[0] <= real(val) <= bounding_box[1]
 * and
 *
 *      bounding_box[2] <= imag(val) <= bounding_box[3]
 *
 * are kept. *N_ptr is the number of values. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * move to the beginning of vals. 
 *
 * If the array rearrange_as_well is not NULL, then the values in there are
 * rearranged together with the values in vals. It is ignored if ==NULL.
 *
 * Returns SUCCESS or an error code.
 * @ingroup misc
 */
FNFT_INT fnft__misc_filter(FNFT_UINT * const N, FNFT_COMPLEX * const vals,
    FNFT_COMPLEX * const rearrange_as_well,
    FNFT_REAL const * const bounding_box);

/**
 * This function filters the array vals. Only values OUTSIDE the bounding box
 * are kept. *N_ptr is the number of values. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * move to the beginning of vals. 
 *
 * If the array rearrange_as_well is not NULL, then the values in there are
 * rearranged together with the values in vals. It is ignored if ==NULL.
 *
 * Returns SUCCESS or an error code.
 * @ingroup misc
 */
FNFT_INT fnft__misc_filter_nonreal(FNFT_UINT *N_ptr, FNFT_COMPLEX * const vals,
    const FNFT_REAL tol_im);

/**
 * This function removes all entries from the array vals with |Im(val)|>tol.
 * *N_ptr is the number of values. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * move to the beginning of vals. 
 *
 * Returns SUCCESS or an error code.
 * @ingroup misc
 */
FNFT_INT fnft__misc_filter_inv(FNFT_UINT * const N_ptr, FNFT_COMPLEX * const vals,
    FNFT_COMPLEX * const rearrange_as_well,
    FNFT_REAL const * const bounding_box);

/**
 * Merges elements in an array with distance lower than tol.
 * @ingroup misc
 */
FNFT_INT fnft__misc_merge(FNFT_UINT *N_ptr, FNFT_COMPLEX * const vals,
    FNFT_REAL tol);

/**
 * Computes a subsampled version of q. The length of q is D>=2. The routine
 * will allocate memory for the subsampled signal qsub and updates the
 * pointer *qsub_ptr such that it points to the newly allocated qsub. The
 * user is responsible to freeing the memory later. The new number of samples
 * Dsub>=2 and the subsampling factor D/Dsub are stored in *Dsub_ptr and
 * *subsampling_factor_ptr. Returns SUCCESS or an error code.
 * @ingroup misc
 */
FNFT_INT fnft__misc_downsample(FNFT_COMPLEX const * const q, const FNFT_UINT D,
    FNFT_COMPLEX ** qsub_ptr, FNFT_UINT * const Dsub_ptr,
    FNFT_UINT * const subsampling_factor_ptr);

/**
 * Sinc function for complex arguments.
 * @ingroup misc
 */
FNFT_COMPLEX fnft__misc_CSINC(FNFT_COMPLEX x);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define misc_print_buf(...) fnft__misc_print_buf(__VA_ARGS__)
#define misc_rel_err(...) fnft__misc_rel_err(__VA_ARGS__)
#define misc_hausdorff_dist(...) fnft__misc_hausdorff_dist(__VA_ARGS__)
#define misc_sech(...) fnft__misc_sech(__VA_ARGS__)
#define misc_l2norm2(...) fnft__misc_l2norm2(__VA_ARGS__)
#define misc_filter(...) fnft__misc_filter(__VA_ARGS__)
#define misc_filter_inv(...) fnft__misc_filter_inv(__VA_ARGS__)
#define misc_filter_nonreal(...) fnft__misc_filter_nonreal(__VA_ARGS__)
#define misc_merge(...) fnft__misc_merge(__VA_ARGS__)
#define misc_downsample(...) fnft__misc_downsample(__VA_ARGS__)
#define misc_CSINC(...) fnft__misc_CSINC(__VA_ARGS__)
#endif

#endif
