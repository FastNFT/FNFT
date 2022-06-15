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
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
* Peter J. Prins (TU Delft) 2021.
*/

/**
 * @file fnft__misc.h
 * @brief Miscellaneous functions used in the FNFT library.
 * @ingroup misc
 */

#ifndef FNFT__MISC_H
#define FNFT__MISC_H

#include "fnft.h"
#include <string.h>

/**
 * @brief Helper function for debugging. Prints an array in MATLAB style.
 *
 * @ingroup misc
 * This function prints an array in MATLAB style.
 * @param[in] len Length of the array to be printed.
 * @param[in] buf Array to be printed.
 * @param[in] varname Name of the array being printed.
 */
void fnft__misc_print_buf(const FNFT_UINT len, FNFT_COMPLEX const * const buf,
                          char const * const varname);

/**
 * @brief Relative l1 error between two vectors.
 *
 * @ingroup misc
 * This function computes the relative l1 error between two vectors.\n
 * \f$err = \frac{\sum_{i=0}^{i=len-1} |vec\_numer[i]-vec\_exact[i]|}{\sum_{i=0}^{i=len-1} |vec\_exact[i]|}\f$.
 * @param[in] len Length of the vectors
 * @param[in] vec_numer Complex array of numerically computed result of length len.
 * @param[in] vec_exact Complex array of exact result of length len.
 * @return Returns the real valued relative error err.
 */
FNFT_REAL fnft__misc_rel_err(const FNFT_UINT len,
    FNFT_COMPLEX const * const vec_numer, FNFT_COMPLEX const * const vec_exact);

/**
 * @brief Hausdorff distance between two vectors.
 *
 * @ingroup misc
 * This function computes the Hausdorff distance between two vectors vecA and vecB.
 * @param[in] lenA Length of vector vecA.
 * @param[in] vecA Complex vector of length lenA.
 * @param[in] lenB length of vector vecB.
 * @param[in] vecB Complex vector of length lenB.
 * @return Returns the real valued Hausdorff distance between the vectors vecA and vecB.
 */
FNFT_REAL fnft__misc_hausdorff_dist(const FNFT_UINT lenA,
    FNFT_COMPLEX const * const vecA, const FNFT_UINT lenB,
    FNFT_COMPLEX const * const vecB);

/**
 * @brief Hyperbolic secant.
 *
 * @ingroup misc
 * This function returns the hyperbolic secant of a \link FNFT_COMPLEX \endlink.
 * @param[in] Z \link FNFT_COMPLEX \endlink argument.
 * @return hyperbolic secant of Z.
 */
FNFT_COMPLEX fnft__misc_sech(FNFT_COMPLEX Z);

/**
 * @brief Squared l2 norm.
 *
 * @ingroup misc
 * This function computes the quantity\n
 * \f$ val = \frac{b-a}{N} \sum_{i=0}^{i=N-1}|Z[i]|^2\f$.
 * @param[in] N Number of elements in the array.
 * @param[in] Z Complex valued array of length N.
 * @param[in] a Real number corresponding to the left boundary of the first element of Z.
 * @param[in] b Real number corresponding to right boundary of the last element of Z.
 * @return Returns the quantity val. Returns NAN if N==0 or a>=b.
 */
FNFT_REAL fnft__misc_l2norm2(const FNFT_UINT N, FNFT_COMPLEX const * const Z,
    const FNFT_REAL a, const FNFT_REAL b);

/**
 * @brief Filters array by retaining elements inside a bounding box
 *
 * @ingroup misc
 * This function filters the array vals. Only values that satisfy
 *
 *      bounding_box[0] <= real(val) <= bounding_box[1]
 *
 * and
 *
 *      bounding_box[2] <= imag(val) <= bounding_box[3]
 *
 * are kept.
 * @param[in,out] N It is the pointer to the number of values to be filtered. On exit *N is overwritten with
 * the number of values that have survived fitering. Their values will be
 * moved to the beginning of vals.
 * @param[in,out] vals Complex valued array with elements to be filtered.
 * @param[in] rearrange_as_well Complex valued array. If the array rearrange_as_well is not NULL,
 * then the values in there are rearranged together with the values in vals.
 * @param[in] bounding_box A real array of 4 elements. The elements determine the corners of the
 * bounding box being used for filtering.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_filter(FNFT_UINT * const N, FNFT_COMPLEX * const vals,
    FNFT_COMPLEX * const rearrange_as_well,
    FNFT_REAL const * const bounding_box);


/**
 * @brief Filter array based on specified tolerance.
 *
 * @ingroup misc
 * This function removes all entries from the array vals with |Im(val)|>tol_im.
 * @param[in,out] N_ptr Pointer to number of values to be filtered. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * moved to the beginning of vals.
 * @param[in,out] vals Complex valued array with elements to be filtered.
 * @param[in] tol_im Real valued tolerance.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_filter_nonreal(FNFT_UINT *N_ptr, FNFT_COMPLEX * const vals,
    const FNFT_REAL tol_im);

/**
 * @brief Filters array by retaining elements outside a bounding box.
 *
 * @ingroup misc
 * This function filters the array vals. Only values OUTSIDE the bounding box
 * are kept.
 * @param[in,out] N_ptr It is the pointer to the number of values to be filtered. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * moved to the beginning of vals.
 * @param[in,out] vals Complex valued array with elements to be filtered.
 * @param[in] rearrange_as_well Complex valued array. If the array rearrange_as_well is not NULL,
 * then the values in there are rearranged together with the values in vals.
 * @param[in] bounding_box A real array of 4 elements. The elements determine the corners of the
 * bounding box being used for filtering.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_filter_inv(FNFT_UINT * const N_ptr, FNFT_COMPLEX * const vals,
    FNFT_COMPLEX * const rearrange_as_well,
    FNFT_REAL const * const bounding_box);

/**
 * @brief Merges elements in an array with distance lower than tol.
 *
 * @ingroup misc
 * This function filters an array by merging elements if distance between the elements is less than tol.
 * @param[in,out] N_ptr It is the pointer to the number of elements to be filtered. On exit *N_ptr is overwritten with
 * the number of values that have survived fitering. Their values will be
 * moved to the beginning of vals.
 * @param[in,out] vals Complex valued array with elements to be filtered.
 * @param[in] tol Real valued tolerance.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_merge(FNFT_UINT *N_ptr, FNFT_COMPLEX * const vals,
    FNFT_REAL tol);

/**
 * @brief Downsamples an array.
 *
 * @ingroup misc
 * Computes a subsampled version of q. The length of q is D>=2. The routine
 * will allocate memory for the subsampled signal qsub and updates the
 * pointer *qsub_ptr such that it points to the newly allocated qsub. The
 * user is responsible to freeing the memory later.
 * @param[in] D Number of samples in array q.
 * @param[in] q Complex valued array to be subsampled.
 * @param[out] qsub_ptr Pointer to the starting location of subsampled signal.
 * @param[out] Dsub_ptr Pointer to new number of samples. Upon entry, *Dsub_ptr
 *             should contain a desired number of samples. Upon exit, *Dsub_ptr
 *             has been overwritten with the actual number of samples that the
 *             routine has chosen. It is usually close to the desired one.
 * @param[out] first_last_index Vector of length two. Upon exit, it contains
 *             the original index of the first and the last sample in qsub.
 *             That is, qsub[0]=q[first_last_index[0]] and
 *             qsub[Dsub-1]=q[first_last_index[1]].
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_downsample(const FNFT_UINT D, FNFT_COMPLEX const * const q,
    FNFT_UINT * const Dsub_ptr, FNFT_COMPLEX ** qsub_ptr,
    FNFT_UINT * const first_last_index);

/**
 * @brief Sinc function for complex arguments.
 *
 * @ingroup misc
 * Function computes the sinc function sin(x)/x for \link FNFT_COMPLEX \endlink argument.
 * If x is close to 0, the calculation is approximated with
 * sinc(x) = cos(x/sqrt(3)) + O(x^4)
 * @param[in] x \link FNFT_COMPLEX \endlink argument.
 * @return sinc(x).
 */
static inline FNFT_COMPLEX fnft__misc_CSINC(FNFT_COMPLEX x)
{
    const FNFT_REAL sinc_th=1.0E-8;

    if (FNFT_CABS(x)>=sinc_th)
        return FNFT_CSIN(x)/x;
    else
        return FNFT_CCOS(x/FNFT_CSQRT(3));
}

/**
 * @brief Derivative of sinc function for complex arguments.
 *
 * @ingroup misc
 * Function computes the derivative of the sinc function sin(x)/x for
 * \link FNFT_COMPLEX \endlink argument.
 * This derivative can be expressed analytically as
 * d/dx sinc(x) = 0 if x=0, and (cos(x) - sinc(x)) / x otherwise.
 * However, if x is close to 0, a numerical implementation like that results in
 * catastrophic cancellation. Therefore we use its Taylor approximation instead:
 * \f$ \frac{d}{dx} \text{sinc}(x) = \sum_{n=0}^\infty \frac{x^{2n+1} (-1)^{n+1}}{(2n+3)(2*n+1)!} )\f$
 * of which we use the first 9 terms, i.e. a 18-th order Taylor approximation.
 * This polynomial is computed with Horner's method.
 * @param[in] x \link FNFT_COMPLEX \endlink argument.
 * @return d/dx sinc(x).
 */
static inline FNFT_COMPLEX fnft__misc_CSINC_derivative(FNFT_COMPLEX x)
{
    const FNFT_REAL sinc_derivative_th=1.0;
    FNFT_COMPLEX returnvalue;

    if (FNFT_CABS(x)>=sinc_derivative_th)
        returnvalue = (FNFT_CCOS(x)-fnft__misc_CSINC(x))/x;
    else {
        returnvalue = 0.0;
        FNFT_COMPLEX x2 = x * x;
        for (FNFT_UINT n=9; n-->0; )
            returnvalue = (( n%2 ? 1.0 : -1.0) + x2/(2*n+2) * returnvalue ) / (2*n+3);
        returnvalue *= x;
    }
    return returnvalue;
}

/**
 * @brief Closest larger or equal number that is a power of two.
 *
 * @ingroup misc
 * @param [in] number
 * @return min{r >= number : exists d such that r = 2^d}
 */
FNFT_UINT fnft__misc_nextpowerof2(const FNFT_UINT number);

/**
 * @brief Resamples an array.
 *
 * @ingroup misc
 * Computes a resampled version of q. The length of q is D>=2. Performs
 * bandlimited interpolation to obtain signal samples at new locations
 * \f$q(t_{n}+\delta)\f$ from samples given at \f$q(t_{n})\f$ for \f$n=0,1,...,D-1\f$.
 * This is required for the CF4_2, CF4_3, CF5_3, CF6_4, 4SPLIT4A and 4SPLIT4B discretizations.
 * The routine checks the difference between the l2-norm of the complete spectrum 
 * and approximately 90% of the spectrum. The following warning is issued if
 * the difference is high. "Signal does not appear to be bandlimited. Interpolation 
 * step may be inaccurate. Try to reduce the step size, or switch to a discretization 
 * that does not require interpolation".
 * @param[in] D Number of samples in array q.
 * @param[in] eps_t Step-size of t.
 * @param[in] q Complex valued array to be resampled.
 * @param[in] delta Real valued shift to be applied during resampling.
 * @param[out] q_new Complex valued array of resampled signal.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_resample(const FNFT_UINT D, const FNFT_REAL eps_t, FNFT_COMPLEX const * const q,
    const FNFT_REAL delta, FNFT_COMPLEX *const q_new);

/**
 * @brief Multiples two complex valued matrices of any compatible size.
 *
 * @ingroup  misc
 *
 * Right-multiples an nxm matrix A by an mxp matrix B. The resulting nxp matrix
 * is stored in C.
 * @param[in] n Positive integer that is the number of rows of A and C.
 * @param[in] m Positive integer that is the number of columns of A and the
 * number of rows of B.
 * @param[in] p Positive integer that is the number of columns of B and C.
 * @param[in] A Pointer to the first element of complex valued nxm matrix A.
 * @param[in] B Pointer to the first element of complex valued mxp matrix B.
 * @param[out] C Pointer to the first element of complex valued nxp matrix C
 * which will contain the result C=AB on return. The user should allocate
 * memory for C, different from the memory that contains A and B.
 *
 */
static inline void fnft__misc_matrix_mult(const FNFT_UINT n,
                                          const FNFT_UINT m,
                                          const FNFT_UINT p,
                                          FNFT_COMPLEX const * const A,
                                          FNFT_COMPLEX const * const B,
                                          FNFT_COMPLEX * const C){
    for (FNFT_UINT cn = 0; cn < n; cn++) {
        for (FNFT_UINT cp = 0; cp < p; cp++) {
            C[ cn*p + cp] = 0.0;
            for (FNFT_UINT cm = 0; cm < m; cm++)
                C[ cn*p + cp ] += A[ cn*m + cm ] * B[ cm*p + cp ];
        }
    }
    
    return;
}

/**
 * @brief Multiples two square matrices of size 2.
 *
 * @ingroup  misc
 *
 * Multiples two square matrices U and T of size 2. T is replaced by the
 * result U*T.
 * @param[in] U Pointer to the first element of complex values 2x2 matrix U.
 * @param[in,out] T Pointer to the first element of complex values 2x2 matrix T.
 * Contains the result U*T on return.
 */
static inline void fnft__misc_mat_mult_2x2(FNFT_COMPLEX * const U,
        FNFT_COMPLEX *const T){

    FNFT_COMPLEX TM[4] = { 0 };
    fnft__misc_matrix_mult(2, 2, 2, U, T, &TM[0]);
    memcpy(T,TM,4*sizeof(FNFT_COMPLEX));
    
    return;
}

/**
 * @brief Multiples two square matrices of size 4.
 *
 * @ingroup  misc
 *
 * Multiples two square matrices U and T of size 4. T is replaced by the
 * result U*T.
 * @param[in] U Pointer to the first element of complex values 4x4 matrix U.
 * @param[in,out] T Pointer to the first element of complex values 4x4 matrix T.
 * Contains the result U*T on return.
 */
static inline void fnft__misc_mat_mult_4x4(FNFT_COMPLEX * const U,
        FNFT_COMPLEX *const T){

    FNFT_COMPLEX TM[16] = { 0 };
    fnft__misc_matrix_mult(4, 4, 4, U, T, &TM[0]);
    memcpy(T,TM,16*sizeof(FNFT_COMPLEX));
    
    return;
}

/**
 * @brief This routine returns the nth degree Legendre polynomial at x.
 *
 * @ingroup  misc
 *
 * Calculates the the nth degree Legendre polynomial at x using a recursive 
 * relation (<a href="https://en.wikipedia.org/wiki/Legendre_polynomials#Definition_via_generating_function">Online, Accessed July 2020</a>)
 * @param[in] n Positive integer that is the order of the Legendre polynomial.
 * @param[in] x Real scalar value at which the value of the polynomial is to be calculated.
 * @return Returns the value of nth degree Legendre polynomial at x.
 */
static inline FNFT_REAL fnft__misc_legendre_poly(const FNFT_UINT n, const FNFT_REAL x){
    FNFT_UINT  i;
    FNFT_REAL P, P_1, P_2;
    if (n == 0)
        P = 1;
    else if (n == 1)
        P = x;
    else{
        P_1 = x;
        P_2 = 1;
        for (i = 2; i <= n; i++) {
            P = (2.0*i-1)*x*P_1/i -(i-1.0)*P_2/i;
            P_2 = P_1;
            P_1 = P;
        }
    }
    return P;
}

/**
 * @brief Minimum number in an array
 *
 * @ingroup misc
 * This finds the smallest element in the vector vecA.
 * @param[in] lenA Length of vector vecA.
 * @param[in] vecA Real vector of length lenA.
 * @return Returns the smallest element in the vector vecA.
 */
FNFT_REAL fnft__misc_min(const FNFT_UINT lenA,
    FNFT_REAL const * const vecA);

/**
 * @brief Maximum number in an array
 *
 * @ingroup misc
 * This finds the largest element in the vector vecA.
 * @param[in] lenA Length of vector vecA.
 * @param[in] vecA Real vector of length lenA.
 * @return Returns the largest element in the vector vecA.
 */
FNFT_REAL fnft__misc_max(const FNFT_UINT lenA,
    FNFT_REAL const * const vecA);

/**
 * @brief This routine calcuates the argument and quadrant of complex values.
 *
 * @ingroup  misc
 *
 * Calculates the argument from 0 to 2*Pi and the quadrant a complex point
 * lies in.
 * @param[in] N Length of vals array.
 * @param[in] vals Complex array.
 * @param[out] quadrants Integer valued array. If the complex point is in
 * first quadrant then it will be 1. It will be 2 for second qudrant and so on.
 * @param[out] arguments Real valued array containing the argument of the 
 * complex value between 0 to 2*PI rather than the usual -PI to PI.
 * @return Returns SUCCESS or an error code.
 */
FNFT_INT fnft__misc_quadrant(FNFT_UINT N, FNFT_COMPLEX * const vals, 
        FNFT_UINT * const quadrants, FNFT_REAL * const arguments);

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
#define misc_CSINC_derivative(...) fnft__misc_CSINC_derivative(__VA_ARGS__)
#define misc_nextpowerof2(...) fnft__misc_nextpowerof2(__VA_ARGS__)
#define misc_resample(...) fnft__misc_resample(__VA_ARGS__)
#define misc_matrix_mult(...) fnft__misc_matrix_mult(__VA_ARGS__)
#define misc_mat_mult_2x2(...) fnft__misc_mat_mult_2x2(__VA_ARGS__)
#define misc_mat_mult_4x4(...) fnft__misc_mat_mult_4x4(__VA_ARGS__)
#define misc_legendre_poly(...) fnft__misc_legendre_poly(__VA_ARGS__)
#define misc_quadrant(...) fnft__misc_quadrant(__VA_ARGS__)
#define misc_min(...) fnft__misc_min(__VA_ARGS__)
#define misc_max(...) fnft__misc_max(__VA_ARGS__)

#endif

#endif
