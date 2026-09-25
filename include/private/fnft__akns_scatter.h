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
* Sander Wahls (TU Delft) 2017-2018, 2023.
* Shrinivas Chimmalgi (TU Delft) 2017-2020.
* Peter J. Prins (TU Delft) 2021.
* Sander Wahls (KIT), 2023.
* Igor Chekhovskoy (NSU, FRC ICT) 2026.
*/

/**
 * @file fnft__akns_scatter.h
 * @brief Slow forward scattering.
 * @ingroup akns
 */

#ifndef FNFT__AKNS_SCATTER_H
#define FNFT__AKNS_SCATTER_H

#include "fnft__akns_discretization.h"
#include <float.h>
#include <stdio.h>
#include <string.h> // for memcpy
#include "fnft__errwarn.h"
#include "fnft__misc.h"// for square_matrix_mult
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#include <mathimf.h>
#endif

/* Internal exact ES scalar kernel. It is inline so the production scatterer
 * and its direct numerical oracle exercise one implementation without adding
 * an exported private symbol. */
static inline void fnft__akns_scatter_exact_scalar_pair(
        FNFT_COMPLEX const root, FNFT_COMPLEX * const f,
        FNFT_COMPLEX * const g)
{
    const FNFT_REAL root_real = FNFT_CREAL(root);
    const FNFT_REAL root_imag = FNFT_CIMAG(root);
    const FNFT_REAL safe_real_limit = FNFT_LOG(DBL_MAX)-8.0;

    if (isfinite(root_real) && isfinite(root_imag)
            && FNFT_FABS(root_real) <= safe_real_limit) {
        FNFT_REAL sin_imag, cos_imag, sinh_real, cosh_real;
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
        sincos(root_imag,&sin_imag,&cos_imag);
        sinhcosh(root_real,&sinh_real,&cosh_real);
#else
#if defined(__MINGW32__) && defined(__GNUC__) && !defined(__clang__)
        /* MinGW's paired function uses x87 argument reduction. Restrict
         * it to small arguments and retain libm's full-range functions. */
        if (FNFT_FABS(root_imag) <= 1.0) {
            __builtin_sincos(root_imag,&sin_imag,&cos_imag);
        } else
#endif
        {
            sin_imag = FNFT_SIN(root_imag);
            cos_imag = FNFT_COS(root_imag);
        }
        /* Both hyperbolic functions share one exponential. Using expm1
         * retains relative accuracy near zero; the product below avoids
         * squaring an exponentially large intermediate value. */
        const FNFT_REAL em1 = expm1(FNFT_FABS(root_real));
        const FNFT_REAL ratio = em1/(1.0+em1);
        sinh_real = copysign(0.5*(em1+ratio),root_real);
        cosh_real = 1.0+0.5*em1*ratio;
#endif
        const FNFT_COMPLEX f_candidate = cosh_real*cos_imag
                +(FNFT_COMPLEX)I*(sinh_real*sin_imag);
        FNFT_COMPLEX g_candidate;

        if (root_real == 0.0 && root_imag != 0.0) {
            /* Reuse sin/cos, but retain the signed zeros of
             * sin(i*root)/(i*root) on the imaginary axis. */
            const FNFT_COMPLEX rotated = (FNFT_COMPLEX)I*root;
            const FNFT_REAL parts[2] = {
                -sin_imag,cos_imag*FNFT_CIMAG(rotated)
            };
            FNFT_COMPLEX numerator;
            memcpy(&numerator,parts,sizeof(numerator));
            g_candidate = numerator/rotated;
        } else {
            g_candidate = (sinh_real*cos_imag
                    +(FNFT_COMPLEX)I*(cosh_real*sin_imag))/root;
        }

        if (isfinite(FNFT_CREAL(f_candidate))
                && isfinite(FNFT_CIMAG(f_candidate))
                && isfinite(FNFT_CREAL(g_candidate))
                && isfinite(FNFT_CIMAG(g_candidate))) {
            *f = f_candidate;
            *g = g_candidate;
            return;
        }
    }
    *f = FNFT_CCOSH(root);
    *g = fnft__misc_CSINC((FNFT_COMPLEX)I*root);
}


/**
 * @brief Computes the scattering matrix and its derivative.
 * 
 * The function computes the scattering matrix and the derivative of the scattering matrix with 
 * respect to \f$\lambda\f$. The function performs slow direct scattering and is primarily based on the reference 
 * Boffetta and Osborne 
 * (<a href="http://dx.doi.org/10.1016/0021-9991(92)90370-E">J. Comput. Physics 1992 </a>).
 * 
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q_0, q_1, \dots, q_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in,out] r Array of length D, contains samples \f$ r_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ r_0, r_1, \dots, r_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in] eps_t Step-size, eps_t \f$= (T[1]-T[0])/(D-1) \f$.
 * @param[in] K Number of values of \f$\lambda\f$.
 * @param[in] lambda Array of length K, contains the values of \f$\lambda\f$.
 * @param[out] result Array of length K*8 or K*4,  If derivative_flag=0 returns
 * [S11(lambda[0]) S12(lambda[0]) S21(lambda[0]) S22(lambda[0]) S11(lambda[1]) ... S22(lambda[K-1]] in result where S = [S11, S12; S21, S22] is the
 * scattering matrix computed using the chosen discretization.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in 
 * result where S11' is the derivative of S11 w.r.t to \f$\lambda\f$.
 * Should be preallocated with size 4*K or 8*K accordingly.
 * @param[in,out] W Pass an array of size K. Upon exit, it contains scaling factors that
 * arise due to an internal normalization of the scattering matrix (to deal with potential
 * overflow issues). The result has to be scaled by the power of 2 to obtain the final
 * values. No normalization is carried out if NULL is passed.
 * @param[in] discretization The type of discretization to be used. Should be of type 
 * \link fnft__akns_discretization_t \endlink. Not all akns_discretization_t 
 * discretizations are supported. Check \link fnft__akns_discretization_t \endlink 
 * for list of supported types.
 * @param[in] PDE Partial differential equation of type \link fnft__akns_pde_t \endlink.
 * @param[in] vanilla_flag For calculations for the KdV equation, pass 1 for the original mapping to the AKNS framework with r=-1. Pass 0 for the alternative mapping with q=-1. Unused for NSE.
 * @param[in] derivative_flag Should be set to either 0 or 1. If set to 1
 * the derivatives [S11' S12' S21' S22'] are calculated. result should be
 * preallocated with size 8*K if flag is set to 1.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup akns
 */
FNFT_INT fnft__akns_scatter_matrix(FNFT_UINT const D,
                                   FNFT_COMPLEX const * const q,
                                   FNFT_COMPLEX const * const r,
                                   FNFT_REAL const eps_t,
                                   FNFT_UINT const K,
                                   FNFT_COMPLEX const * const lambda,
                                   FNFT_COMPLEX * const result,
                                   FNFT_INT * const W,
                                   fnft__akns_discretization_t discretization,
                                   fnft__akns_pde_t const PDE,
                                   FNFT_UINT const vanilla_flag,
                                   FNFT_UINT const derivative_flag);

/**
 * @brief Padé variant of \link fnft__akns_scatter_matrix \endlink for
 * ES4/ES6/ES8.
 *
 * @param[in] pade_degree Degree of the diagonal Padé approximant, from 1
 * through 7.
 */
FNFT_INT fnft__akns_scatter_matrix_pade(FNFT_UINT const D,
                                   FNFT_COMPLEX const * const q,
                                   FNFT_COMPLEX const * const r,
                                   FNFT_REAL const eps_t,
                                   FNFT_UINT const K,
                                   FNFT_COMPLEX const * const lambda,
                                   FNFT_COMPLEX * const result,
                                   FNFT_INT * const W,
                                   fnft__akns_discretization_t discretization,
                                   fnft__akns_pde_t const PDE,
                                   FNFT_UINT const vanilla_flag,
                                   FNFT_UINT const pade_degree,
                                   FNFT_UINT const derivative_flag);

/**
 * @brief Computes the first column of a slow ES scattering matrix.
 *
 * This internal routine is restricted to the NSE reduction, real spectral
 * points and ES4/ES6/ES8. The Padé degree is zero for the exact exponential
 * and between one and seven for a diagonal Padé approximation.
 */
FNFT_INT fnft__akns_scatter_matrix_first_column(FNFT_UINT const D,
                                   FNFT_COMPLEX const * const q,
                                   FNFT_COMPLEX const * const r,
                                   FNFT_REAL const eps_t,
                                   FNFT_UINT const K,
                                   FNFT_COMPLEX const * const lambda,
                                   FNFT_COMPLEX * const H11,
                                   FNFT_COMPLEX * const H21,
                                   FNFT_INT * const W,
                                   fnft__akns_discretization_t discretization,
                                   FNFT_UINT const pade_degree);


/**
 * @brief Computes \f$a(\lambda)\f$, \f$ a'(\lambda) = \frac{\partial a(\lambda)}{\partial \lambda}\f$
 * and \f$b(\lambda)\f$ for complex values \f$\lambda\f$ assuming that they are very close to the true
 * bound-states.
 *
 * The function performs slow direct scattering and is primarily based on the references
 *      - Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot; Computation of the direct scattering transform for the nonlinear Schroedinger  equation,&quot;</a> J. Comput. Phys. 102(2), 1992.
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot; Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *      - Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2932256">&quot; Soliton Phase Shift Calculation for the Korteweg–De Vries Equation,&quot;</a> IEEE Access, vol. 7, pp. 122914--122930, July 2019.
 *
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ q_0, q_1, \dots, q_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in] r Array of length D, contains samples \f$ r_n\f$ for \f$n=0,1,\dots,D-1\f$
 * in ascending order (i.e., \f$ r_0, r_1, \dots, r_{D-1} \f$). The values
 * should be specifically precalculated based on the chosen discretization.
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be T[0]<T[1].
 * @param[in] K Number of bound-states.
 * @param[in] bound_states Array of length K, contains the bound-states \f$\lambda\f$.
 * @param[out] a_vals Array of length K, contains the values of \f$a(\lambda)\f$.
 * @param[out] aprime_vals Array of length K, contains the values of
 * \f$ a'(\lambda) = \frac{\partial a(\lambda)}{\partial \lambda}\f$.
 * @param[out] b_vals Array of length K, contains the values of \f$b(\lambda)\f$.
 * The \f$b(\lambda)\f$ are calculated using the criterion from
 * Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2932256">&quot;
 * Soliton Phase Shift Calculation for the Korteweg–De Vries Equation,&quot;</a>.
 * @param[in,out] Ws Pass an array of size K. Upon exit, it contains scaling factors that
 * arise due to an internal normalization of the scattering process (to deal with potential
 * overflow issues). The returned values for a and a_prime still have to be multiplied with
 * corresponding power of two (i.e. POW(2, Ws[i])) to obtain the final values. Note
 * that this is not required for the values of b. No normalization is carried out if NULL
 * is passed.
 * @param[in] discretization The type of discretization to be used. Should be of type
 * \link fnft__akns_discretization_t \endlink. Not all akns_discretization_t discretizations are supported.
 * Check \link fnft_nse_discretization_t \endlink for list of supported types.
 * @param[in] PDE The partial differential equation for which the calculation
 * has to be done. Should be of type \link fnft__akns_pde_t \endlink.
 * @param[in] vanilla_flag For calculations for the KdV equation, pass 1 for the
 * original mapping to the AKNS framework with r=-1. Pass 0 for the alternative
 * mapping with q=-1. Unused for NSE.
 * @param[in] skip_b_flag If set to 1 the routine will not compute \f$b(\lambda)\f$.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup akns
 */
FNFT_INT akns_scatter_bound_states(FNFT_UINT const D,
                                   FNFT_COMPLEX const * const q,
                                   FNFT_COMPLEX const * const r,
                                   FNFT_REAL const *const T,
                                   FNFT_UINT const K,
                                   FNFT_COMPLEX const * const bound_states,
                                   FNFT_COMPLEX * const a_vals,
                                   FNFT_COMPLEX * const aprime_vals,
                                   FNFT_COMPLEX * const b_vals,
                                   FNFT_INT * Ws,
                                   fnft__akns_discretization_t const discretization,
                                   fnft__akns_pde_t const PDE,
                                   FNFT_UINT const vanilla_flag,
                                   FNFT_UINT const skip_b_flag);

/**
 * @brief Padé variant of \link akns_scatter_bound_states \endlink for
 * ES4/ES6/ES8.
 *
 * @param[in] pade_degree Degree of the diagonal Padé approximant, from 1
 * through 7.
 */
FNFT_INT fnft__akns_scatter_bound_states_pade(FNFT_UINT const D,
                                   FNFT_COMPLEX const * const q,
                                   FNFT_COMPLEX const * const r,
                                   FNFT_REAL const *const T,
                                   FNFT_UINT const K,
                                   FNFT_COMPLEX const * const bound_states,
                                   FNFT_COMPLEX * const a_vals,
                                   FNFT_COMPLEX * const aprime_vals,
                                   FNFT_COMPLEX * const b_vals,
                                   FNFT_INT * Ws,
                                   fnft__akns_discretization_t const discretization,
                                   fnft__akns_pde_t const PDE,
                                   FNFT_UINT const vanilla_flag,
                                   FNFT_UINT const pade_degree,
                                   FNFT_UINT const skip_b_flag);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_scatter_matrix(...) fnft__akns_scatter_matrix(__VA_ARGS__)
#define akns_scatter_matrix_pade(...) fnft__akns_scatter_matrix_pade(__VA_ARGS__)
#define akns_scatter_matrix_first_column(...) \
    fnft__akns_scatter_matrix_first_column(__VA_ARGS__)
#define akns_scatter_bound_states_pade(...) fnft__akns_scatter_bound_states_pade(__VA_ARGS__)
#endif

#endif
