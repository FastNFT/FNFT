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
 * Sander Wahls (KIT) 2026
 * Fabian Fischer (Hiwi KIT) 2026
 */

/**
 * @file fnft_kdvv_inverse.h
 * @brief Fast inverse nonlinear Fourier transform for the vanishing
 *  Korteweg-de Vries equation.
 * @ingroup fnft_inverse
 */

#ifndef FNFT_KDVV_INVERSE_H
#define FNFT_KDVV_INVERSE_H

#include "fnft__errwarn.h"
#include "fnft__misc.h"
#include "fnft_numtypes.h"

/**
 * @brief Fast inverse nonlinear Fourier transform for the 
 *  Korteweg-de Vries equation with vanishing boundary conditions.
 * 
 * This routine computes the inverse nonlinear Fourier transform for the
 * Korteweg-de Vries equation
 * \f[ q_x + 6qq_{t} + q_{ttt}=0, \quad  q=q(x,t), \f]
 * of Gardner et al. (<a href="https://doi.org/10.1103/PhysRevLett.19.1095">
 * Phys. Rev. Lett., 1967</a>)
 * 
 * The Fast inverse nonlinear Fourier transform for the Korteweg-de Vries equation
 * uses the crum transformation. The main references are:
 *      - Prins and Wahls, <a href="https://doi.org/10.1016/j.cnsns.2021.105782">&quot;An accurate O(N^2) floating point algorithm for the Crum transform of the KdV equation,&quot;</a> Communications in Nonlinear Science and Numerical Simulation 102, Article 105782, 2021.
 *      - The matlab project of P. Prins to the Crum Transformation
 *
 * @param[in] M Number of samples of the continuous spectrum.
 * @param[in,out] contspec Array of length M, contains samples
 *  \f$ \hat{q}(\xi_n) \f$, where \f$ \xi_n = XI[0] + n(XI[1]-XI[0])/(M-1) \f$
 *  and \f$n=0,1,\dots,M-1\f$, of the to-be-inverted continuous spectrum in
 *  ascending order (i.e.,
 *  \f$ \hat{q}(\xi_0), \hat{q}(\xi_1), \dots, \hat{q}(\xi_{M-1}) \f$).
 *  Note: contspec functionality currently out of function (state 04/2026)! 
 *  It will be implicitly assumed to \f$ 0 \f$ for all \f$ \xi \f$.
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum.
 *  Note: contspec functionality currently out of function (state 04/2026)! It is 
 *  equal which values are chosen at the moment.
 * @param[in] K Number of discrete spectrum points.
 * @param[in] bound_states Complex array of length K. Bound states have to be positive, 
 *  purely imaginary numbers (lie on the upper half of the imaginary axis). The bound 
 *  states have to be in descending order. To add a solition with height \f$ h_i \f$ 
 *  the bound state have to be \f$ \gamma_i = \sqrt{ h_i/2} \f$. 
 * @param[in] normconsts_or_residues Complex array of length K. Values of
 *  either the norming constants \f$ b(\xi) \f$ or the residues
 *  \f$ \frac{b(\xi)}{\partial{a(\xi)}/\partial{\xi}}\f$ at the values bound_states.
 *  The signs of the norming constants have to alternate regards to the order of the 
 *  bound states. The sign of the normconst for the biggest eigenvalue has to be positive
 *  Note: currently this array is always interpreted as norming constants (state 04/2026)! 
 *  Residues functionality is not implemented yet!
 * @param[in] D Number of samples of the to be generated signal q. Should be a
 *  positive power of two.
 * @param[out] q Array of length D. Is filled with samples
 *  \f$ q(t_n) \f$, where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$
 *  and \f$n=0,1,\dots,D-1\f$, of the to-be-generated signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$).
 *  Has to be preallocated by the user.
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample of q. It should be \f$ T[0]<T[1] \f$.
 * @param[in] opts_ptr Pointer to a \link fnft_kdvv_inverse_opts_t \endlink
 *  object. The object  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_kdvv_inverse_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
 *  Note: opts_ptr is currently out of function (state 04/2026)!
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft_inverse
 */
FNFT_INT fnft_kdvv_inverse(
    const FNFT_UINT M,
    FNFT_COMPLEX * const contspec,
    FNFT_REAL const * const XI,
    FNFT_UINT const K,
    FNFT_COMPLEX const * const bound_states,
    FNFT_COMPLEX const * const normconsts_or_residues,
    const FNFT_UINT D,
    FNFT_COMPLEX * const q,
    FNFT_REAL const * const T,
    void *opts_ptr);

#endif
