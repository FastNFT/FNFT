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

/**
 * @brief Fast inverse nonlinear Fourier transform for the 
 *  Korteweg-de Vries equation with vanishing boundary conditions.
 *
 * @param[in] M Number of samples of the continuous spectrum.
 * @param[in,out] contspec Array of length M, contains samples
 *  \f$ \hat{q}(\xi_n) \f$, where \f$ \xi_n = XI[0] + n(XI[1]-XI[0])/(M-1) \f$
 *  and \f$n=0,1,\dots,M-1\f$, of the to-be-inverted continuous spectrum in
 *  ascending order (i.e.,
 *  \f$ \hat{q}(\xi_0), \hat{q}(\xi_1), \dots, \hat{q}(\xi_{M-1}) \f$).
 *  Please note that the routine currently might overwrite this array.
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum.
 * @param[in] K Number of discrete spectrum points.
 * @param[in] bound_states Complex array of length K. Complex roots of
 *  \f$ a(\xi) \f$ in the upper half of the complex-plane.
 * @param[in] normconsts_or_residues Complex array of length K. Values of
 *  either the norming constants \f$ b(\xi) \f$ or the residues
 *  \f$ \frac{b(\xi)}{\partial{a(\xi)}/\partial{\xi}}\f$ at the values bound_states.
 * @param[in] D Number of samples of the to be generated signal q. Should be a
 *  positive power of two.
 * @param[out] q Array of length D. Is filled with samples
 *  \f$ q(t_n) \f$, where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$
 *  and \f$n=0,1,\dots,D-1\f$, of the to-be-generated signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$).
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample of q. It should be T[0]<T[1].
 * @param[in] opts_ptr Pointer to a \link fnft_kdvv_inverse_opts_t \endlink
 *  object. The object  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_nsev_inverse_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
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
