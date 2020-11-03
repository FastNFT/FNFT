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
* Peter J Prins (TU Delft) 2017-2018.
*/

/**
 * @file fnft_kdvv.h
 * @brief Fast nonlinear Fourier transform for the vanishing Korteweg-de Vries
 * equation.
 *
 * @ingroup fnft
 */

#ifndef FNFT_KDVV_H
#define FNFT_KDVV_H

#include "fnft_kdv_discretization_t.h"

/**
 * @struct fnft_kdvv_opts_t
 * @brief Stores additional options for the routine \link fnft_kdvv \endlink. 
 * @ingroup fnft
 * @ingroup data_types
 * Use the \link fnft_kdvv_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_kdvv_opts_t::discretization
 *  Controls which discretization is applied to the continuous-time scattering
 *  problem. See \link fnft_kdv_discretization_t \endlink.
 */
typedef struct {
    fnft_kdv_discretization_t discretization;
} fnft_kdvv_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_kdvv \endlink with
 * default settings.
 *
 * @returns A \link fnft_kdvv_opts_t \endlink object with default options.
 *
 * @ingroup fnft
 */
fnft_kdvv_opts_t fnft_kdvv_default_opts();

/**
 * @brief Fast nonlinear Fourier transform for the Korteweg-de Vries
 * equation with vanishing boundary conditions.
 *
 * This routine computes the nonlinear Fourier transform for the
 * Korteweg-de Vries equation
 * \f[ q_x + 6qq_{t} + q_{ttt}=0, \quad  q=q(x,t), \f]
 * of Gardner et al. (<a href="https://doi.org/10.1103/PhysRevLett.19.1095">
 * Phys. Rev. Lett., 1967</a>)
 * for initial conditions with vanishing boundaries
 * \f[ \lim_{t\to \pm \infty }q(x_0,t) = 0 \text{ sufficiently rapidly.} \f]
 *
 * @param[in] D Number of samples.
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 * where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 * the to-be-transformed signal in ascending order
 * (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$).
 * @param[in] T Array of length 2, contains the position in time of the first
 * and of the last sample. It should be T[0]<T[1].
 * @param[in] M Number of points at which the continuous spectrum (aka
 * reflection coefficient) should be computed.
 * @param[out] contspec Array of length M in which the routine will store the
 * desired samples \f$ R(\xi_m) \f$ of the continuous spectrum (aka
 * reflection coefficient) in ascending order,
 * where \f$ \xi_m = XI[0]+m(XI[1]-XI[0])/(M-1) \f$ and \f$m=0,1,\dots,M-1\f$.
 * Has to be preallocated by the user.
 * @param[in] XI Array of length 2, contains the position of the first and the
 * last sample of the continuous spectrum. It should be XI[0]<XI[1]. Can also be
 * NULL if contspec==NULL.
 * @param[in,out] K_ptr Not yet implemented. Pass NULL.
 * @param[out] bound_states Not yet implemented. Pass NULL.
 * @param[out] normconsts_or_residues Not yet implemented. Pass NULL.
 * @param[in] opts_ptr Pointer to a \link fnft_kdvv_opts_t \endlink object. The
 * object can be used to modify the behavior of the routine. Use
 * the routine \link fnft_kdvv_default_opts \endlink
 * to generate such an object and modify as desired. It is also possible to
 * pass NULL, in which case the routine will use the default options. The
 * user is reponsible to freeing the object after the routine has returned.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */

FNFT_INT fnft_kdvv(const FNFT_UINT D, FNFT_COMPLEX * const q, 
    FNFT_REAL const * const T, const FNFT_UINT M, 
    FNFT_COMPLEX * const contspec, FNFT_REAL const * const XI,
    FNFT_UINT * const K_ptr, FNFT_COMPLEX * const bound_states,
    FNFT_COMPLEX * const normconsts_or_residues,
    fnft_kdvv_opts_t * opts_ptr); 

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdvv_opts_t fnft_kdvv_opts_t
#endif

#endif
