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
 * Sander Wahls (TU Delft) 2018.
 */

/**
 * @file fnft_nsev_inverse.h
 * @brief Fast inverse nonlinear Fourier transform for the vanishing
 *  nonlinear Schroedinger equation.
 * @ingroup fnft_inverse
 */

#ifndef FNFT_NSEV_INVERSE_H
#define FNFT_NSEV_INVERSE_H

#include "fnft.h"
#include "fnft_nse_discretization_t.h"

/**
 * Enum that specifies in which form the continuous spectrum is provided.
 * Used in \link fnft_nsev_inverse_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT: The array contspec contains
 *  samples of \f$ b(\xi)/a(\xi) \f$ on the grid specified through the input
 *  XI to \link fnft_nsev_inverse \endlink. \n \n
 *  fnft_nsev_inverse_cstype_B_OF_TAU: The array contspec contains samples of the
 *  inverse Fourier transform \f${ B(\tau) = \frac{1}{2\pi}
 *  \int_{-\infty}^\infty b(\xi) e^{j \xi \tau} d\tau }\f$ of \f$ b(\xi) \f$ at
 *  the locations \f$ \tau_n = 2t_n \f$, \f$ n=0,1,\dots,D-1 \f$, where the
 *  \f$ t_n \f$ are the locations at which \f$ q(t) \f$ is computed. The default
 *  (and currently only implemented) method for this type of spectrum is described
 *  in <a href="https://doi.org/10.1109/ECOC.2017.8346231">[Wahls 2017]</a>.
 *  It is currently REQUIRED that the time window is symmetric, T[0]=-T[1].
 *
 */
typedef enum {
    fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT,
    fnft_nsev_inverse_cstype_B_OF_TAU
} fnft_nsev_inverse_cstype_t;

/**
 * Enum that specifies which algorithm is used to inverst the continuous
 * spectrum. Used in \link fnft_nsev_inverse_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsev_inverse_csmethod_DEFAULT: \link fnft_nsev_inverse \endlink
 *  chooses a default method based on the type of spectrum and wether
 *  we are in the defocusing case or not.\n\n
 *  fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_REFL_COEFF: This is
 *  essentially the algorithm in Section II of
 *  <a href="https://doi.org/10.1109/3.903065">[Skaar et al, J Quantum
 *  Electron 37(2), 2001]</a>, only that we use fast inverse scattering as
 *  in <a href="https://doi.org/10.1109/ISIT.2015.7282741">[Wahls & Poor,
 *  Proc. IEEE ISIT 2015]</a>, <a href="https://doi.org/10.1190/1.1441417">
 *  [McClary, Geophysics 48(10), 1983] </a> instead of conventional layer
 *  peeling.\n\n
 *  fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER: This is
 *  Algorithm 1 in the unpublished preprint http://arxiv.org/abs/1607.01305v2 .
 *  The maximum number of iterations can be controlled using the max_iter
 *  field in \link fnft_nsev_inverse_opts_t \endlink. Requires M=D. Defocusing
 *  case only.
 */
typedef enum {
    fnft_nsev_inverse_csmethod_DEFAULT,
    fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_REFL_COEFF,
    fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER
} fnft_nsev_inverse_csmethod_t;

/**
 * @struct fnft_nsev_inverse_opts_t
 * @brief Stores additional options for the routine
 * \link fnft_nsev_inverse \endlink.
 * @ingroup fnft_inverse
 * @ingroup data_types
 *
 * Use the \link fnft_nsev_inverse_default_opts \endlink routine in order
 * to generate a new variable of this type with default options and modify
 * as needed.
 *
 * @var fnft_nsev_inverse_opts_t::discretization
 *  Controls which discretization is applied to the continuous-time Zakharov-
 *  Shabat scattering problem. See \link fnft_nse_discretization_t \endlink.
 *  Currently, only 2SPLIT2A and 2SPLIT2_MODAL are supported.
 *
 * @var fnft_nsev_inverse_opts_t::contspec_type
 *  Controls how \link fnft_nsev_inverse \endlink interprets the values in
 *  the array contspec. \n
 *  Should be of type \link fnft_nsev_inverse_cstype_t \endlink.
 *
 * @var fnft_nsev_inverse_opts_t::contspec_inversion_method
 *  Determines which algorithm \link fnft_nsev_inverse \endlink uses to
 *  invert the continuous spectrum. \n
 *  Should be of type \link fnft_nsev_inverse_csmethod_t \endlink
 *
 * @var fnft_nsev_inverse_opts_t::max_iter
 *  Determines the maximum number of iterations in iterative methods
 *  for the inversion of the continuous spectrum.
 *
 * @var fnft_nsev_inverse_opts_t::oversampling_factor
 *  The oversampling factor used when calling
 * \link fnft__poly_specfact \endlink internally. (Not all methods
 * compute a spectral factorization). Should be at least one.
 */
typedef struct {
    fnft_nse_discretization_t discretization;
    fnft_nsev_inverse_cstype_t contspec_type;
    fnft_nsev_inverse_csmethod_t contspec_inversion_method;
    FNFT_UINT max_iter;
    FNFT_UINT oversampling_factor;
} fnft_nsev_inverse_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_nsev_inverse \endlink
 * with default settings.
 *
 * @returns A \link fnft_nsev_inverse_opts_t \endlink object with the
 *  following options.\n\n
 *  discretization = fnft_nse_discretization_2SPLIT2A\n\n
 *  contspec_type = fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT\n\n
 *  contspec_inverse_method = fnft_nsev_inverse_csmethod_DEFAULT
 *
 * @ingroup fnft_inverse
 */
fnft_nsev_inverse_opts_t fnft_nsev_inverse_default_opts();

/**
 * @brief Determines the locations of the first and last sample of
 * the grid that has to be used when providing a continuous spectrum
 * to \link fnft_nsev_inverse \endlink.
 *
 * @param[in] D Number of samples of the to be generated signal q.
 * @param[in] T Array of length two. Contains the desired location
 *  of the first and last sample of the signal q.
 * @param[in] M Number of samples in the nonlinear frequency domain.
 * @param[out] XI Array of length 2. Will be overwritten with the
 *  desired values specifying the grid for the continuous spectrum.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup fnft_inverse
 */
FNFT_INT fnft_nsev_inverse_XI(
    const FNFT_UINT D,
    FNFT_REAL const * const T,
    const FNFT_UINT M,
    FNFT_REAL * const XI,
    const fnft_nse_discretization_t discretization);

/**
 * @brief Fast inverse nonlinear Fourier transform for the nonlinear
 *  Schroedinger equation with vanishing boundary conditions.
 *
 * This routine inverts the nonlinear Fourier transform for the nonlinear
 * Schroedinger equation \f[ iq_x + q_{tt} \pm 2q|q|^2=0, \quad  q=q(x,t), \f]
 * of Zakharov and Shabat (<a href="http://jetp.ac.ru/cgi-bin/e/index/e/34/1/p62?a=list">Soviet. Phys. JTEP 31(1), 1972</a>)
 * for initial conditions with vanishing boundaries
 * \f[ \lim_{t\to \pm \infty }q(x_0,t) = 0 \text{ sufficiently rapidly.} \f]
 *
 * This routine is the counter part to \link fnft_nsev \endlink. The inversion
 * of discrete specturm is not yet supported. The main references are\n
 *
 * - McClary, <a href="https://doi.org/10.1190/1.1441417">&quot;Fast
 *   Seismic Inversion&quot;</a>, Geophysis 48(10), 1983.
 * - Skaar et al, <a href="">&quot; On the synthesis of fiber Bragg
 *   gratings by layer peeling&quot;</a>, IEEE J Quantum Electon.
 *   37(2), 2001.
 * - Wahls and Poor, <a href="https://doi.org/10.1109/ISIT.2015.7282741">
 *   &quot;Fast Inverse Nonlinear Fourier Transform For Generating
 *   Multi-Solitons In Optical Fiber&quot;</a>, Proc. IEEE ISIT 2015.
 * - Wahls and Vaibhav, <a href="https://arxiv.org/abs/1607.01305v2">
 *   &quot;Fast Inverse Nonlinear Fourier Transforms for Continuous
 *   Spectra of Zakharov-Shabat Type&quot;</a>, Unpublished Preprint,
 *   arXiv:1607.01305v2 [cs.IT], Dec. 2016.
 * - Wahls, <a href="https://doi.org/10.1109/ECOC.2017.8346231">&quot;Generation of Time-Limited Signals in the Nonlinear Fourier Domain via b-Modulation&quot;,</a> Proc. ECOC 2017.
 *
 * @param[in] M Number of samples of the continuous spectrum.
 * @param[in,out] contspec Array of length M, contains samples
 *  \f$ \hat{q}(\xi_n) \f$, where \f$ \xi_n = XI[0] + n(XI[1]-XI[0])/(M-1) \f$
 *  and \f$n=0,1,\dots,M-1\f$, of the to-be-inverted continuous spectrum in
 *  ascending order (i.e.,
 *  \f$ \hat{q}(\xi_0), \hat{q}(\xi_1), \dots, \hat{q}(\xi_{M-1}) \f$).
 *  Please note that the routine currently might overwrite this array.
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum. Currently, the positions returned
 *  by \link fnft_nsev_inverse_XI \endlink MUST be used. It is NOT checked
 *  whether the user adheres to this requirement.
 * @param[in] K Ignored as inversion of discrete spectrum is not yet implemented.
 * @param[in] bound_states Ignored as inversion of discrete spectrum is not yet
 *  implemented.
 * @param[in] normconsts_or_residues Ignored as inversion of discrete spectrum is
 *  not yet implemented.
 * @param[in] D Number of samples of the to be generated signal q.
 * @param[out] q Array of length D. Is filled with samples
 *  \f$ q(t_n) \f$, where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$
 *  and \f$n=0,1,\dots,D-1\f$, of the to-be-generated signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$).
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample of q. It should be T[0]<T[1].
 * @param[in] kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one
 * @param[in] opts_ptr Pointer to a \link fnft_nsev_inverse_opts_t \endlink
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
FNFT_INT fnft_nsev_inverse(
    const FNFT_UINT M,
    FNFT_COMPLEX * const contspec,
    FNFT_REAL const * const XI,
    FNFT_UINT const K,
    FNFT_COMPLEX const * const bound_states,
    FNFT_COMPLEX const * const normconsts_or_residues,
    const FNFT_UINT D,
    FNFT_COMPLEX * const q,
    FNFT_REAL const * const T,
    const FNFT_INT kappa,
    fnft_nsev_inverse_opts_t *opts_ptr);

#endif
