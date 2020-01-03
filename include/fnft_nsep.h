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
* Sander Wahls (TU Delft) 2017-2018, 2020.
*/

/**
 * @file fnft_nsep.h
 * @brief Fast nonlinear Fourier transform for the periodic nonlinear
 *  Schroedinger equation.
 * @ingroup fnft
 */

#ifndef FNFT_NSEP_H
#define FNFT_NSEP_H

#include "fnft_nse_discretization_t.h"

/**
 * Enum that controls how spectrum is localized. Used in
 * \link fnft_nsep_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsep_opts_loc_SUBSAMPLE_AND_REFINE: Similar approach as for
 *  fnft_nsev_opts_dsloc_SUBSAMPLE_AND_REFINE (see
 *  \link fnft_nsev_opts_t::bound_state_localization\endlink.)\n\n
 *  fnft_nsep_opts_loc_GRIDSEARCH: Uses a grid search to localize roots. Can
 *  only find main and auxiliary spectrum points on the real axis. In the
 *  defocusing case, the main spectrum is always real.\n\n
 *  fnft_nsep_opts_loc_MIXED: Uses the SUBSAMPLE_AND_REFINE method to find the
 *  non-real parts of the spectra and the GRIDSEARCH method to find the real
 *  parts.
 */
typedef enum {
    fnft_nsep_loc_SUBSAMPLE_AND_REFINE,
    fnft_nsep_loc_GRIDSEARCH,
    fnft_nsep_loc_MIXED
} fnft_nsep_loc_t;

/**
 * Enum that controls how spectrum is filtered. Used in
 * \link fnft_nsep_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsep_opts_filt_NONE: No filtering.\n\n
 *  fnft_nsep_opts_filt_MANUAL: Only points within the specified
 *  \link fnft_nsep_opts_t::bounding_box \endlink are kept. \n\n
 *  fnft_nsep_opts_filt_AUTO: As above, but the boundary box is chosen by
 *  the routine.
 */
typedef enum {
    fnft_nsep_filt_NONE,
    fnft_nsep_filt_MANUAL,
    fnft_nsep_filt_AUTO
} fnft_nsep_filt_t;

/**
 * @struct fnft_nsep_opts_t
 * @brief Stores additional options for the routine \link fnft_nsep \endlink.
 *
 * @ingroup fnft
 * @ingroup data_types
 *
 * Use the \link fnft_nsep_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_nsep_opts_t::localization
 *  Controls how both the main spectrum and the auxiliary main spectrum points
 *  are localized. \n
 *  Should be of type \link fnft_nsep_loc_t \endlink.
 *
 * @var fnft_nsep_opts_t::filtering
 *  Controls how both the main spectrum and the auxiliary spectrum are
 *  filtered. It is recommended to use MANUAL filtering as accurate filtering
 *  reduces the runtime of the routine.\n
 *  Should be of type \link fnft_nsep_filt_t \endlink.
 *
 * @var fnft_nsep_opts_t::bounding_box
 *  Array of four reals. Defines a box in the complex plane that is used for
 *  filtering: \n
 *  bounding_box[0] <= real(lambda) <= bounding_box[1] \n
 *  bounding_box[2] <= imag(lambda) <= bounding_box[3] \n
 *
 * @var fnft_nsep_opts_t::max_evals
 *  Maximum number of function evaluations per root allowed for the refinement
 *  step during localization. Note that this is not necessarily the number
 *  of Newton iterations as Newton's method for higher order roots is used for
 *  some parts of the spectra.
 *
 * @var fnft_nsep_opts_t::discretization
 *  See \link fnft_nsev_opts_t::discretization \endlink.
 *
 * @var fnft_nsep_opts_t::normalization_flag
 *  See \link fnft_nsev_opts_t::normalization_flag \endlink.
 *
 * @var fnft_nsep_opts_t::floquet_range
 *   Array of two reals. The mainspec variable will contain the z that solve
 *   Delta(z)=rhs, where Delta(z)=0.5 trace{monodromy matrix(z)} and rhs is
 *   varied over a equidistant grid of
 *   \link fnft_nsep_opts_t::floquet_nvals \endlink values with the first grid
 *   point being floquet_range[0] and last grid point being floquet_range[1].
 *   By default, floquet_range={-1,1} and floquet_nvals=2, which corresponds to
 *   the conventional main spectrum. By choosing larger nvals, one can determine
 *   the spines (or bands) connecting the points in the main spectrum.
 *
 * @var fnft_nsep_opts_t::floquet_nvals
 *  See \link fnft_nsep_opts_t::floquet_range \endlink.
 */
typedef struct {
    fnft_nsep_loc_t localization;
    fnft_nsep_filt_t filtering;
    FNFT_REAL bounding_box[4];
    FNFT_UINT max_evals;
    fnft_nse_discretization_t discretization;
    FNFT_INT normalization_flag;
    FNFT_REAL floquet_range[2];
    FNFT_UINT floquet_nvals;
} fnft_nsep_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_nsep \endlink with
 * default settings.
 *
 * @ingroup fnft
 * @returns A \link fnft_nsep_opts_t \endlink object with the following options.\n
 *  localization = fnft_nsep_loc_MIXED\n
 *  filtering = fnft_nsep_filt_AUTO\n
 *  max_evals = 20\n
 *  bounding_box[0] = -FNFT_INF\n
 *  bounding_box[1] = FNFT_INF\n
 *  bounding_box[2] = -FNFT_INF\n
 *  bounding_box[3] = FNFT_INF\n
 *  normalization_flag = 1\n
 *  discretization = fnft_nse_discretization_2SPLIT2A\n
 */
fnft_nsep_opts_t fnft_nsep_default_opts();

/**
 * @brief Fast nonlinear Fourier transform for the nonlinear Schroedinger
 *  equation with periodic boundary conditions.
 *
 * @ingroup fnft
 * \n \n
 * This routine computes the nonlinear Fourier transform for the nonlinear
 * Schroedinger equation \f[ iq_x + q_{tt} \pm 2q|q|^2=0, \quad  q=q(x,t), \f]
 * of Kotlyarov and Its (<a href="https://arxiv.org/abs/1401.4445">Problems
 * of Mathemetical Physics and Functional Analysis, 1976)</a> for initial
 * conditions with periodic boundary conditions,
 * \f[ q(x_0,t + L) = q(x_0,t). \f]
 * \n
 * The main references for the numerical algorithm are:
 *      - Wahls and Poor, <a href="http://dx.doi.org/10.1109/TIT.2015.2485944">&quot;Fast numerical nonlinear Fourier transforms,&quot;</a> IEEE Trans. Inform. Theor. 61(12), 2015.
 *      - Prins and Wahls, &quot;Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot; to appear in Proc. ICASSP 2018.
 *
 * @param[in] D Number of samples. Has to be a power of two.
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n*L/D \f$, where L=T[2]-T[1] is the period and
 *  \f$n=0,1,\dots,D-1\f$, of the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] T Array of length 2. T[0] is the position in time of the first
 *  sample. T[2] is the beginning of the next period. (The location of the last
 *  sample is thus t_{D-1}=T[2]-L/D.) It should be T[0]<T[1].
 * @param[in,out] K_ptr Upon entry, *K_ptr should contain the length of the array
 *  main_spec. Upon return, *K_ptr contains the number of actually detected
 *  points in the main spectrum. If the length of the array main_spec was not
 *  sufficient to store all of the detected points in the main spectrum, a
 *  warning is printed and as many points as possible are returned instead.
 *  Note that in order to skip the computation of the main spectrum completely,
 *  it is not sufficient to pass *K_ptr==NULL. Instead, pass main_spec==NULL.
 * @param[out] main_spec Array. Upon return, the routine has stored the detected
 *  main specrum points (i.e., the points for which the trace of the monodromy
 *  matrix is either +2 or -2) in the first *K_ptr entries of this array.
 *  If NULL is passed instead, the main spectrum will not be computed.
 *  Has to be preallocated by the user. The user can choose an arbitrary
 *  length. Typically, D is a good choice.
 * @param[in,out] M_ptr Upon entry, *M_ptr should contain the length of the array
 *  aux_spec. Upon return, *M_ptr contains the number of actually detected
 *  points in the auxiliary spectrum. If the length of the array aux_spec was
 *  not sufficient to store all of the detected points in the auxiliary
 *  spectrum, a warning is printed and as many points as possible are returned
 *  instead. Note that in order to skip the computation of the auxiliary
 *  spectrum completely, it is not sufficient to pass *M_ptr==NULL. Instead,
 *  pass aux_spec==NULL.
 * @param[out] aux_spec Array. Upon return, the routine has stored the detected
 *  auxiliary specrum points (i.e., the points for which the upper right
 *  element of the monodromy matrix is zero) in the first *M_ptr entries of
 *  this array. If NULL is passed instead, the auxiliary spectrum will not be
 *  computed. Has to be preallocated by the user. The user can choose an
 *  arbitrary length. Typically, D is a good choice.
 * @param[in] sheet_indices Not yet implemented. Pass NULL.
 * @param[in] kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one
 * @param[in] opts Pointer to a \link fnft_nsep_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_nsep_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft_nsep(const FNFT_UINT D, FNFT_COMPLEX const * const q,
    FNFT_REAL const * const T, FNFT_UINT * const K_ptr,
    FNFT_COMPLEX * const main_spec, FNFT_UINT * const M_ptr,
    FNFT_COMPLEX * const aux_spec, FNFT_REAL * const sheet_indices,
    const FNFT_INT kappa, fnft_nsep_opts_t * opts);

#endif
