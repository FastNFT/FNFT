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
* Peter J Prins (TU Delft) 2020-2021.
* Sander Wahls (KIT) 2023.
*/

/**
 * @file fnft_nsev.h
 * @brief Fast nonlinear Fourier transform for the vanishing nonlinear
 *  Schroedinger equation.
 * @ingroup fnft
 */

#ifndef FNFT_NSEV_H
#define FNFT_NSEV_H

#include "fnft__nse_fscatter.h"
#include "fnft__nse_scatter.h"
#include "fnft__misc.h" // for l2norm
#include <string.h> // for memcpy
#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_chirpz.h"

/**
 * Enum that specifies how the bound states are filtered. Used in
 * \link fnft_nsev_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsev_bsfilt_NONE: All detected roots of \f$ a(\lambda) \f$ are returned. \n \n
 *  fnft_nsev_bsfilt_BASIC: Only roots in the upper halfplane are returned and roots very close
 *  to each other are merged. \n \n
 *  fnft_nsev_bsfilt_FULL: Bound states in physically implausible regions and
 *  outside the region based on the step-size of the supplied samples are
 *  rejected. \n \n
 *  fnft_nsev_opts_filt_MANUAL: Only points within the specified
 *  \link fnft_nsev_opts_t::bounding_box \endlink are kept. \n\n
 */
typedef enum {
    fnft_nsev_bsfilt_NONE,
    fnft_nsev_bsfilt_BASIC,
    fnft_nsev_bsfilt_FULL,
    fnft_nsev_bsfilt_MANUAL
} fnft_nsev_bsfilt_t;

/**
 * Enum that specifies how the bound states are localized. Used in
 * \link fnft_nsev_opts_t \endlink. \n \n
 * @ingroup data_types
 *  fnft_nsev_bsloc_FAST_EIGENVALUE: A root finding routine due to Aurentz et al. (see
 *  https://arxiv.org/abs/1611.02435 and https://github.com/eiscor/eiscor)
 *  with \f$ O(D^2) \f$ complexity is used to detect the roots of
 *  \f$ a(\lambda) \f$. (Note: FNFT incorporates a development version of this
 *  routine as no release was available yet.) This method is relatively slow. \n \n
 *  fnft_nsev_bsloc_NEWTON: Newton's method is used to refine a given set of initial guesses.
 *  The discretization used for the the refinement is one of the base methods \link fnft_nse_discretization_t.h \endlink.
 *  The number of iterations is specified through the field \link fnft_nsev_opts_t::niter
 *  \endlink.
 *  The array bound_states passed to \link fnft_nsev \endlink
 *  should contain the initial guesses and *K_ptr should specify the number of
 *  initial guesses. It is sufficient if bound_states and normconst_or_residues
 *  are of length *K_ptr in this case. This method can be very fast if good
 *  initial guesses for the bound states are available. The complexity is
 *  \f$ O(niter (*K\_ptr) D) \f$. \n \n
 *  fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE: This method offers a good compromise
 *  between the previous two. The method automatically finds initial guesses for
 *  the NEWTON method by first applying the FAST_EIGENVALUE method to a
 *  subsampled version of the signal. Second these initial guesses are refined
 *  using the NEWTON method. The number of samples of the subsampled signal can
 *  be controlled using the parameter Dsub in \link fnft_nsev_opts_t \endlink.
 *  If Dsub=0, the routine automatically chooses this number such that the
 *  complexity is \f$ O(D \log^2 D + niter K D) \f$, where \f$ K \f$ is the
 *  number of bound states that survived the filtering operation of the initial
 *  call to the fnft_nsev_bsloc_FAST_EIGENVALUE method w.r.t. the subsampled
 *  signal. By choosing Dsub between 2 and D, the user can be request a
 *  different number of samples. Note that algorithm uses this value only as an
 *  indication.
 */
typedef enum {
    fnft_nsev_bsloc_FAST_EIGENVALUE,
    fnft_nsev_bsloc_NEWTON,
    fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE
} fnft_nsev_bsloc_t;

/**
 * Enum that specifies the type of the discrete spectrum computed by the
 * routine. Used in \link fnft_nsev_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsev_dstype_NORMING_CONSTANTS: The array is filled with the norming constants
 *  \f$ b_k \f$. \n\n
 *  fnft_nsev_dstype_RESIDUES: The array is filled with the residues (aka spectral amplitudes)
 *  \f$ b_k\big/\frac{da(\lambda_k)}{d\lambda} \f$. \n \n
 *  fnft_nsev_dstype_BOTH: The array contains both, first the norming constants and then the
 *  residues. Note that the length of the array passed by the user has to be 2*(*K_ptr) in this case.
 */
typedef enum {
    fnft_nsev_dstype_NORMING_CONSTANTS,
    fnft_nsev_dstype_RESIDUES,
    fnft_nsev_dstype_BOTH
} fnft_nsev_dstype_t;

/**
 * Enum that specifies the type of the continuous spectrum computed by the
 * routine. Used in \link fnft_nsev_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_nsev_cstype_REFLECTION_COEFFICIENT: The array is filled with the values of
 *  \f$ b(\xi)/a(\xi) \f$ on the grid specified in the description of
 *  \link fnft_nsev \endlink. \n\n
 *  fnft_nsev_cstype_AB: The array is filled with the values of \f$a(\xi)\f$ on the grid
 *  specified in the description of \link fnft_nsev \endlink, followed by
 *  the values of \f$ b(\xi) \f$ on the same grid. Note that the length of the
 *  array contspec passed by the user has to be 2*M in this case.\n\n
 *  fnft_nsev_cstype_BOTH: The first M values of the array are filled with the
 *  values returned by the REFLECTION_COEFFICIENT method. They are followed by
 *  the 2*M values returned by the AB method. Note that the length of the array
 *  passed by the user has to be 3*M in this case.
 */
typedef enum {
    fnft_nsev_cstype_REFLECTION_COEFFICIENT,
    fnft_nsev_cstype_AB,
    fnft_nsev_cstype_BOTH
} fnft_nsev_cstype_t;

/**
 * @struct fnft_nsev_opts_t
 * @brief Stores additional options for the routine \link fnft_nsev \endlink.
 * @ingroup fnft
 * @ingroup data_types
 *
 * Use the \link fnft_nsev_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_nsev_opts_t::bound_state_filtering
 *  Controls how \link fnft_nsev \endlink decide whether a numerically found
 *  root of \f$ a(\lambda) \f$ is an actual bound state or not. \n
 *  Should be of type \link fnft_nsev_bsfilt_t \endlink.
 *
 * @var fnft_nsev_opts_t::bound_state_localization
 *  Controls how \link fnft_nsev \endlink localizes bound states. \n
 * Should be of type \link fnft_nsev_bsloc_t \endlink.
 *
 * @var fnft_nsev_opts_t::Dsub
 *   Controls how many samples are used after subsampling when bound states are
 *   localized using the fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE method. See
 *   \link fnft_nsev_bsloc_t \endlink for details.
 *
 * @var fnft_nsev_opts_t::niter
 *  For fnft_nsev_bsloc_NEWTON and fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE: Maximum
 *  number of Newton iterations to be carried out.
 * 
 *  @var fnft_nsev_opts_t::tol
 *  Some bound state localization methods such as Newton have tolerance
 *  parameters. If this value is negative, these parameters will be chosen
 *  automatically. Set to a non-negative value to chose the tolerance manually.
 *
 * @var fnft_nsev_opts_t::discspec_type
 *  Controls how \link fnft_nsev \endlink fills the array
 *  normconsts_or_residues. \n
 * Should be of type \link fnft_nsev_dstype_t \endlink.
 *
 * @var fnft_nsev_opts_t::contspec_type
 *  Controls how \link fnft_nsev \endlink fills the array
 *  contspec. \n
 * Should be of type \link fnft_nsev_cstype_t \endlink.
 *
 * @var fnft_nsev_opts_t::normalization_flag
 *  Controls whether intermediate results during the fast forward scattering
 *  step are normalized. This takes a bit longer but sometimes increases the
 *  accuracy of the results. By default, normalization is enabled (i.e., the
 *  flag is one). To disable, set the flag to zero.\n\n
 *
 * @var fnft_nsev_opts_t::discretization
 *  Controls which discretization is applied to the continuous-time Zakharov-
 *  Shabat scattering problem. See \link fnft_nse_discretization_t \endlink.\n\n
 *
 *  @var fnft_nsev_opts_t::richardson_extrapolation_flag
 *  Controls whether Richardson extrapolation is applied to try and improve
 *  the accuracy of the computed spectrum. First approximation is computed
 *  as usual using all the supplied samples. A second approximation is computed
 *  using only half the samples and it is combined with the first approximation
 *  which should ideally result in a better approximation. See Chimmalgi,
 *  Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;
 *  Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential
 *  Integrators,&quot;</a> IEEE Access 7, 2019. Note that in certain situations
 *  such as discontinuous signals, applying Richardson extrapolation may result in
 *  worse accuracy compared to the first approximation.
 *  By default, Richardson extrapolation is disabled (i.e., the
 *  flag is zero). To enable, set the flag to one.\n\n
 *
 * @var fnft_nsev_opts_t::bounding_box
 *  Array of four reals. Defines a box in the complex plane that is used for
 *  manual filtering: \n
 *  bounding_box[0] <= real(lambda) <= bounding_box[1] \n
 *  bounding_box[2] <= imag(lambda) <= bounding_box[3] \n
 */
typedef struct {
    fnft_nsev_bsfilt_t bound_state_filtering;
    fnft_nsev_bsloc_t bound_state_localization;
    FNFT_UINT niter;
    FNFT_REAL tol;
    FNFT_UINT Dsub;
    fnft_nsev_dstype_t discspec_type;
    fnft_nsev_cstype_t contspec_type;
    FNFT_INT normalization_flag;
    fnft_nse_discretization_t discretization;
    FNFT_UINT richardson_extrapolation_flag;
    FNFT_REAL bounding_box[4];
} fnft_nsev_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_nsev \endlink with
 * default settings.
 *
 * @returns A \link fnft_nsev_opts_t \endlink object with the following options.\n
 *  bound_state_filtering = fnft_nsev_bsfilt_FULL\n
 *  bound_state_localization = fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE\n
 *  niter = 100\n
 *  tol = -1.0\n
 *  discspec_type = fnft_nsev_dstype_NORMING_CONSTANTS\n
 *  contspec_type = fnft_nsev_cstype_REFLECTION_COEFFICIENT\n
 *  normalization_flag = 1\n
 *  discretization = fnft_nse_discretization_2SPLIT4B\n
 *  richardson_extrapolation_flag = 0\n
 *  bounding_box = {NAN, NAN, NAN, NAN}\n
 *
 * @ingroup fnft
 */
fnft_nsev_opts_t fnft_nsev_default_opts();

/**
 * @brief Returns the maximum number of bound states that can be detected by
 * fnft_nsev.
 *
 * @param[in] D Number of samples that will be passed to
 * \link fnft_nsev \endlink. Should be
 *  larger than zero.
 * @param[in] opts Options that will be passed to fnft_nsev. If NULL is passed,
 *  the default options will be used.
 * @return Returns the maximum number of bound states or zero on error.
 *
 * @ingroup fnft
 */
FNFT_UINT fnft_nsev_max_K(const FNFT_UINT D,
    fnft_nsev_opts_t const * const opts);

/**
 * @brief Nonlinear Fourier transform for the nonlinear Schroedinger equation
 * with vanishing boundary conditions. Fast algorithms are used if the discretization supports it.
 *
 * This routine computes the nonlinear Fourier transform for the nonlinear
 * Schroedinger equation \f[ iq_x + q_{tt} \pm 2q|q|^2=0, \quad  q=q(x,t), \f]
 * of Zakharov and Shabat (<a href="http://jetp.ac.ru/cgi-bin/e/index/e/34/1/p62?a=list">Soviet. Phys. JTEP 31(1), 1972</a>)
 * for initial conditions with vanishing boundaries
 * \f[ \lim_{t\to \pm \infty }q(x_0,t) = 0 \text{ sufficiently rapidly.} \f]
 * \n
 * The main references are:
 *      - Wahls and Poor,<a href="http://dx.doi.org/10.1109/ICASSP.2013.6638772">&quot;Introducing the fast nonlinear Fourier transform,&quot;</a> Proc. ICASSP 2013.
 *      - Wahls and Poor, <a href="http://dx.doi.org/10.1109/TIT.2015.2485944">&quot;Fast numerical nonlinear Fourier transforms,&quot;</a> IEEE Trans. Inform. Theor. 61(12), 2015.
 *      - Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot; Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot; </a>Proc. ICASSP 2018, pp. 4524-4528
 *
 * The routine also utilizes ideas from the following papers:
 *      - Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;Computation of the direct scattering transform for the nonlinear Schroedinger equation,&quot;</a> J. Comput. Phys. 102(2), 1992.
 *      - Aref, <a href="https://arxiv.org/abs/1605.06328">&quot;Control and Detection of Discrete Spectral Amplitudes in Nonlinear Fourier Spectrum,&quot;</a> Preprint, arXiv:1605.06328 [math.NA], May 2016.
 *      - Hari and Kschischang, <a href="https://doi.org/10.1109/JLT.2016.2577702">&quot;Bi-Directional Algorithm for Computing Discrete Spectral Amplitudes in the NFT,&quot; </a>J. Lightwave Technol. 34(15), 2016.
 *      - Aref et al., <a href="https://doi.org/10.1109/JLT.2018.2794475">"Modulation Over Nonlinear Fourier Spectrum: Continuous and Discrete Spectrum"</a>, J. Lightwave Technol. 36(6), 2018.
 *      - Aurentz et al., <a href="https://arxiv.org/abs/1611.02435">&quot;Roots of Polynomials: on twisted QR methods for companion matrices and pencils,&quot;</a> Preprint, arXiv:1611.02435 [math.NA]</a>, Dec. 2016.
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2932256">&quot; Soliton Phase Shift Calculation for the Korteweg–De Vries Equation,&quot;</a> IEEE Access, vol. 7, pp. 122914--122930, July 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *
 * The routine supports all discretizations of type \link fnft_nse_discretization_t \endlink. The following discretizations use fast
 * algorithms which have a computational complexity of \f$ \mathcal{O}(D\log^2 D)\f$ for \f$ D\f$ point continuous spectrum given \f$ D\f$ samples:
 *       - fnft_nse_discretization_2SPLIT1A
 *       - fnft_nse_discretization_2SPLIT1B
 *       - fnft_nse_discretization_2SPLIT2A
 *       - fnft_nse_discretization_2SPLIT2B
 *       - fnft_nse_discretization_2SPLIT2S
 *       - fnft_nse_discretization_2SPLIT2_MODAL
 *       - fnft_nse_discretization_2SPLIT3A
 *       - fnft_nse_discretization_2SPLIT3B
 *       - fnft_nse_discretization_2SPLIT3S
 *       - fnft_nse_discretization_2SPLIT4A
 *       - fnft_nse_discretization_2SPLIT4B
 *       - fnft_nse_discretization_2SPLIT5A
 *       - fnft_nse_discretization_2SPLIT5B
 *       - fnft_nse_discretization_2SPLIT6A
 *       - fnft_nse_discretization_2SPLIT6B
 *       - fnft_nse_discretization_2SPLIT7A
 *       - fnft_nse_discretization_2SPLIT7B
 *       - fnft_nse_discretization_2SPLIT8A
 *       - fnft_nse_discretization_2SPLIT8B
 *       - fnft_nse_discretization_4SPLIT4A
 *       - fnft_nse_discretization_4SPLIT4B
 *
 * The following discretizations use classical algorithms which have a computational
 * complexity of \f$ \mathcal{O}(D^2)\f$ for \f$ D\f$ point continuous spectrum given \f$ D\f$ samples:
 *       - fnft_nse_discretization_BO
 *       - fnft_nse_discretization_CF4_2
 *       - fnft_nse_discretization_CF4_3
 *       - fnft_nse_discretization_CF5_3
 *       - fnft_nse_discretization_CF6_4
 *       - fnft_nse_discretization_ES4
 *       - fnft_nse_discretization_TES4
 *
 * The accuray of the computed quantities for a given signal depends primarily on the number of samples \f$ D\f$ and the numerical method. When the exact spectrum is
 * is know, the accuracy can be quantified by defining a suitable error. The error usually decreases with increasing \f$ D\f$ assuming everthing else remains the same.
 * The rate at which the error decreases with increase in \f$ D\f$ is known as the order of the method. The orders of the various discretizations can be found at \link fnft_nse_discretization_t \endlink.
 * The orders of the discretizations which use exponential splitting schemes should be the same as their base methods but can deviate when accuracy of the splitting scheme is low.
 *
 * Application of one step Richardson extrapolation should in theory increase the order by one. However, it can deviate both ways. In case of some smooth signals the order
 * may increase by two instead of one. On the other hand for discontinuous signals it maybe deterimental to apply Richardson extrapolation. In the following case the order of the method has been observed to be less than expected:
 *       - Focusing case CF5_3, residues and norming constants have order five instead of six.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be \f$T[0]<T[1]\f$.
 * @param[in] M Number of points at which the continuous spectrum (aka
 *  reflection coefficient) should be computed.
 * @param[out] contspec Array of length M in which the routine will store the
 *  desired samples \f$ r(\xi_m) \f$ of the continuous spectrum (aka
 *  reflection coefficient) in ascending order,
 *  where \f$ \xi_m = XI[0]+m(XI[1]-XI[0])/(M-1) \f$ and \f$m=0,1,\dots,M-1\f$.
 *  Has to be preallocated by the user. If NULL is passed instead, the
 *  continuous spectrum will not be computed. By changing the options, it is
 *  also possible to compute the values of \f$ a(\xi) \f$ and \f$ b(\xi) \f$
 *  instead. In that case, twice the amount of memory has to be allocated.
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum. It should be \f$XI[0]<XI[1]\f$. Can also be
 *  NULL if contspec==NULL.
 * @param[in,out] K_ptr Upon entry, *K_ptr should contain the length of the array
 *  bound_states. Upon return, *K_ptr contains the number of actually detected
 *  bound states. If the length of the array bound_states was not sufficient
 *  to store all of the detected bound states, a warning is printed and as many
 *  bound states as possible are returned instead. Note that in order to skip
 *  the computation of the bound states completely, it is not sufficient to pass
 *  *K_ptr==0. Instead, one needs to pass bound_states==NULL.
 * @param[in,out] bound_states Array. Upon entry, only if the option
 *  `bound_state_localization=fnft_nsev_bsloc_NEWTON` is passed, this array should
 *  contain initial guesses for the bound states. Otherwise the input values are
 *  ignored. Upon return, the routine has stored the detected
 *  bound states (aka eigenvalues) in the first *K_ptr entries of this array.
 *  If NULL is passed instead, the discrete spectrum will not be computed.
 *  Has to be preallocated by the user. The user can choose an arbitrary
 *  length. Typically, D is a good choice.
 * @param[out] normconsts_or_residues Array of the same length as bound_states. Upon
 *  return, the routine has stored the residues (aka spectral amplitudes)
 *  \f$\rho_k = b_k\big/ \frac{da(\lambda_k)}{d\lambda}\f$ in the
 *  first *K_ptr entries of this array. By passing a proper opts, it is also
 *  possible to store the norming constants (the '\f$ b_k \f$') or both. Has to
 *  be pre-allocated by the user. If NULL is passed instead, the residues
 *  will not be computed.
 * @param[in] kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one.
 * @param[in] opts Pointer to a \link fnft_nsev_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_nsev_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_nsev(const FNFT_UINT D, FNFT_COMPLEX const * const q,
    FNFT_REAL const * const T, const FNFT_UINT M,
    FNFT_COMPLEX * const contspec, FNFT_REAL const * const XI,
    FNFT_UINT * const K_ptr, FNFT_COMPLEX * const bound_states,
    FNFT_COMPLEX * const normconsts_or_residues, const FNFT_INT kappa,
    fnft_nsev_opts_t * opts);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define nsev_bsfilt_NONE fnft_nsev_bsfilt_NONE
#define nsev_bsfilt_BASIC fnft_nsev_bsfilt_BASIC
#define nsev_bsfilt_FULL fnft_nsev_bsfilt_FULL
#define nsev_bsfilt_MANUAL fnft_nsev_bsfilt_MANUAL
#define nsev_bsloc_FAST_EIGENVALUE fnft_nsev_bsloc_FAST_EIGENVALUE
#define nsev_bsloc_NEWTON fnft_nsev_bsloc_NEWTON
#define nsev_bsloc_SUBSAMPLE_AND_REFINE fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE
#define nsev_dstype_NORMING_CONSTANTS fnft_nsev_dstype_NORMING_CONSTANTS
#define nsev_dstype_RESIDUES fnft_nsev_dstype_RESIDUES
#define nsev_dstype_BOTH fnft_nsev_dstype_BOTH
#define nsev_cstype_REFLECTION_COEFFICIENT fnft_nsev_cstype_REFLECTION_COEFFICIENT
#define nsev_cstype_AB fnft_nsev_cstype_AB
#define nsev_cstype_BOTH fnft_nsev_cstype_BOTH
#endif

#endif
