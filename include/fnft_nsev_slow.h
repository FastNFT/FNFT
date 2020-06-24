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
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
*/

/**
 * @file fnft_nsev_slow.h
 * @brief Slow nonlinear Fourier transform for the vanishing nonlinear
 *  Schroedinger equation.
 * @ingroup fnft
 */

#ifndef FNFT_NSEV_SLOW_H
#define FNFT_NSEV_SLOW_H

#include "fnft_nsev.h"

/**
 * @struct fnft_nsev_slow_opts_t
 * @brief Stores additional options for the routine \link fnft_nsev_slow \endlink.
 * @ingroup fnft
 * @ingroup data_types
 *
 * Use the \link fnft_nsev_slow_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_nsev_slow_opts_t::bound_state_filtering
 *  Controls how \link fnft_nsev_slow \endlink decide whether a numerically found
 *  root of \f$ a(\lambda) \f$ is an actual bound state or not. \n
 *  Should be of type \link fnft_nsev_bsfilt_t \endlink.
 *
 * @var fnft_nsev_slow_opts_t::bound_state_localization
 *  Controls how \link fnft_nsev_slow \endlink localizes bound states. \n
 *  Currently only fnft_nsev_bsloc_NEWTON of type \link fnft_nsev_bsloc_t \endlink
 *  is supported.
 *
 * @var fnft_nsev_slow_opts_t::niter
 *  Number of Newton iterations to be carried out when the
 *  fnft_nsev_bsloc_NEWTON method is used.
 *
 * @var fnft_nsev_slow_opts_t::discspec_type
 *  Controls how \link fnft_nsev \endlink fills the array
 *  normconsts_or_residues. \n
 *  Should be of type \link fnft_nsev_dstype_t \endlink.
 *
 * @var fnft_nsev_slow_opts_t::contspec_type
 *  Controls how \link fnft_nsev \endlink fills the array
 *  contspec. \n
 *  Should be of type \link fnft_nsev_cstype_t \endlink.
 *
 * @var fnft_nsev_slow_opts_t::discretization
 *  Controls which discretization is applied to the continuous-time Zakharov-
 *  Shabat scattering problem. See \link fnft_nse_discretization_t \endlink.
 *
 * @var fnft_nsev_slow_opts_t::richardson_extrapolation_flag
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
 *  flag is zero). To enable, set the flag to one.
 */
typedef struct {
    fnft_nsev_bsfilt_t bound_state_filtering;
    fnft_nsev_bsloc_t bound_state_localization;
    FNFT_UINT niter;
    fnft_nsev_dstype_t discspec_type;
    fnft_nsev_cstype_t contspec_type;
    fnft_nse_discretization_t discretization;
    FNFT_UINT richardson_extrapolation_flag;
} fnft_nsev_slow_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_nsev_slow \endlink with
 * default settings.
 *
 * @returns A \link fnft_nsev_slow_opts_t \endlink object with the following options.\n
 *  bound_state_filtering = fnft_nsev_bsfilt_FULL\n
 *  bound_state_localization = fnft_nsev_bsloc_NEWTON\n
 *  niter = 10\n
 *  discspec_type = fnft_nsev_dstype_NORMING_CONSTANTS\n
 *  contspec_type = fnft_nsev_cstype_REFLECTION_COEFFICIENT\n
 *  discretization = nse_discretization_BO\n
 *  richardson_extrapolation_flag = 0
 *
  * @ingroup fnft
 */
fnft_nsev_slow_opts_t fnft_nsev_slow_default_opts();

/**
 * @brief Slow nonlinear Fourier transform for the nonlinear Schroedinger
 *  equation with vanishing boundary conditions.
 *
 * This routine computes the nonlinear Fourier transform for the nonlinear
 * Schroedinger equation \f[ iq_x + q_{tt} \pm 2q|q|^2=0, \quad  q=q(x,t), \f]
 * of Zakharov and Shabat (<a href="http://jetp.ac.ru/cgi-bin/e/index/e/34/1/p62?a=list">Soviet. Phys. JTEP 31(1), 1972</a>)
 * for initial conditions with vanishing boundaries
 * \f[ \lim_{t\to \pm \infty }q(x_0,t) = 0 \text{ sufficiently rapidly.} \f]
 * \n
 * The main references are:
 *      - Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot; Computation of the direct scattering transform for the nonlinear Schroedinger  equation,&quot;</a> J. Comput. Phys. 102(2), 1992.
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot; Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *      - Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2932256">&quot; Soliton Phase Shift Calculation for the Korteweg–De Vries Equation,&quot;</a> IEEE Access, vol. 7, pp. 122914--122930, July 2019.
 *
 * The accuray of the computed quantities for a given signal depends primarily on the number of samples \f$ D\f$ and the numerical method. When the exact spectrum is 
 * is know, the accuracy can be quantified by defining a suitable error. The error usually decreases with increasing \f$ D\f$ assuming everthing else remains the same. 
 * The rate at which the error decreases with increase in \f$ D\f$ is known as the order of the method. The orders of the various discretizations can be found at \link fnft_nse_discretization_t \endlink.
 * In the following cases the orders of the methods has been observed to be less than expected:
 *       - Focusing case CF4_3, residues have order two instead of four.
 *       - Focusing case CF6_4, residues have order three instead of six.
 *       - Focusing case TES4, residues have order two instead of four.
 *
 * Application of one step Richardson extrapolation should in theory increase the order by one. However, it can deviate both ways. In case of some smooth signals the order
 * may increase by two instead of one. On the other hand for discontinuous signals it maybe deterimental to apply Richardson extrapolation. In the following cases the orders of the methods has been observed to be less than expected:
 *       - Focusing case CF4_3, residues have order two instead of five.
 *       - Focusing case CF5_3, residues have order five instead of six.
 *       - Focusing case CF6_4, residues have order four instead of seven.
 *       - Focusing case TES4, residues have order two instead of five.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$)
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be \f$T[0]<T[1]\f$).
 * @param[in] M Number of points at which the continuous spectrum (aka
 *  reflection coefficient) should be computed.
 * @param[out] contspec Array of length M in which the routine will store the
 *  desired samples \f$ r(\xi_m) \f$ of the continuous spectrum (aka
 *  reflection coefficient) in ascending order,
 *  where \f$ \xi_m = XI[0]+(XI[1]-XI[0])/(M-1) \f$ and \f$m=0,1,\dots,M-1\f$.
 *  Has to be preallocated by the user. If NULL is passed instead, the
 *  continuous spectrum will not be computed. By changing the options, it is
 *  also possible to compute the values of \f$ a(\xi) \f$ and \f$ b(\xi) \f$
 *  instead. In that case, twice the amount of memory has to be allocated.
 * @param[in] XI Array of length 2, contains the position of the first and the last
 *  sample of the continuous spectrum. It should be \f$XI[0]<XI[1]\f$. Can also be
 *  NULL if contspec==NULL.
 * @param[in,out] K_ptr Upon entry, *K_ptr should contain the length of the array
 *  bound_states. Upon return, *K_ptr contains the number of actually detected
 *  bound states. Note that in order to skip
 *  the computation of the bound states completely, it is not sufficient to pass
 *  *K_ptr==0. Instead, one needs to pass bound_states==NULL.
 * @param[out] bound_states Array. Upon return, the routine has stored the detected
 *  bound states (aka eigenvalues) in the first *K_ptr entries of this array.
 *  If NULL is passed instead, the discrete spectrum will not be computed.
 *  Has to be preallocated by the user. The user can choose an arbitrary
 *  length. Currently only NEWTON method is supported hence bound_states
 *  should contain initial guesses of the bound states.
 * @param[out] normconsts_or_residues Array of the same length as bound_states. Upon
 *  return, the routine has stored the residues (aka spectral amplitudes)
 *  \f$\rho_k = b_k\big/ \frac{da(\lambda_k)}{d\lambda}\f$ in the
 *  first *K_ptr entries of this array. By passing a proper opts, it is also
 *  possible to store the norming constants (the '\f$ b_k \f$') or both. Has to
 *  be pre-allocated by the user. If NULL is passed instead, the residues
 *  will not be computed.
 * @param[in] kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one
 * @param[in] opts Pointer to a \link fnft_nsev_slow_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_nsev_slow_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_nsev_slow(const FNFT_UINT D, FNFT_COMPLEX * const q,
    FNFT_REAL const * const T, const FNFT_UINT M,
    FNFT_COMPLEX * const contspec, FNFT_REAL const * const XI,
    FNFT_UINT * const K_ptr, FNFT_COMPLEX * const bound_states,
    FNFT_COMPLEX * const normconsts_or_residues, const FNFT_INT kappa,
    fnft_nsev_slow_opts_t *opts);


#endif
