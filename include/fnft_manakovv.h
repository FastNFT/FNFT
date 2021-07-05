#ifndef FNFT_MANAKOVV_H
#define FNFT_MANAKOVV_H

/**
 * @file fnft_manakovv.h
 * @brief Fast nonlinear Fourier transform for the vanishing Manakov equation.
 * @ingroup fnft
 */
#include "fnft__manakov_fscatter.h"
#include "fnft__manakov_scatter.h"
#include "fnft__misc.h" // for l2norm
#include <string.h> // for memcpy
#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_chirpz.h"
#include "fnft_manakov_discretization_t.h"

typedef enum{
    fnft_manakovv_bsfilt_NONE,
    fnft_manakovv_bsfilt_BASIC,
    fnft_manakovv_bsfilt_FULL
} fnft_manakovv_bsfilt_t;


typedef enum{
    fnft_manakovv_bsloc_FAST_EIGENVALUE,
    fnft_manakovv_bsloc_NEWTON,
    fnft_manakovv_bsloc_SUBSAMPLE_AND_REFINE
} fnft_manakovv_bsloc_t;

typedef enum{
    fnft_manakovv_dstype_NORMING_CONSTANTS,
    fnft_manakovv_dstype_RESIDUES,
    fnft_manakovv_dstype_BOTH
} fnft_manakovv_dstype_t;

/**
 * Enum that specifies the type of the continuous spectrum computed by the
 * routine. Used in \link fnft_manakovv_opts_t \endlink.\n \n
 * @ingroup data_types
 *  fnft_manakovv_cstype_REFLECTION_COEFFICIENT: The array is filled with the values of
 *  \f$ b1(\xi)/a(\xi) \f$ followed by the values of \f$ b2(\xi)/a(\xi) \f$ on 
 * the grid specified in the description of \link fnft_manakovv \endlink. The length of contspec should be 2*M in this case. \n\n
 *  fnft_manakovv_cstype_AB: The array is filled with the values of \f$a(\xi)\f$ on the grid
 *  specified in the description of \link fnft_manakovv \endlink, followed by
 *  the values of \f$ b1(\xi) \f$ and \f$ b2(\xi) \f$ on the same grid. Note that the length of the
 *  array contspec passed by the user has to be 3*M in this case.\n\n
 *  fnft_manakovv_cstype_BOTH: The first 2*M values of the array are filled with the
 *  values returned by the REFLECTION_COEFFICIENT method. They are followed by
 *  the 3*M values returned by the AB method. Note that the length of the array
 *  passed by the user has to be 5*M in this case.
 */
typedef enum {
	fnft_manakovv_cstype_REFLECTION_COEFFICIENT,
	fnft_manakovv_cstype_AB,
	fnft_manakovv_cstype_BOTH
} fnft_manakovv_cstype_t;

/**
 * @struct fnft_manakovv_opts_t
 * @brief Stores additional options for the routine \link fnft_manakovv \endlink.
 * @ingroup fnft
 * @ingroup data_types
 *
 * Use the \link fnft_manakovv_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_manakovv_opts_t::bound_state_filtering
 *  Controls how \link fnft_manakovv \endlink decide whether a numerically found
 *  root of \f$ a(\lambda) \f$ is an actual bound state or not. \n
 *  Should be of type \link fnft_manakovv_bsfilt_t \endlink.
 *
 * @var fnft_manakovv_opts_t::bound_state_localization
 *  Controls how \link fnft_manakovv \endlink localizes bound states. \n
 * Should be of type \link fnft_manakovv_bsloc_t \endlink.
 *
 * @var fnft_manakovv_opts_t::Dsub
 *   Controls how many samples are used after subsampling when bound states are
 *   localized using the fnft_manakovv_bsloc_SUBSAMPLE_AND_REFINE method. See
 *   \link fnft_manakovv_bsloc_t \endlink for details.
 *
 * @var fnft_manakovv_opts_t::niter
 *  Number of Newton iterations to be carried out when either the
 *  fnft_manakovv_bsloc_NEWTON or the fnft_manakovv_bsloc_SUBSAMPLE_AND_REFINE method
 *  is used.
 *
 * @var fnft_manakovv_opts_t::discspec_type
 *  Controls how \link fnft_manakovv \endlink fills the array
 *  normconsts_or_residues. \n
 * Should be of type \link fnft_manakovv_dstype_t \endlink.
 *
 * @var fnft_manakovv_opts_t::contspec_type
 *  Controls how \link fnft_manakovv \endlink fills the array
 *  contspec. \n
 * Should be of type \link fnft_manakovv_cstype_t \endlink.
 *
 * @var fnft_manakovv_opts_t::normalization_flag
 *  Controls whether intermediate results during the fast forward scattering
 *  step are normalized. This takes a bit longer but sometimes increases the
 *  accuracy of the results. By default, normalization is enabled (i.e., the
 *  flag is one). To disable, set the flag to zero.\n\n
 *
 * @var fnft_manakovv_opts_t::discretization
 *  Controls which discretization is applied to the continuous-time Zakharov-
 *  Shabat scattering problem. See \link fnft_manakov_discretization_t \endlink.\n\n
 *
 *  @var fnft_manakovv_opts_t::richardson_extrapolation_flag
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
	fnft_manakovv_bsfilt_t bound_state_filtering;
	fnft_manakovv_bsloc_t bound_state_localization;
	FNFT_UINT niter;
	FNFT_UINT Dsub;
	fnft_manakovv_dstype_t discspec_type;
	fnft_manakovv_cstype_t contspec_type;
	FNFT_INT normalization_flag;
	fnft_manakov_discretization_t discretization;
	FNFT_UINT richardson_extrapolation_flag;
} fnft_manakovv_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_manakovv \endlink with
 * default settings.
 *
 * @returns A \link fnft_manakovv_opts_t \endlink object with the following options.\n
 *  bound_state_filtering = fnft_manakovv_bsfilt_FULL\n
 *  bound_state_localization = fnft_manakovv_bsloc_SUBSAMPLE_AND_REFINE\n
 *  niter = 10\n
 *  Dsub = 0\n
 *  discspec_type = fnft_manakovv_dstype_NORMING_CONSTANTS\n
 *  contspec_type = fnft_manakovv_cstype_REFLECTION_COEFFICIENT\n
 *  normalization_flag = 1\n
 *  discretization = fnft_manakov_discretization_2SPLIT4B\n
 *  richardson_extrapolation_flag = 0\n
 *
  * @ingroup fnft
 */
fnft_manakovv_opts_t fnft_manakovv_default_opts();

/**
 * @brief Nonlinear Fourier transform for the Manakov equation
 * with vanishing boundary conditions. Fast algorithms are used if the discretization supports it.
 *
 * This routine computes the nonlinear Fourier transform for the Manakov
 * equation (vector nonlinear Schrödinger equation)\f[ iq_x + q_{tt} \pm 2q|q|^2=0, \quad  q=[q1(x,t); q2(x,t)] \f]
 * of S.V. Manakov (<a href="http://www.jetp.ac.ru/cgi-bin/dn/e_038_02_0248.pdf">Soviet Physics-JETP 38.2 (1974): 248-253</a>)
 * for initial conditions with vanishing boundaries
 * \f[ \lim_{t\to \pm \infty }q(x_0,t) = 0 \text{ sufficiently rapidly.} \f]
 * \n
 * The main references are:
 *      - Wahls and Poor,<a href="http://dx.doi.org/10.1109/ICASSP.2013.6638772">&quot;Introducing the fast nonlinear Fourier transform,&quot;</a> Proc. ICASSP 2013.
 *      - Wahls and Poor, <a href="http://dx.doi.org/10.1109/TIT.2015.2485944">&quot;Fast numerical nonlinear Fourier transforms,&quot;</a> IEEE Trans. Inform. Theor. 61(12), 2015.
 *      - Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot; Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot; </a>Proc. ICASSP 2018, pp. 4524-4528
 * 		- Chimmalgi, Prins and Wahls, <a href="https://ieeexplore-ieee-org.tudelft.idm.oclc.org/abstract/document/8856211">&quot;Fast nonlinear Fourier transform algorithms using higher order exponential integrators,&quot; <a/>IEEE Access 7 (2019): 145161-145176
 *
 * The routine also utilizes ideas from the following papers:
 *      - Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;Computation of the direct scattering transform for the nonlinear Schroedinger equation,&quot;</a> J. Comput. Phys. 102(2), 1992.
 *      - Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019.
 *      - Medvedev, Vaseva, Chekhovskoy and  Fedoruk, <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.
 *
 * The routine supports all discretizations of type \link fnft_manakov_discretization_t \endlink. The following discretizations use fast
 * algorithms which have a computational complexity of \f$ \mathcal{O}(D\log^2 D)\f$ for \f$ D\f$ point continuous spectrum given \f$ D\f$ samples:
 *       - fnft_manakov_discretization_2SPLIT3A
 *       - fnft_manakov_discretization_2SPLIT3B
 *       - fnft_manakov_discretization_2SPLIT4A
 *       - fnft_manakov_discretization_2SPLIT4B
 *       - fnft_manakov_discretization_2SPLIT6B
 *       - fnft_manakov_discretization_4SPLIT4A
 *       - fnft_manakov_discretization_4SPLIT4B
 * 	     - fnft_manakov_discretization_4SPLIT6B
 * 		 - fnft_manakov_discretization_FTES4_4A
 * 		 - fnft_manakov_discretization_FTES4_4B
 * 		 - fnft_manakov_discretization_FTES4_suzuki
 *
 * The following discretizations use classical algorithms which have a computational
 * complexity of \f$ \mathcal{O}(D^2)\f$ for \f$ D\f$ point continuous spectrum given \f$ D\f$ samples:
 *       - fnft_manakov_discretization_BO
 *       - fnft_manakov_discretization_CF4_2
 *
 * The accuray of the computed quantities for a given signal depends primarily on the number of samples \f$ D\f$ and the numerical method. When the exact spectrum is
 * is know, the accuracy can be quantified by defining a suitable error. The error usually decreases with increasing \f$ D\f$ assuming everthing else remains the same.
 * The rate at which the error decreases with increase in \f$ D\f$ is known as the order of the method. The orders of the various discretizations can be found at \link fnft_manakov_discretization_t \endlink.
 * The orders of the discretizations which use exponential splitting schemes should be the same as their base methods but can deviate when accuracy of the splitting scheme is low.
 *
 * Application of one step Richardson extrapolation should in theory increase the order by one. However, it can deviate both ways. In case of some smooth signals the order
 * may increase by two instead of one. On the other hand for discontinuous signals it maybe detrimental to apply Richardson extrapolation. In the following case the error with Richardson extrapolation has been observed to be higher than without:
 *       - Focusing case 4split6B
 *
 * @param[in] D Number of samples
 * @param[in] q1 First element of the potential function. Array of length D, contains samples \f$ q1(t_n)=q1(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q1(t_0), q1(t_1), \dots, q1(t_{D-1}) \f$)
 * * @param[in] q2 Second element of the potential function, analoguous to q1
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be \f$T[0]<T[1]\f$.
 * @param[in] M Number of points at which the continuous spectrum (aka
 *  reflection coefficient and/or NFT coefficient a, b1, b2) should be computed.
 * @param[out] contspec Array of length 2*M, 3*M or 5*M dependent on desired contspec of type \link fnft_manakovv_cstype_t \endlink
 *  in which the routine will store the reflection coefficients \f$ \rho \f$ and/or the NFT coefficients a, b1, b2
 *  in ascending order: \f$ \rho_1(\xi_m) = b1(\xi_m)/a(\xi_m) \f$
 *  followed by \f$ \rho_1(\xi_m) = b1(\xi_m)/a(\xi_m) \f$ if only the reflection coefficients are requested,
 *  where \f$ \xi_m = XI[0]+m(XI[1]-XI[0])/(M-1) \f$ and \f$m=0,1,\dots,M-1\f$.
 *  \f$ a(\xi_m) \f$, \f$ b_1(\xi_m) \f$ followed by \f$ b_2(\xi_m) \f$ if the NFT coefficients are requested and
 *  both these arrays if both are requested. If NULL is passed instead, the
 *  continuous spectrum will not be computed.
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
 *  `bound_state_localization=fnft_manakovvv_bsloc_NEWTON` is passed, this array should
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
 * @param[in] kappa =+1 for the focusing Manakov equation,
 *  =-1 for the defocusing one.
 * @param[in] opts Pointer to a \link fnft_manakovv_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_manakovv_default_opts \endlink
 *  to generate such an object and modify as desired. It is also possible to
 *  pass NULL, in which case the routine will use the default options. The
 *  user is reponsible to freeing the object after the routine has returned.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_manakovv(const FNFT_UINT D, FNFT_COMPLEX const* const q1, FNFT_COMPLEX const* const q2,
	FNFT_REAL const* const T, const FNFT_UINT M,
	FNFT_COMPLEX* const contspec, FNFT_REAL const* const XI,
	FNFT_UINT* const K_ptr, FNFT_COMPLEX* const bound_states,
	FNFT_COMPLEX* const normconsts_or_residues, const FNFT_INT kappa,
	fnft_manakovv_opts_t* opts);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakovv_bsfilt_NONE fnft_manakovv_bsfilt_NONE
#define manakovv_bsfilt_BASIC fnft_manakovv_bsfilt_BASIC
#define manakovv_bsfilt_FULL fnft_manakovv_bsfilt_FULL
#define manakovv_bsloc_FAST_EIGENVALUE fnft_manakovv_bsloc_FAST_EIGENVALUE
#define manakovv_bsloc_NEWTON fnft_manakovv_bsloc_NEWTON
#define manakovv_bsloc_SUBSAMPLE_AND_REFINE fnft_manakovv_bsloc_SUBSAMPLE_AND_REFINE
#define manakovv_dstype_NORMING_CONSTANTS fnft_manakovv_dstype_NORMING_CONSTANTS
#define manakovv_dstype_RESIDUES fnft_manakovv_dstype_RESIDUES
#define manakovv_dstype_BOTH fnft_manakovv_dstype_BOTH
#define manakovv_cstype_REFLECTION_COEFFICIENT fnft_manakovv_cstype_REFLECTION_COEFFICIENT
#define manakovv_cstype_AB fnft_manakovv_cstype_AB
#define manakovv_cstype_BOTH fnft_manakovv_cstype_BOTH
#endif

#endif      // FNFT_MANAKOVV_H