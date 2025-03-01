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
* Sander Wahls (KIT) 2023, 2025.
*/

/**
 * @file fnft_kdvp.h
 * @brief Fast nonlinear Fourier transform for the periodic
 *  Korteweg-de Vries equation.
 * @ingroup fnft
 */

#ifndef FNFT_KDVV_H
#define FNFT_KDVV_H

#include "fnft__kdv_discretization.h"

/**
 * @struct fnft_kdvp_opts_t
 * @brief Stores additional options for the routines \link fnft_kdvp \endlink and \link fnft_kdvp_floquet \endlink.
 * @ingroup fnft
 * @ingroup data_types
 *
 * Use the \link fnft_kdvp_default_opts \endlink routine in order to generate
 * a new variable of this type with default options and modify as needed.
 *
 * @var fnft_kdvp_opts_t::normalization_flag
 *  Controls whether intermediate results during the forward scattering
 *  step are normalized. This takes a bit longer but sometimes increases the
 *  accuracy of the results. By default, normalization is enabled (i.e., the
 *  flag is one). To disable, set the flag to zero.\n\n
 *
 * @var fnft_kdvp_opts_t::discretization
 *  Controls which discretization is used to compute the monodromy matrix.
 *  See \link fnft_kdv_discretization_t \endlink. Currently,
 *  only the BO and BO_VANILLA discretizations are allowed.\n\n
 *
 * @var fnft_kdvp_opts_t::grid_spacing
 *   Grid spacing parameter for the localization of the main and auxiliary spectrum.
 *   This option currently has to be set manually. The default value returned by
 *   \link fnft_kdvp_default_opts \endlink will lead to an error. That is because
 *   currently grid_spacing is a crucial parameter that should be selected with care,
 *   but the long term plan is to add alternative methods for finding the spectra
 *   that do not require grid search.
 *
 * @var fnft_kdvp_opts_t::niter
 *  Maximum number of Newton iterations to be carried out when the initial main and auxiliary
 *  spectrum points found by grid search are refined. Can be zero or positive.
 *
 * @var fnft_kdvp_opts_t::tol
 *  Tolerance that controls when the refinement of the initial main and auxiliary
 *  spectrum points found by grid search is stopped. If <=0, a default value is used.
 */
typedef struct {
    FNFT_INT normalization_flag;
    FNFT_INT keep_degenerate_flag;
    fnft_kdv_discretization_t discretization;
    FNFT_REAL grid_spacing;
    FNFT_UINT niter;
    FNFT_REAL tol;
} fnft_kdvp_opts_t;

/**
 * @brief Creates a new options variable for \link fnft_kdvp \endlink and
 * \link fnft_kdvp_floquet \endlink with default settings.
 *
 * @returns A \link fnft_kdvp_opts_t \endlink object with the following options.\n
 *    normalization_flag = 1\n
 *    keep_degenerate_flag = 0\n
 *    discretization = kdv_discretization_BO\n
 *    grid_spacing = 0\n
 *    niter = 100\n
 *    tol = -1.0\n
 *
 * @ingroup fnft
 */
fnft_kdvp_opts_t fnft_kdvp_default_opts();

/**
 * @brief Nonlinear Fourier transform for the Korteweg-de Vries
 * equation with periodic boundary conditions.
 *
 * This routine computes the nonlinear Fourier transform for the
 * Korteweg-de Vries equation
 * \f[ q_x + 6qq_{t} + q_{ttt}=0, \quad  q=q(x,t), \f]
 * of Gardner et al. (<a href="https://doi.org/10.1103/PhysRevLett.19.1095">
 * Phys. Rev. Lett., 1967</a>)
 * for initial conditions with periodic boundaries:
 * \f[ q(x_0,t) = q(x_0,t+P), P>0. \f]
 * Currently, fast algorithms are NOT used. The complexity is O(D*L), where
 * D is the number of signal samples and L is the number of grid points needed
 * to achieve a grid spacing of at most opts_ptr->grid_spacing on a spectral
 * interval given by the user.\n
 * 
 * The definition of the NFT for the periodic KdV equation can be found in the paper
 *      - Osborne, <a href="https://doi.org/10.1016/0378-4754(94)00029-8">&quot;Automatic algorithm for the numerical inverse scattering transform of the Korteweg–de Vries equation&quot;</a> Math. Comput. Simul. 37(4-5), 1994.
 *
 * This routine however does NOT implement the automatic algorithm proposed in that paper at the moment.
 * Instead, fnft_kdvp first combines simple grid searches with regula falsi to obtain initial guesses for the
 * main and auxliary spectrum, respectively, which are then further refined using Newton's method. Since the grid_spacing
 * parameter is of crucial importance, the user must set it manually (via the opts_ptr parameter). The default
 * value will lead to an error. For the future, it is planned to use an "automatic" algorithm that does not rely
 * on grid search.
 *
 * The routine utilizes \link fnft__kdv_scatter_matrix \endlink to compute the monodromy matrix. Currently,
 * only two discretizations of the type \link fnft_kdv_discretization_t \endlink are supported:
 *       - fnft_kdv_discretization_BO(_VANILLA)
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$). The imaginary part of
 *  every sample must be 0.
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be \f$T[0]<T[1]\f$.
 * @param[in] E Array of length 2, specifies the interval of the spectral parameter
 *  in which the algorithm looks for points in the main and auxiliary spectrum. It
 *  should be \f$E[0]<E[1]\f$.
 * @param[in,out] K_ptr Initially, 2*(*K_ptr) is size of the array main_spec provided by the
 *  user. Later, the routine updates *K_ptr to the detected number of points in the main
 *  spectrum. If the length of the arrays was not sufficient to store all of the detected
 *  main spectrum points, an error is raised.
 * @param[out] main_spec Array of length 2*(*K_ptr) in which the routine will store the
 *  desired main spectrum in the form \f$ E_1, s_1, E_2, s_2, ..., E_K, s_K \f$, where
 *  the E_i are the zeros of Delta(E+-1) and Delta(E) is the Floquet discriminant.
 *  It is s_i=1 if Delta(E_i-1)=0 and s_i=-1 if Delta(E_i+1)=0.
 *  The array has to be preallocated by the user. 
 * @param[in,out] M_ptr Initially, *M_ptr is size of the arrays aux_spec and sheet_indices
 *  provided by the user (both have the same size). Later, the routine updates *M_ptr to the
 *  detected number of points in the auxiliary spectrum. If the length of the arrays was not
 *  sufficient to store all of the detected auxiliary spectrum points, an error is raised.
 * @param[out] aux_spec Array of length *M_ptr in which the routine will store the
 *  desired auxiliary spectrum in the form \f$ \mu_1, \mu_2, ..., \mu_M \f$, where
 *  the \f$\mu_i\f$ are the zeros of \f$\alpha_{21}(E)\f$, which is an element of the monodromy matrix.
 *  The array has to be preallocated by the user, or an error will occur. 
 * @param[out] sheet_indices Array of length *M_ptr in which the routine will store the
 *  sheet indices corresponding to detected auxiliary spectrum points.
 *  Has to be preallocated by the user. 
 * @param[in] opts_ptr Pointer to a \link fnft_kdvp_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_kdvp_default_opts \endlink
 *  to generate such an object and modify as desired. Currently, the grid_spacing option
 *  MUST be set by the user.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_kdvp( const FNFT_UINT D,
                    FNFT_COMPLEX * const q,
                    FNFT_REAL const * const T,
                    FNFT_REAL * const E,
                    FNFT_UINT * const K_ptr,
                    FNFT_REAL * const main_spec, 
                    FNFT_UINT * const M_ptr,
                    FNFT_REAL * const aux_spec,
                    FNFT_REAL * const sheet_indices,
                    fnft_kdvp_opts_t * opts_ptr);
/**
 * @brief Floquet discriminant and the element alpha_21 of the monodromy
 * matrix used to define the nonlinear Fourier transform for the Korteweg-de
 * Vries equation with periodic boundary conditions.
 *
 * This routine computes the Floquet discriminent and the function alpha_21(E)
 * as defined e.g. in 
 *  
 * The definition of the NFT for the periodic KdV equation can be found in the paper
 *      - Osborne, <a href="https://doi.org/10.1016/0378-4754(94)00029-8">&quot;Automatic algorithm for the numerical inverse scattering transform of the Korteweg–de Vries equation&quot;</a> Math. Comput. Simul. 37(4-5), 1994.
 *
 * Both are intermediate quantities that are normally not interesting for end
 * users it iself. One use of this routine is to choose the spectral interval
 * and grid spacing to enable the use of \link fnft_kdvp \endlink and
 * \link fnft_kdvp_ampmodfreq \endlink. The underlying numerical method is
 * described at \link fnft_kdvp \endlink.
 *
 * @param[in] D Number of samples
 * @param[in] q Array of length D, contains samples \f$ q(t_n)=q(x_0, t_n) \f$,
 *  where \f$ t_n = T[0] + n(T[1]-T[0])/(D-1) \f$ and \f$n=0,1,\dots,D-1\f$, of
 *  the to-be-transformed signal in ascending order
 *  (i.e., \f$ q(t_0), q(t_1), \dots, q(t_{D-1}) \f$). The imaginary part of
 *  every sample must be 0.
 * @param[in] T Array of length 2, contains the position in time of the first and
 *  of the last sample. It should be \f$T[0]<T[1]\f$.
 * @param[in] E Array of length 2, specifies the interval of the spectral parameter
 *  in which the algorithm looks for points in the main and auxiliary spectrum. It
 *  should be \f$E[0]<E[1]\f$.
 * @param[in] L Number of points in the interval [E[0], E[1]] at which the Floquet
 *  discriminant and alpha_21(E) will be evaluated.
 * @param[out] DEL Array of length L in which the routine will store the
 *  desired values of the Floquet discriminant at the locations \f$E_i = E[0]+i dE\f$, where
 *  \f$ dE = (E[1]-E[0])/(L-1) \f$ and \f$i=0,1,2,\dots, L-1\f$.
 *  The array has to be preallocated by the user. 
 * @param[out] al21 Array of length L in which the routine will store the
 *  desired values of al21 at the locations \f$E_i = E[0]+i dE\f$, where
 *  \f$ dE = (E[1]-E[0])/(L-1) \f$ and \f$i=0,1,2,\dots, L-1\f$.
 *  The array has to be preallocated by the user. 
 * @param[in] opts_ptr Pointer to a \link fnft_kdvp_opts_t \endlink object. The object
 *  can be used to modify the behavior of the routine. Use
 *  the routine \link fnft_kdvp_default_opts \endlink
 *  to generate such an object and modify as desired. Currently, the grid_spacing option
 *  MUST be set by the user.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_kdvp_floquet( const FNFT_UINT D,
                            FNFT_COMPLEX * const q,
                            FNFT_REAL const * const T,
                            FNFT_REAL * const E,
                            const FNFT_UINT L,
                            FNFT_REAL * const DEL, 
                            FNFT_REAL * const al21,
                            fnft_kdvp_opts_t * opts_ptr);
/**
 * @brief Convert main spectra into amplitudes, moduli and frequencies.
 *
 * This routine takes a main spectrum computed by \link fnft_kdvp \endlink
 * and determines the amplitudes, moduli and frequencies of the corresponding
 * hyperelliptic modes. The idea goes back to Osborne and coworkers.
 * See, e.g., 
 *      - Osborne and Bergamasco, <a href="https://doi.org/10.1016/0167-2789(86)90160-0">&quot;The solitons of Zabusky and Kruskal revisited: Perspective in terms of the periodic spectral transform&quot;</a>, Physica D 18(1-3), 1986.
 * Various slightly different definitions can be found in the literature. This
 * algorithm uses the ones described in
 *      - Brühl et al., <a href="https://doi.org/10.1016/j.wavemoti.2022.102905">&quot;Comparative analysis of bore propagation over long distances using conventional linear and KdV-based nonlinear Fourier transform&quot;</a>, Wave Motion 111, 2022.
 * with the small exception that the frequencies are scaled to be in Hz, and
 * the ampltiudes of the radiation modes are scaled by two (like the single-sided FFT).
 *
 * @param[in] K_ptr
 * @param[in] K_ptr Initially, 2*(*K_ptr) is size of the array main_spec provided by the
 *  user. Later, the routine updates *K_ptr to the detected number of hyperelliptic modes.
 * @param[out] main_spec Main spectrum computed by \link fnft_kdvp \endlink
 * @param[out] ampmodspec Contains the amplitudes A_i, moduli m_i and frequencies f_i
 *  of the hyperelliptic modes, where i=0,1,...,Kout-1 with Kout being the value of *K_ptr
 *  upon exit. The ordering of the data is A_1, m_1, f_1, ..., A_Kout, m_Kout, f_Kout.
 *  Needs to be preallocated by the user and to be of size 3*Kin, where Kin is the value
 *  *K_ptr upon entry. 
 */
FNFT_INT fnft_kdvp_ampmodfreq( FNFT_UINT * const K_ptr,
                               FNFT_REAL const * const main_spec,
                               FNFT_REAL * const ampmodfreq);

#endif
