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
* Sander Wahls (TU Delft) 2017.
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
* Igor Chekhovskoy 2026.
 * Irina Vaseva 2026.
*/

/**
 * @file fnft_nse_discretization_t.h
 * @brief Lists discretizations for the nonlinear Schroedinger equation.
 * @ingroup fnft
 */
#ifndef FNFT_NSE_DISCRETIZATION_T_H
#define FNFT_NSE_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations and methods used to compute nonlinear Fourier
 * transforms for the Nonlinear Schroedinger Equation.
 *
 * Each discretization is always related to a particular numerical method and some
 * are additionally related to a splitting-scheme.\n
 * `fnft_nse_discretization_BO` has been taken from Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;
 * Computation of the direct scattering transform for the nonlinear Schroedinger  equation,&quot;</a> J. Comput. Phys. 102(2), 1992.\n
 * `fnft_nse_discretization_CFx_y` are from Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;
 * Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019. 
 * They are higher-order commutator-free exponential integrators with `x` denoting the order of the method
 * and `y` the number of matrix exponentials required per signal sample. \n
 * `fnft_nse_discretization_ES4` and `fnft_nse_discretization_TES4` are fourth-order
 * methods from Medvedev, Vaseva, Chekhovskoy and  Fedoruk
 * <a href="https://doi.org/10.1364/OE.377140">&quot;
 * Exponential fourth order schemes for direct Zakharov-Shabat problem,&quot;</a> Optics Express, vol. 28, pp. 20--39, 2020.\n
 * `fnft_nse_discretization_CT4` is the conservative fourth-order method from
 * S. Medvedev, I. Vaseva, I. Chekhovskoy and M. Fedoruk,
 * <a href="https://doi.org/10.1364/OL.44.002264">&quot;Numerical algorithm with
 * fourth-order accuracy for the direct Zakharov-Shabat problem,&quot;</a> Optics
 * Letters 44(9), 2264--2267 (2019).\n
 * `fnft_nse_discretization_ES6` is the slow sixth-order exponential scheme
 * from S. Medvedev, I. Chekhovskoy, I. Vaseva and M. Fedoruk,
 * <a href="https://doi.org/10.1016/j.jcp.2021.110764">&quot;Fast sixth-order
 * algorithm based on the generalized Cayley transform for the Zakharov-Shabat
 * system associated with nonlinear Schrodinger equation,&quot;</a> J. Comput. Phys.
 * 448, 110764 (2022).\n
 * `fnft_nse_discretization_ES8` is the slow eighth-order exponential scheme
 * from S. Medvedev, I. Chekhovskoy, I. Vaseva and M. Fedoruk,
 * <a href="https://doi.org/10.48550/arXiv.2608.11892">&quot;Fast Eighth-Order
 * Padé Schemes Based on Chebyshev Polynomials for the Direct
 * Zakharov-Shabat Problem,&quot;</a> arXiv:2608.11892v1 [math.NA], preprint
 * (2026).\n
 * All above discretizations only support Newton method based bound states
 * localization (see fnft_nsev_bsloc_NEWTON of type \link fnft_nsev_bsloc_t \endlink) in \link fnft_nsev \endlink. \n 
 * The exponential spliting schemes, defined in
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528 have been applied to the second-order method by Boffetta and Osborne
 * and to the fourth-order CF4_2 method to obtain other discretizations.\n
 * The `fnft_nse_discretization_2SPLIT2_MODAL` discretization is an exception. It is the normalized Ablowitz-Ladik 
 * discretization Eq. 25 in Wahls and Vaibhav<a href="https://arxiv.org/pdf/1607.01305v2.pdf">&quot;
 * Fast Inverse Nonlinear Fourier Transforms for Continuous Spectra of Zakharov-Shabat Type
 * ,&quot;</a> Unpublished.\n 
 * All other discretizations have the notation `xSPLITyz`, where `x` is the error order 
 * of the base numerical method and `y` is the order of accuracy of splitting scheme. `z` is type of splitting and
 * can be `A`, `B` or `S`, with `A` standing for schemes implemented as defined in Prins and Wahls, 
 * <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528. `B` type of splitting are the same as `A` with the positions of the 
 * two terms in the splitting interchanged. `S` is for splittings not mentioned in above reference.\n
 * `fnft_nse_discretization_FTES4_4A` and `fnft_nse_discretization_FTES4_4B` are fast
 * versions of TES4. The TES4 correction is from the Optics Express reference above;
 * the fourth-order 4A and 4B splittings are from the Prins and Wahls ICASSP 2018 reference.
 * They support the bound-state localization methods available for fast discretizations.\n
 * `fnft_nse_discretization_FTES4_suzuki` is the conservative fast TES4 scheme based on
 * Suzuki factorization from S. Medvedev, I. Chekhovskoy, I. Vaseva and M. Fedoruk,
 * <a href="https://doi.org/10.1364/OL.387436">&quot;Conservative multi-exponential scheme
 * for solving the direct Zakharov-Shabat scattering problem,&quot;</a> Optics Letters 45(7),
 * 2082-2085 (2020).\n
 * `fnft_nse_discretization_FES4_PADE` and `fnft_nse_discretization_FES6_PADE`
 * are fast fourth- and sixth-order exponential schemes based on diagonal Padé
 * approximants from S. Medvedev, I. Chekhovskoy, I. Vaseva and M. Fedoruk,
 * <a href="https://doi.org/10.1016/j.jcp.2021.110764">&quot;Fast sixth-order
 * algorithm based on the generalized Cayley transform for the Zakharov-Shabat
 * system associated with nonlinear Schrodinger equation,&quot;</a> J. Comput. Phys.
 * 448, 110764 (2022). Their Padé degree and
 * linear-fractional-map scale are selected through the corresponding options
 * structure.\n
 * `fnft_nse_discretization_FES8_PADE` is the direct Cayley fast variant of
 * the eighth-order exponential scheme from S. Medvedev, I. Chekhovskoy,
 * I. Vaseva and M. Fedoruk,
 * <a href="https://doi.org/10.48550/arXiv.2608.11892">&quot;Fast Eighth-Order
 * Padé Schemes Based on Chebyshev Polynomials for the Direct Zakharov-Shabat
 * Problem,&quot;</a> arXiv:2608.11892v1 [math.NA], preprint (2026). The article
 * reports continuous-spectrum experiments for Padé degrees 3--6 and finds
 * the direct Cayley variants less accurate than the slow and Chebyshev-based
 * variants. In particular, evaluation of the high-degree global power-basis
 * polynomials can become ill-conditioned as the grid is refined. This direct
 * reference variant is therefore limited to the continuous spectrum; the
 * Chebyshev representation is the practical high-grid variant. Degree 7 is
 * provided by the same general Padé mechanism.\n
 * `-2S` is from G. Strang,<a href="https://link.springer.com/content/pdf/10.1007/BF00281235.pdf">&quot;
 * Accurate partial difference methods I: Linear Cauchy problems,&quot;</a> 
 * in Archive for Rational Mechanics and Analysis, 12(1), 392-402, Jan 1963. It is also
 * known as the Symmetric Weighted Sequential Splitting scheme (SWSS).\n
 *`-3S` is from Eq. 14.4 in S. Brustein and A. Mirin,<a href="https://doi.org/10.1016/0021-9991(70)90080-X">&quot;
 * Third Order Difference Methods for Hyperbolic Equations,&quot;</a> 
 * J. Comput. Phys., 5, 547-571, 1970.\n
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-2S` `-4B`, `-6B` and `-8B`.
 * The choice between these is a trade-off between speed and accuracy.
 *
 * `fnft_nse_discretization_2SPLIT1A`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 1\n
 * `fnft_nse_discretization_2SPLIT1B`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 1\n
 * `fnft_nse_discretization_2SPLIT2A`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 2\n
 * `fnft_nse_discretization_2SPLIT2B`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 2\n
 * `fnft_nse_discretization_2SPLIT2S`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 2\n
 * `fnft_nse_discretization_2SPLIT2_MODAL`: Order of base method = 2, Degree = 1, Order of accuracy of splitting-scheme = 2\n
 * `fnft_nse_discretization_2SPLIT3A`: Order of base method = 2, Degree = 3, Order of accuracy of splitting-scheme = 3\n
 * `fnft_nse_discretization_2SPLIT3B`: Order of base method = 2, Degree = 3, Order of accuracy of splitting-scheme = 3\n
 * `fnft_nse_discretization_2SPLIT3S`: Order of base method = 2, Degree = 2, Order of accuracy of splitting-scheme = 3\n
 * `fnft_nse_discretization_2SPLIT4A`: Order of base method = 2, Degree = 4, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_2SPLIT4B`: Order of base method = 2, Degree = 2, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_2SPLIT5A`: Order of base method = 2, Degree = 15, Order of accuracy of splitting-scheme = 5\n
 * `fnft_nse_discretization_2SPLIT5B`: Order of base method = 2, Degree = 15, Order of accuracy of splitting-scheme = 5\n
 * `fnft_nse_discretization_2SPLIT6A`: Order of base method = 2, Degree = 12, Order of accuracy of splitting-scheme = 6\n
 * `fnft_nse_discretization_2SPLIT6B`: Order of base method = 2, Degree = 6, Order of accuracy of splitting-scheme = 6\n
 * `fnft_nse_discretization_2SPLIT7A`: Order of base method = 2, Degree = 105, Order of accuracy of splitting-scheme = 7\n
 * `fnft_nse_discretization_2SPLIT7B`: Order of base method = 2, Degree = 105, Order of accuracy of splitting-scheme = 7\n
 * `fnft_nse_discretization_2SPLIT8A`: Order of base method = 2, Degree = 24, Order of accuracy of splitting-scheme = 8\n
 * `fnft_nse_discretization_2SPLIT8B`: Order of base method = 2, Degree = 12, Order of accuracy of splitting-scheme = 8\n
 * `fnft_nse_discretization_4SPLIT4A`: Order of base method = 4, Degree = 4, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_4SPLIT4B`: Order of base method = 4, Degree = 2, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_FTES4_4A`: Order of base method = 4, Degree = 4, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_FTES4_4B`: Order of base method = 4, Degree = 2, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_FTES4_suzuki`: Order of base method = 4, Degree = 7, Order of accuracy of splitting-scheme = 4
 * `fnft_nse_discretization_FES4_PADE`: Order of base method = 4, Padé degree = 2--7, Order of accuracy = 4\n
 * `fnft_nse_discretization_FES6_PADE`: Order of base method = 6, Padé degree = 3--7, Order of accuracy = 6
 * `fnft_nse_discretization_CT4`: Non-polynomial slow method, order of accuracy = 4\n
 * `fnft_nse_discretization_ES6`: Non-polynomial slow method, order of accuracy = 6\n
 * `fnft_nse_discretization_ES8`: Non-polynomial slow method, order of accuracy = 8\n
 * `fnft_nse_discretization_FES8_PADE`: Eighth-order base method, Padé degree = 3--7, polynomial degree = 10 times the Padé degree, order of accuracy = 6 for degree 3 and 8 for degrees 4--7\n
 *
 * Used in \link fnft_nsev_opts_t \endlink, \link fnft_nsep_opts_t \endlink
 *  and \link fnft_nsev_inverse_opts_t \endlink.
 *
 * @ingroup data_types
 */
typedef enum {
    fnft_nse_discretization_2SPLIT2_MODAL,
    fnft_nse_discretization_BO,
    fnft_nse_discretization_2SPLIT1A,
    fnft_nse_discretization_2SPLIT1B,
    fnft_nse_discretization_2SPLIT2A,
    fnft_nse_discretization_2SPLIT2B,
    fnft_nse_discretization_2SPLIT2S,
    fnft_nse_discretization_2SPLIT3A,
    fnft_nse_discretization_2SPLIT3B,
    fnft_nse_discretization_2SPLIT3S,
    fnft_nse_discretization_2SPLIT4A,
    fnft_nse_discretization_2SPLIT4B,
    fnft_nse_discretization_2SPLIT5A,
    fnft_nse_discretization_2SPLIT5B,
    fnft_nse_discretization_2SPLIT6A,
    fnft_nse_discretization_2SPLIT6B,
    fnft_nse_discretization_2SPLIT7A,
    fnft_nse_discretization_2SPLIT7B,
    fnft_nse_discretization_2SPLIT8A,
    fnft_nse_discretization_2SPLIT8B,
    fnft_nse_discretization_4SPLIT4A,
    fnft_nse_discretization_4SPLIT4B,
    fnft_nse_discretization_CF4_2,
    fnft_nse_discretization_CF4_3,
    fnft_nse_discretization_CF5_3,
    fnft_nse_discretization_CF6_4,
    fnft_nse_discretization_ES4,
    fnft_nse_discretization_TES4,
    fnft_nse_discretization_FTES4_4A,
    fnft_nse_discretization_FTES4_4B,
    fnft_nse_discretization_FTES4_suzuki,
    fnft_nse_discretization_FES4_PADE,
    fnft_nse_discretization_FES6_PADE,
    fnft_nse_discretization_CT4,
    fnft_nse_discretization_ES6,
    fnft_nse_discretization_ES8,
    fnft_nse_discretization_FES8_PADE
} fnft_nse_discretization_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_discretization_2SPLIT2_MODAL fnft_nse_discretization_2SPLIT2_MODAL
#define nse_discretization_BO fnft_nse_discretization_BO
#define nse_discretization_t fnft_nse_discretization_t
#define nse_discretization_2SPLIT1A fnft_nse_discretization_2SPLIT1A
#define nse_discretization_2SPLIT1B fnft_nse_discretization_2SPLIT1B
#define nse_discretization_2SPLIT2A fnft_nse_discretization_2SPLIT2A
#define nse_discretization_2SPLIT2B fnft_nse_discretization_2SPLIT2B
#define nse_discretization_2SPLIT2S fnft_nse_discretization_2SPLIT2S
#define nse_discretization_2SPLIT3A fnft_nse_discretization_2SPLIT3A
#define nse_discretization_2SPLIT3B fnft_nse_discretization_2SPLIT3B
#define nse_discretization_2SPLIT3S fnft_nse_discretization_2SPLIT3S
#define nse_discretization_2SPLIT4A fnft_nse_discretization_2SPLIT4A
#define nse_discretization_2SPLIT4B fnft_nse_discretization_2SPLIT4B
#define nse_discretization_2SPLIT5A fnft_nse_discretization_2SPLIT5A
#define nse_discretization_2SPLIT5B fnft_nse_discretization_2SPLIT5B
#define nse_discretization_2SPLIT6A fnft_nse_discretization_2SPLIT6A
#define nse_discretization_2SPLIT6B fnft_nse_discretization_2SPLIT6B
#define nse_discretization_2SPLIT7A fnft_nse_discretization_2SPLIT7A
#define nse_discretization_2SPLIT7B fnft_nse_discretization_2SPLIT7B
#define nse_discretization_2SPLIT8A fnft_nse_discretization_2SPLIT8A
#define nse_discretization_2SPLIT8B fnft_nse_discretization_2SPLIT8B
#define nse_discretization_4SPLIT4A fnft_nse_discretization_4SPLIT4A
#define nse_discretization_4SPLIT4B fnft_nse_discretization_4SPLIT4B
#define nse_discretization_CF4_2 fnft_nse_discretization_CF4_2
#define nse_discretization_CF4_3 fnft_nse_discretization_CF4_3
#define nse_discretization_CF5_3 fnft_nse_discretization_CF5_3
#define nse_discretization_CF6_4 fnft_nse_discretization_CF6_4
#define nse_discretization_ES4 fnft_nse_discretization_ES4
#define nse_discretization_TES4 fnft_nse_discretization_TES4
#define nse_discretization_CT4 fnft_nse_discretization_CT4
#define nse_discretization_ES6 fnft_nse_discretization_ES6
#define nse_discretization_ES8 fnft_nse_discretization_ES8
#define nse_discretization_FTES4_4A fnft_nse_discretization_FTES4_4A
#define nse_discretization_FTES4_4B fnft_nse_discretization_FTES4_4B
#define nse_discretization_FTES4_suzuki fnft_nse_discretization_FTES4_suzuki
#define nse_discretization_FES4_PADE fnft_nse_discretization_FES4_PADE
#define nse_discretization_FES6_PADE fnft_nse_discretization_FES6_PADE
#define nse_discretization_FES8_PADE fnft_nse_discretization_FES8_PADE


#endif

#endif
