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
* Shrinivas Chimmalgi (TU Delft) 2018.
*/

/**
 * @file fnft__akns_discretization_t.h
 * @brief Lists discretizations of AKNS system.
 *
 * @ingroup akns
 */
#ifndef FNFT__AKNS_DISCRETIZATION_T_H
#define FNFT__AKNS_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations used to compute nonlinear Fourier
 * transforms of systems that fit the AKNS model.
 *
 * All discretizations have the notation `xSPLITyz`, where `x` is the error
 * order of the unsplit scheme and `y` is the order of accuracy of splitting
 * scheme. `z` is the type of splitting and can be `A`, `B`, `S`, or `_MODAL`,
 * with `A` standing for schemes implemented as defined in P.J. Prins and S.
 * Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;Higher
 * order exponential splittings for the fast non-linear Fourier transform of the
 * KdV equation,&quot;</a> in Proceedings 2018 IEEE International Conference on
 * Acoustics, Speech and Signal Processing (ICASSP) (pp. 4524-4528). Piscataway,
 * NJ, USA: IEEE. `B` type of splitting is the same as `A` with the positions of
 * the two terms in the splitting interchanged, refered to as the dual schemes
 * in the aforementioned reference. 'S' is for two splitting schemes not
 * appearing therein:\n
 * -`2S` is known as Symmetric Weighted Sequential Splitting (SWSS) and appears
 * in G. Strang,
 * <a href="https://link.springer.com/content/pdf/10.1007/BF00281235.pdf">
 * &quot;Accurate partial difference methods I: Linear Cauchy
 * problems,&quot;</a>
 * in Archive for Rational Mechanics and Analysis, 12(1), 392-402, Jan 1963;\n
 * -`3S` is obtained from Eq. 14.4 in S. Brustein and A. Mirin,
 * <a href="https://doi.org/10.1016/0021-9991(70)90080-X">&quot;Third Order
 * Difference Methods for Hyperbolic Equations,&quot;</a>
 * J. Comput. Phys., 5, 547-571, 1970.\n
 * `fnft__akns_discretization_2SPLIT2_MODAL` is the normalized Ablowitz-Ladik
 * discretization, Eq. 25 in Wahls and Vaibhav
 * <a href="https://arxiv.org/pdf/1607.01305v2.pdf">&quot;Fast Inverse Nonlinear
 * Fourier Transforms for Continuous Spectra of Zakharov-Shabat Type,&quot;</a>
 * Unpublished.
 *
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-2S` `-4B`, `-6B` and `-8B`.
 * The choice between these is a trade-off between speed and accuracy.
 *
 * `fnft__akns_discretization_2SPLIT1A`: Degree = 1, Order of accuracy = 1;\n
 * `fnft__akns_discretization_2SPLIT1B`: Degree = 1, Order of accuracy = 1;\n
 * `fnft__akns_discretization_2SPLIT2A`: Degree = 1, Order of accuracy = 2;\n
 * `fnft__akns_discretization_2SPLIT2B`: Degree = 1, Order of accuracy = 2;\n
 * `fnft__akns_discretization_2SPLIT2S`: Degree = 1, Order of accuracy = 2;\n
 * `fnft__akns_discretization_2SPLIT2_MODAL`: Degree = 1, Order of accuracy =
 * 2;\n
 * `fnft__akns_discretization_2SPLIT3A`: Degree = 3, Order of accuracy = 3;\n
 * `fnft__akns_discretization_2SPLIT3B`: Degree = 3, Order of accuracy = 3;\n
 * `fnft__akns_discretization_2SPLIT3S`: Degree = 2, Order of accuracy = 3;\n
 * `fnft__akns_discretization_2SPLIT4A`: Degree = 4, Order of accuracy = 4;\n
 * `fnft__akns_discretization_2SPLIT4B`: Degree = 2, Order of accuracy = 4;\n
 * `fnft__akns_discretization_2SPLIT5A`: Degree = 15, Order of accuracy = 5;\n
 * `fnft__akns_discretization_2SPLIT5B`: Degree = 15, Order of accuracy = 5;\n
 * `fnft__akns_discretization_2SPLIT6A`: Degree = 12, Order of accuracy = 6;\n
 * `fnft__akns_discretization_2SPLIT6B`: Degree = 6, Order of accuracy = 6;\n
 * `fnft__akns_discretization_2SPLIT7A`: Degree = 105, Order of accuracy = 7;\n
 * `fnft__akns_discretization_2SPLIT7B`: Degree = 105, Order of accuracy = 7;\n
 * `fnft__akns_discretization_2SPLIT8A`: Degree = 24, Order of accuracy = 8;\n
 * `fnft__akns_discretization_2SPLIT8B`: Degree = 12, Order of accuracy = 8.
 *
 * `fnft_nse_discretization_BO` has been taken from Boffetta and Osborne,
 * <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;Computation of
 * the direct scattering transform for the nonlinear Schroedinger
 * equation,&quot;</a> J. Comput. Phys. 102(2), 1992.
 * It is supported by \link fnft__nse_scatter.h \endlink.
 *
 * Used in \link fnft__akns_fscatter\endlink and
 * \link fnft__akns_scatter.h\endlink.
 *
 * @ingroup data_types
 */
typedef enum {
    fnft__akns_discretization_2SPLIT2_MODAL,
    fnft__akns_discretization_2SPLIT1A,
    fnft__akns_discretization_2SPLIT1B,
    fnft__akns_discretization_2SPLIT2A,
    fnft__akns_discretization_2SPLIT2B,
    fnft__akns_discretization_2SPLIT2S,
    fnft__akns_discretization_2SPLIT3A,
    fnft__akns_discretization_2SPLIT3B,
    fnft__akns_discretization_2SPLIT3S,
    fnft__akns_discretization_2SPLIT4A,
    fnft__akns_discretization_2SPLIT4B,
    fnft__akns_discretization_2SPLIT5A,
    fnft__akns_discretization_2SPLIT5B,
    fnft__akns_discretization_2SPLIT6A,
    fnft__akns_discretization_2SPLIT6B,
    fnft__akns_discretization_2SPLIT7A,
    fnft__akns_discretization_2SPLIT7B,
    fnft__akns_discretization_2SPLIT8A,
    fnft__akns_discretization_2SPLIT8B,
    fnft__akns_discretization_BO
} fnft__akns_discretization_t;


#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_discretization_2SPLIT2_MODAL fnft__akns_discretization_2SPLIT2_MODAL
#define akns_discretization_2SPLIT1A fnft__akns_discretization_2SPLIT1A
#define akns_discretization_2SPLIT1B fnft__akns_discretization_2SPLIT1B
#define akns_discretization_2SPLIT2A fnft__akns_discretization_2SPLIT2A
#define akns_discretization_2SPLIT2B fnft__akns_discretization_2SPLIT2B
#define akns_discretization_2SPLIT2S fnft__akns_discretization_2SPLIT2S
#define akns_discretization_2SPLIT3A fnft__akns_discretization_2SPLIT3A
#define akns_discretization_2SPLIT3B fnft__akns_discretization_2SPLIT3B
#define akns_discretization_2SPLIT3S fnft__akns_discretization_2SPLIT3S
#define akns_discretization_2SPLIT4A fnft__akns_discretization_2SPLIT4A
#define akns_discretization_2SPLIT4B fnft__akns_discretization_2SPLIT4B
#define akns_discretization_2SPLIT5A fnft__akns_discretization_2SPLIT5A
#define akns_discretization_2SPLIT5B fnft__akns_discretization_2SPLIT5B
#define akns_discretization_2SPLIT6A fnft__akns_discretization_2SPLIT6A
#define akns_discretization_2SPLIT6B fnft__akns_discretization_2SPLIT6B
#define akns_discretization_2SPLIT7A fnft__akns_discretization_2SPLIT7A
#define akns_discretization_2SPLIT7B fnft__akns_discretization_2SPLIT7B
#define akns_discretization_2SPLIT8A fnft__akns_discretization_2SPLIT8A
#define akns_discretization_2SPLIT8B fnft__akns_discretization_2SPLIT8B
#define akns_discretization_BO fnft__akns_discretization_BO
#define akns_discretization_t fnft__akns_discretization_t
#endif

#endif
