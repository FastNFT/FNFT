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
* Shrinivas Chimmalgi (TU Delft) 2019.
*/

/**
 * @file fnft_kdv_discretization_t.h
 * @brief Lists discretizations for the Korteweg-de Vries equation.
 *
 * @ingroup fnft
 */
#ifndef FNFT_KDV_DISCRETIZATION_T_H
#define FNFT_KDV_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations used to compute nonlinear Fourier
 * transforms for the Korteweg-de Vries equation.
 * 
 * These discretizations are based on exponential spliting schemes, defined in
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528.\n
 * All discretizations have the notation `xSPLITyz`, where `x` is the error order 
 * of unsplit scheme and `y` is the order of accuracy of splitting scheme. `z` is type of splitting and
 * can be `A`, `B` or `S`, with `A` standing for schemes implemented as defined in Prins and Wahls, 
 * <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528. `B` type of splitting are the same as 'A' with the positions of the 
 * two terms in the splitting interchanged. 'S' is for splittings not mentioned in above reference.\n
 * `-2S` is from G. Strang,<a href="https://link.springer.com/content/pdf/10.1007/BF00281235.pdf">&quot;
 * Accurate partial difference methods I: Linear Cauchy problems,&quot;</a> 
 * in Archive for Rational Mechanics and Analysis, 12(1), 392-402, Jan 1963. It is also
 * known as the Symmetric Weighted Sequential Splitting scheme (SWSS).\n
 *`-3S` is from Eq. 14.4 in S. Brustein and A. Mirin,<a href="https://doi.org/10.1016/0021-9991(70)90080-X">&quot;
 * Third Order Difference Methods for Hyperbolic Equations,&quot;</a> 
 * J. Comput. Phys., 5, 547-571, 1970.\n
 * `fnft_kdv_discretization_BO` has been taken from Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;
 * Computation of the direct scattering transform for the nonlinear Schroedinger  equation,&quot;</a> J. Comput. Phys. 102(2), 1992. 
 * It is supported by \link fnft__kdv_scatter.h \endlink.\n 
 * `fnft_kdv_discretization_CFx_y` are from Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;
 * Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators,&quot;</a> IEEE Access 7, 2019. 
 * They are higher-order commutator-free exponential integrators with `x` denoting the order of the method
 * and `y` the number of matrix exponentials required per signal sample. 
 * They are supported by \link fnft__kdv_scatter.h \endlink.\n 
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-2S` `-4B`, `-6B` and `-8B`.
 * The choice between these is a trade-off between speed and accuracy.
 *
 * `fnft_kdv_discretization_2SPLIT1A`: Degree = 1, Order of accuracy = 1\n
 * `fnft_kdv_discretization_2SPLIT1B`: Degree = 1, Order of accuracy = 1\n
 * `fnft_kdv_discretization_2SPLIT2A`: Degree = 1, Order of accuracy = 2\n
 * `fnft_kdv_discretization_2SPLIT2B`: Degree = 1, Order of accuracy = 2\n
 * `fnft_kdv_discretization_2SPLIT2S`: Degree = 1, Order of accuracy = 2\n
 * `fnft_kdv_discretization_2SPLIT3A`: Degree = 3, Order of accuracy = 3\n
 * `fnft_kdv_discretization_2SPLIT3B`: Degree = 3, Order of accuracy = 3\n
 * `fnft_kdv_discretization_2SPLIT3S`: Degree = 2, Order of accuracy = 3\n
 * `fnft_kdv_discretization_2SPLIT4A`: Degree = 4, Order of accuracy = 4\n
 * `fnft_kdv_discretization_2SPLIT4B`: Degree = 2, Order of accuracy = 4\n
 * `fnft_kdv_discretization_2SPLIT5A`: Degree = 15, Order of accuracy = 5\n
 * `fnft_kdv_discretization_2SPLIT5B`: Degree = 15, Order of accuracy = 5\n
 * `fnft_kdv_discretization_2SPLIT6A`: Degree = 12, Order of accuracy = 6\n
 * `fnft_kdv_discretization_2SPLIT6B`: Degree = 6, Order of accuracy = 6\n
 * `fnft_kdv_discretization_2SPLIT7A`: Degree = 105, Order of accuracy = 7\n
 * `fnft_kdv_discretization_2SPLIT7B`: Degree = 105, Order of accuracy = 7\n
 * `fnft_kdv_discretization_2SPLIT8A`: Degree = 24, Order of accuracy = 8\n
 * `fnft_kdv_discretization_2SPLIT8B`: Degree = 12, Order of accuracy = 8
 *
 * Used in \link fnft_kdvv_opts_t \endlink.
 *
 * @ingroup data_types
 */
typedef enum {
    fnft_kdv_discretization_2SPLIT1A,
    fnft_kdv_discretization_2SPLIT1B,
    fnft_kdv_discretization_2SPLIT2A,
    fnft_kdv_discretization_2SPLIT2B,
    fnft_kdv_discretization_2SPLIT2S,
    fnft_kdv_discretization_2SPLIT3A,
    fnft_kdv_discretization_2SPLIT3B,
    fnft_kdv_discretization_2SPLIT3S,
    fnft_kdv_discretization_2SPLIT4A,
    fnft_kdv_discretization_2SPLIT4B,
    fnft_kdv_discretization_2SPLIT5A,
    fnft_kdv_discretization_2SPLIT5B,
    fnft_kdv_discretization_2SPLIT6A,
    fnft_kdv_discretization_2SPLIT6B,
    fnft_kdv_discretization_2SPLIT7A,
    fnft_kdv_discretization_2SPLIT7B,
    fnft_kdv_discretization_2SPLIT8A,
    fnft_kdv_discretization_2SPLIT8B,
    fnft_kdv_discretization_4SPLIT4A,
    fnft_kdv_discretization_4SPLIT4B,
    fnft_kdv_discretization_BO,
    fnft_kdv_discretization_CF4_2,
    fnft_kdv_discretization_CF4_3,
    fnft_kdv_discretization_CF5_3,
    fnft_kdv_discretization_CF6_4
} fnft_kdv_discretization_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_discretization_2SPLIT1A fnft_kdv_discretization_2SPLIT1A
#define kdv_discretization_2SPLIT1B fnft_kdv_discretization_2SPLIT1B
#define kdv_discretization_2SPLIT2A fnft_kdv_discretization_2SPLIT2A
#define kdv_discretization_2SPLIT2B fnft_kdv_discretization_2SPLIT2B
#define kdv_discretization_2SPLIT2S fnft_kdv_discretization_2SPLIT2S
#define kdv_discretization_2SPLIT3A fnft_kdv_discretization_2SPLIT3A
#define kdv_discretization_2SPLIT3B fnft_kdv_discretization_2SPLIT3B
#define kdv_discretization_2SPLIT3S fnft_kdv_discretization_2SPLIT3S
#define kdv_discretization_2SPLIT4A fnft_kdv_discretization_2SPLIT4A
#define kdv_discretization_2SPLIT4B fnft_kdv_discretization_2SPLIT4B
#define kdv_discretization_2SPLIT5A fnft_kdv_discretization_2SPLIT5A
#define kdv_discretization_2SPLIT5B fnft_kdv_discretization_2SPLIT5B
#define kdv_discretization_2SPLIT6A fnft_kdv_discretization_2SPLIT6A
#define kdv_discretization_2SPLIT6B fnft_kdv_discretization_2SPLIT6B
#define kdv_discretization_2SPLIT7A fnft_kdv_discretization_2SPLIT7A
#define kdv_discretization_2SPLIT7B fnft_kdv_discretization_2SPLIT7B
#define kdv_discretization_2SPLIT8A fnft_kdv_discretization_2SPLIT8A
#define kdv_discretization_2SPLIT8B fnft_kdv_discretization_2SPLIT8B
#define kdv_discretization_4SPLIT4A fnft_kdv_discretization_4SPLIT4A
#define kdv_discretization_4SPLIT4B fnft_kdv_discretization_4SPLIT4B
#define kdv_discretization_BO fnft_kdv_discretization_BO
#define kdv_discretization_CF4_2 fnft_kdv_discretization_CF4_2
#define kdv_discretization_CF4_3 fnft_kdv_discretization_CF4_3
#define kdv_discretization_CF5_3 fnft_kdv_discretization_CF5_3
#define kdv_discretization_CF6_4 fnft_kdv_discretization_CF6_4
#define kdv_discretization_t fnft_kdv_discretization_t
#endif

#endif
