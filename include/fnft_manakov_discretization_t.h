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
* Contributor:
* Lianne de Vries (TU Delft student) 2021.
* Sander Wahls (TU Delft) 2021.
*/

/**
 * @file fnft_manakov_discretization_t.h
 * @brief Lists discretizations for the Manakov Equation.
 * @ingroup fnft
 */

#ifndef FNFT__MANAKOV_DISCRETIZATION_T_H    
#define FNFT__MANAKOV_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations and methods used to compute nonlinear Fourier
 * transforms for the Manakov Equation.
 *
 * Each discretization is always related to a particular numerical method and some
 * are additionally related to a splitting-scheme.\n
 * FTES4_YZ and XSPLITYZ are methods to which an exponential splitting scheme defined by
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot; Higher Order Exponential
 * Splittings for the Fast Non-Linear Fourier Transform of the Korteweg-De Vries Equation &quot;</a>, Proc. ICASSP 2018,
 * has been applied.\n
 * For FTES4_YZ, the base method is the TES4 method from Medvedev, Vaseva, Chekhovskoy and Fedoruk
 * <a href="https://doi.org/10.1364/OE.377140">&quot; Exponential fourth order schemes for direct Zakharov-Shabat problem &quot;</a>, Optics Express 28(1) 2020.\n
 * For XSPLITYZ X denotes the order of the base method. For X=2 the second order method by Boffetta and Osborne
 * from <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;Computation of the direct scattering transform for the nonlinear Schroedinger equation &quot;</a>
 * J. Comput. Phys. 102(2), 1992, has been chosen. For X=4 the base method
 * is the fourth order
 * method CF4_2 from Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;Fast Nonlinear Fourier Transform Algorithms Using Higher Order Exponential Integrators&quot;</a>, IEEE Access (7) 2019.\n
 * For FTES4_YZ and XSPLITYZ, Y denotes the order of the chosen exponential splitting scheme from
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. `z` is type of splitting and
 * can be `A` or `B, with `A` standing for schemes implemented as defined in Prins and Wahls, 
 * <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528. `B` type of splitting are the same as `A` with the positions of the 
 * two terms in the splitting interchanged.\n
 * FTES4_suzuki is the FTES4 method where suzuki factorization (<a href="https://journals.jps.jp/doi/10.1143/JPSJ.61.3015"> J. Phys. Soc. Jpn. (1992)</a>)
 * was used to split the exponential term, as carried out by Medvedev, Chekhovskoy, Vaseva and Fedoruk in <a href="https://www.osapublishing.org/ol/abstract.cfm?uri=ol-45-7-2082">
 * &quot; Conservative multi-exponential scheme for solving the direct Zakharov–Shabat scattering problem&quot;</a>, Optics Letters 45(7) 2020.\n
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-4B`, `-6B`.
 * The choice between these is a trade-off between speed and accuracy.\n
 *
 * Fast discretizations implemented for the Manakov equation: (degree denotes the polynomial degree of a single transition matrix)\n
 * `fnft_manakov_discretization_2SPLIT3A`: Order of base method = 2, Degree = 6, Order of accuracy splitting-scheme = 3\n
 * `fnft_manakov_discretization_2SPLIT3B`: Order of base method = 2, Degree = 6, Order of accuracy splitting-scheme = 3\n
 * `fnft_manakov_discretization_2SPLIT4A`: Order of base method = 2, Degree = 8, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_2SPLIT4B`: Order of base method = 2, Degree = 4, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_2SPLIT6B`: Order of base method = 2, Degree = 12, Order of accuracy splitting-scheme = 6\n
 * `fnft_manakov_discretization_4SPLIT4A`: Order of base method = 4, Degree = 8, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_4SPLIT4B`: Order of base method = 4, Degree = 4, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_4SPLIT6B`: Order of base method = 4, Degree = 12, Order of accuracy splitting-scheme = 6\n
 * `fnft_manakov_discretization_FTES4_4A`: Order of base method = 4, Degree = 8, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_FTES4_4B`: Order of base method = 4, Degree = 8, Order of accuracy splitting-scheme = 4\n
 * `fnft_manakov_discretization_FTES4_suzuki`: Order of base method = 4, Degree = 14, Order of accuracy splitting-scheme = 4\n
 * 
 * Slow methods implemented for the Manakov equation:\n
 * `fnft_manakov_discretization_CF4_2`: base method for the 4SPLITYZ methods\n
 * `fnft_manakov_discretization_BO`: base method for the 2SPLITYZ methods\n
 *
 * Used in \link fnft_manakovv_opts_t \endlink.
 *
 * @ingroup data_types
 */

typedef enum {
    fnft_manakov_discretization_2SPLIT3A,
    fnft_manakov_discretization_2SPLIT3B,
    fnft_manakov_discretization_2SPLIT4A,
    fnft_manakov_discretization_2SPLIT4B,
	fnft_manakov_discretization_2SPLIT6B,
    fnft_manakov_discretization_4SPLIT4A,
    fnft_manakov_discretization_4SPLIT4B,
    fnft_manakov_discretization_4SPLIT6B,
    fnft_manakov_discretization_FTES4_4A,
    fnft_manakov_discretization_FTES4_4B,
	fnft_manakov_discretization_FTES4_suzuki,
    fnft_manakov_discretization_CF4_2,
    fnft_manakov_discretization_BO,
} fnft_manakov_discretization_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_discretization_t fnft_manakov_discretization_t
#define manakov_discretization_2SPLIT3A fnft_manakov_discretization_2SPLIT3A
#define manakov_discretization_2SPLIT3B fnft_manakov_discretization_2SPLIT3B
#define manakov_discretization_2SPLIT4A fnft_manakov_discretization_2SPLIT4A
#define manakov_discretization_2SPLIT4B fnft_manakov_discretization_2SPLIT4B
#define manakov_discretization_2SPLIT6B fnft_manakov_discretization_2SPLIT6B
#define manakov_discretization_4SPLIT4A fnft_manakov_discretization_4SPLIT4A
#define manakov_discretization_4SPLIT4B fnft_manakov_discretization_4SPLIT4B
#define manakov_discretization_4SPLIT6B fnft_manakov_discretization_4SPLIT6B
#define manakov_discretization_FTES4_4A fnft_manakov_discretization_FTES4_4A
#define manakov_discretization_FTES4_4B fnft_manakov_discretization_FTES4_4B
#define manakov_discretization_FTES4_suzuki fnft_manakov_discretization_FTES4_suzuki
#define manakov_discretization_CF4_2 fnft_manakov_discretization_CF4_2
#define manakov_discretization_BO fnft_manakov_discretization_BO
#endif

#endif
