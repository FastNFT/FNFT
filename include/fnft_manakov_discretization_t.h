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
*/

/**
 * @file fnft_manakov_discretization_t.h
 * @brief Lists discretizations for the Manakov Equation.
 * @ingroup fnft
 */

#ifndef FNFT__MANAKOV_DISCRETIZATION_T_H    //header guard
#define FNFT__MANAKOV_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations and methods used to compute nonlinear Fourier
 * transforms for the Manakov Equation.
 *
 * Each discretization is always related to a particular numerical method and some
 * are additionally related to a splitting-scheme.\n
 * FTES4_YZ and XSPLITYZ are methods to which an exponential splitting scheme defined by
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot; has been applied.\n
 * For FTES4_YZ, the base method is the TES4 method from Medvedev, Vaseva, Chekhovskoy and Fedoruk
 * <a href="https://doi.org/10.1364/OE.377140">&quot;.\n
 * For XSPLITYZ X denotes the order of the base method. For X=2 the second order method by Boffetta and Osborne
 * from <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot; has been chosen. For X=4 the base method
 * is the fourth order
 * method CF4_2 from Chimmalgi, Prins and Wahls, <a href="https://doi.org/10.1109/ACCESS.2019.2945480">&quot;\n
 * For FTES4_YZ and XSPLITYZ, Y denotes the order of the chosen exponential splitting scheme from
 * Prins and Wahls, <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. `z` is type of splitting and
 * can be `A` or `B, with `A` standing for schemes implemented as defined in Prins and Wahls, 
 * <a href="https://doi.org/10.1109/ICASSP.2018.8461708">&quot;
 * Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot;
 * </a>Proc. ICASSP 2018, pp. 4524-4528. `B` type of splitting are the same as `A` with the positions of the 
 * two terms in the splitting interchanged.
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-4B`, `-6B` and `-8B`.
 * The choice between these is a trade-off between speed and accuracy.
 *
 * `fnft_nse_discretization_2SPLIT3A`: Order of base method = 2, Degree = 6, Order of accuracy of splitting-scheme = 3\n
 * `fnft_nse_discretization_2SPLIT3B`: Order of base method = 2, Degree = 6, Order of accuracy of splitting-scheme = 3\n
 * `fnft_nse_discretization_2SPLIT4A`: Order of base method = 2, Degree = 8, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_2SPLIT4B`: Order of base method = 2, Degree = 4, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_4SPLIT4A`: Order of base method = 4, Degree = 8, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_4SPLIT4B`: Order of base method = 4, Degree = 4, Order of accuracy of splitting-scheme = 4\n
 * `fnft_nse_discretization_4SPLIT6A`: Order of base method = 4, Degree = 12, Order of accuracy of splitting-scheme = 6\n
 * `fnft_nse_discretization_4SPLIT6B`: Order of base method = 4, Degree = 12, Order of accuracy of splitting-scheme = 6\n
 *
 * Used in \link fnft_manakov_opts_t \endlink.
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
    fnft_manakov_discretization_FTES4_4B
} fnft_manakov_discretization_t;
// This list only includes the enums for discretizations that have been implemented already. Update when new discretizations are implemented

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
#endif

#endif
