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
 * @file fnft__discretization_t.h
 * @brief Lists discretizations.
 *
 * @ingroup fnft
 */
#ifndef FNFT__DISCRETIZATION_T_H
#define FNFT__DISCRETIZATION_T_H

#include "fnft.h"

/**
 * @brief Enum that specifies discretizations used to compute nonlinear Fourier
 * transforms for the Korteweg-de Vries equation.
 *
 * In general, discretizations with a lower degree are faster, while those with
 * a highter order of accuracy are more accurate. Therefore, the best choice is
 * normally among `-2A`, `-2B`, `-2S` `-4B`, `-6B` and `-8B`.
 * The choice between these is a trade-off between speed and accuracy.
 *
 * `fnft__discretization_2SPLIT1A`: Degree = 1, Order of accuracy = 1\n
 * `fnft__discretization_2SPLIT1B`: Degree = 1, Order of accuracy = 1\n
 * `fnft__discretization_2SPLIT2A`: Degree = 1, Order of accuracy = 2\n
 * `fnft__discretization_2SPLIT2B`: Degree = 1, Order of accuracy = 2\n
 * `fnft__discretization_2SPLIT2S`: Degree = 1, Order of accuracy = 2\n
 * `fnft__discretization_2SPLIT3A`: Degree = 3, Order of accuracy = 3\n
 * `fnft__discretization_2SPLIT3B`: Degree = 3, Order of accuracy = 3\n
 * `fnft__discretization_2SPLIT3S`: Degree = 2, Order of accuracy = 3\n
 * `fnft__discretization_2SPLIT4A`: Degree = 4, Order of accuracy = 4\n
 * `fnft__discretization_2SPLIT4B`: Degree = 2, Order of accuracy = 4\n
 * `fnft__discretization_2SPLIT5A`: Degree = 15, Order of accuracy = 5\n
 * `fnft__discretization_2SPLIT5B`: Degree = 15, Order of accuracy = 5\n
 * `fnft__discretization_2SPLIT6A`: Degree = 12, Order of accuracy = 6\n
 * `fnft__discretization_2SPLIT6B`: Degree = 6, Order of accuracy = 6\n
 * `fnft__discretization_2SPLIT7A`: Degree = 105, Order of accuracy = 7\n
 * `fnft__discretization_2SPLIT7B`: Degree = 105, Order of accuracy = 7\n
 * `fnft__discretization_2SPLIT8A`: Degree = 24, Order of accuracy = 8\n
 * `fnft__discretization_2SPLIT8B`: Degree = 12, Order of accuracy = 8
 *
 * These discretizations are based on exponential spliting schemes, defined in
 * Prins and Wahls, &quot;Higher order exponential splittings for the fast
 * non-linear Fourier transform of the KdV equation,&quot;
 * to appear in Proc. ICASSP 2018.
 *
 * Used in \link fnft__discretization_opts_t \endlink.
 *
 * @ingroup data_types
 */
typedef enum {
    fnft__discretization_2SPLIT1A,
    fnft__discretization_2SPLIT1B,
    fnft__discretization_2SPLIT2A,
    fnft__discretization_2SPLIT2B,
    fnft__discretization_2SPLIT2S,
    fnft__discretization_2SPLIT3A,
    fnft__discretization_2SPLIT3B,
    fnft__discretization_2SPLIT3S,
    fnft__discretization_2SPLIT4A,
    fnft__discretization_2SPLIT4B,
    fnft__discretization_2SPLIT5A,
    fnft__discretization_2SPLIT5B,
    fnft__discretization_2SPLIT6A,
    fnft__discretization_2SPLIT6B,
    fnft__discretization_2SPLIT7A,
    fnft__discretization_2SPLIT7B,
    fnft__discretization_2SPLIT8A,
    fnft__discretization_2SPLIT8B
} fnft__discretization_t;

typedef enum{
    fnft__evolution_equation_nse,
    fnft__evolution_equation_kdv
} fnft__evolution_equation_t;

typedef struct {
    fnft__discretization_t discretization;
    fnft__evolution_equation_t evolution_equation;
} fnft__discretization_opts_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define discretization_2SPLIT1A fnft__discretization_2SPLIT1A
#define discretization_2SPLIT1B fnft__discretization_2SPLIT1B
#define discretization_2SPLIT2A fnft__discretization_2SPLIT2A
#define discretization_2SPLIT2B fnft__discretization_2SPLIT2B
#define discretization_2SPLIT2S fnft__discretization_2SPLIT2S
#define discretization_2SPLIT3A fnft__discretization_2SPLIT3A
#define discretization_2SPLIT3B fnft__discretization_2SPLIT3B
#define discretization_2SPLIT3S fnft__discretization_2SPLIT3S
#define discretization_2SPLIT4A fnft__discretization_2SPLIT4A
#define discretization_2SPLIT4B fnft__discretization_2SPLIT4B
#define discretization_2SPLIT5A fnft__discretization_2SPLIT5A
#define discretization_2SPLIT5B fnft__discretization_2SPLIT5B
#define discretization_2SPLIT6A fnft__discretization_2SPLIT6A
#define discretization_2SPLIT6B fnft__discretization_2SPLIT6B
#define discretization_2SPLIT7A fnft__discretization_2SPLIT7A
#define discretization_2SPLIT7B fnft__discretization_2SPLIT7B
#define discretization_2SPLIT8A fnft__discretization_2SPLIT8A
#define discretization_2SPLIT8B fnft__discretization_2SPLIT8B
#define discretization_t fnft__discretization_t
#define discretization_opts_t fnft__discretization_opts_t
#define evolution_equation_t fnft__evolution_equation_t
#define evolution_equation_nse fnft__evolution_equation_nse
#define evolution_equation_kdv fnft__evolution_equation_kdv
#endif

#endif
