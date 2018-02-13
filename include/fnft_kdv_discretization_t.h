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
*/

/**
 * @file fnft_kdv_discretization_t.h
 * @brief Lists discretizations for the Korteweg-de Vries equation.
 * @ingroup fnft
 */
#ifndef FNFT_KDV_DISCRETIZATION_T_H
#define FNFT_KDV_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * Enum that specifies discretizations used to compute nonlinear Fourier
 * transforms for the Korteweg-de Vries equation. Used in
 * \link fnft_kdvv_opts_t \endlink.\n
 * @ingroup data_types
 */
typedef enum {
    fnft_kdv_discretization_2SPLIT1A,
    fnft_kdv_discretization_2SPLIT1B,
    fnft_kdv_discretization_2SPLIT2A,
    fnft_kdv_discretization_2SPLIT2B,
    fnft_kdv_discretization_2SPLIT3A,
    fnft_kdv_discretization_2SPLIT3B,
    fnft_kdv_discretization_2SPLIT4A,
    fnft_kdv_discretization_2SPLIT4B,
    fnft_kdv_discretization_2SPLIT5A,
    fnft_kdv_discretization_2SPLIT5B,
    fnft_kdv_discretization_2SPLIT6A,
    fnft_kdv_discretization_2SPLIT6B,
    fnft_kdv_discretization_2SPLIT7A,
    fnft_kdv_discretization_2SPLIT7B,
    fnft_kdv_discretization_2SPLIT8A,
    fnft_kdv_discretization_2SPLIT8B
} fnft_kdv_discretization_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdv_discretization_2SPLIT1A fnft_kdv_discretization_2SPLIT1A
#define kdv_discretization_2SPLIT1B fnft_kdv_discretization_2SPLIT1B
#define kdv_discretization_2SPLIT2A fnft_kdv_discretization_2SPLIT2A
#define kdv_discretization_2SPLIT2B fnft_kdv_discretization_2SPLIT2B
#define kdv_discretization_2SPLIT3A fnft_kdv_discretization_2SPLIT3A
#define kdv_discretization_2SPLIT3B fnft_kdv_discretization_2SPLIT3B
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
#define kdv_discretization_t fnft_kdv_discretization_t
#endif

#endif
