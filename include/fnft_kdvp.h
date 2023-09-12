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
* Sander Wahls (KIT) 2023.
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

typedef enum {
    fnft_kdvp_mstype_FLOQUET,
    fnft_kdvp_mstype_EDGEPOINTS_AND_SIGNS,
    fnft_kdvp_mstype_OPENBANDS,
    fnft_kdvp_mstype_AMPLITUDES_MODULI_FREQS
} fnft_kdvp_mstype_t;

typedef struct {
    fnft_kdvp_mstype_t mainspec_type;
    FNFT_INT normalization_flag;
    fnft_kdv_discretization_t discretization;
    FNFT_REAL grid_spacing;
} fnft_kdvp_opts_t;

fnft_kdvp_opts_t fnft_kdvp_default_opts();

FNFT_INT fnft_kdvp( const FNFT_UINT D,
                    FNFT_COMPLEX const * const q,
                    FNFT_REAL const * const T,
                    FNFT_REAL const * E,
                    FNFT_UINT * const K_ptr,
                    FNFT_REAL * const main_spec, 
                    FNFT_UINT * const M_ptr,
                    FNFT_REAL * const aux_spec,
                    FNFT_REAL * const sheet_indices,
                    fnft_kdvp_opts_t * opts_ptr);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define kdvp_mstype_FLOQUET fnft_kdvp_mstype_FLOQUET
#define kdvp_mstype_EDGEPOINTS_AND_SIGNS fnft_kdvp_mstype_EDGEPOINTS_AND_SIGNS
#define kdvp_mstype_OPENBANDS fnft_kdvp_mstype_OPENBANDS
#define kdvp_mstype_AMPLITUDES_MODULI_FREQS fnft_kdvp_mstype_AMPLITUDES_MODULI_FREQS
#endif


#endif
