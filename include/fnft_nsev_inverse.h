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

#ifndef FNFT_NSEV_INVERSE_H
#define FNFT_NSEV_INVERSE_H

#include "fnft.h"
#include "fnft_nse_discretization_t.h"

typedef enum {
    fnft_nsev_inverse_contspec_inversion_method_REFL_COEFF,
    fnft_nsev_inverse_contspec_inversion_method_B_FROM_A,
    fnft_nsev_inverse_contspec_inversion_method_B_FROM_A_WO_SPECFACT
} fnft_nsev_inverse_contspec_inversion_method_t;

typedef struct {
    fnft_nse_discretization_t discretization;
    fnft_nsev_inverse_contspec_inversion_method_t contspec_inversion_method;
} fnft_nsev_inverse_opts_t;

fnft_nsev_inverse_opts_t fnft_nsev_inverse_default_opts();

FNFT_INT fnft_nsev_inverse_XI(
    const FNFT_UINT D,
    FNFT_REAL const * const T,
    const FNFT_UINT M,
    FNFT_REAL * const XI,
    const fnft_nse_discretization_t discretization);

FNFT_INT fnft_nsev_inverse(
    const FNFT_UINT M,
    FNFT_COMPLEX * const contspec,
    FNFT_REAL const * const XI,
    FNFT_UINT const K,
    FNFT_COMPLEX const * const bound_states,
    FNFT_COMPLEX const * const normconsts_or_residues,
    const FNFT_UINT D,
    FNFT_COMPLEX * const q,
    FNFT_REAL const * const T,
    const FNFT_INT kappa,
    fnft_nsev_inverse_opts_t *opts_ptr);

#endif
