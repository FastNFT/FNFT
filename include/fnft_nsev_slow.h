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
* Shrinivas Chimmalgi (TU Delft) 2019.
*/

/**
 * @file fnft_nsev_slow.h
 * @brief Slow nonlinear Fourier transform for the vanishing nonlinear
 *  Schroedinger equation.
 * @ingroup fnft
 */


#include "fnft_nse_discretization_t.h"
#include "fnft_nsev.h"

FNFT_INT fnft_nsev_slow(const FNFT_UINT D, FNFT_COMPLEX * const q,
    FNFT_REAL const * const T, const FNFT_UINT M,
    FNFT_COMPLEX * const contspec, FNFT_REAL const * const XI,
    FNFT_UINT * const K_ptr, FNFT_COMPLEX * const bound_states,
    FNFT_COMPLEX * const normconsts_or_residues, const FNFT_INT kappa,
    fnft_nsev_opts_t *opts);


