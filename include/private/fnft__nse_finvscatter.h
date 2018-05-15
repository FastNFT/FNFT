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
 * \file fnft__nse_finvscatter.h
 * @brief Recovers the signal from a scattering matrix.
 * @ingroup nse
 */

#ifndef FNFT__NSE_FINVSCATTER_H
#define FNFT__NSE_FINVSCATTER_H

#include "fnft_nse_discretization_t.h"

FNFT_INT fnft__nse_finvscatter(
    const FNFT_UINT deg,
    FNFT_COMPLEX const * const transfer_matrix,
    FNFT_COMPLEX * const q,
    const FNFT_REAL eps_t,
    const FNFT_INT kappa,
    const fnft_nse_discretization_t discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_finvscatter(...) fnft__nse_finvscatter(__VA_ARGS__)
#endif

#endif
