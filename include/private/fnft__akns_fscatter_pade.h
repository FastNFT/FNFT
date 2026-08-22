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
 * Igor Chekhovskoy 2026.
 */

/**
 * @file fnft__akns_fscatter_pade.h
 * @brief Fast rational Padé scattering for the AKNS system.
 * @ingroup akns
 */

#ifndef FNFT__AKNS_FSCATTER_PADE_H
#define FNFT__AKNS_FSCATTER_PADE_H

#include "fnft__poly_fmult.h"

/**
 * @brief Number of elements required for the numerator matrix buffer.
 *
 * @param[in] D Number of samples.
 * @param[in] method_order Order of the exponential scheme (4 or 6).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @returns Required number of complex elements, or zero for invalid options.
 */
FNFT_UINT fnft__akns_fscatter_pade_numel(FNFT_UINT D,
        FNFT_UINT method_order, FNFT_UINT pade_degree);

/**
 * @brief Number of elements required for the scalar denominator buffer.
 *
 * @param[in] D Number of samples.
 * @param[in] method_order Order of the exponential scheme (4 or 6).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @returns Required number of complex elements, or zero for invalid options.
 */
FNFT_UINT fnft__akns_fscatter_pade_den_numel(FNFT_UINT D,
        FNFT_UINT method_order, FNFT_UINT pade_degree);

/**
 * @brief Computes the rational polynomial scattering matrix for a diagonal
 * Padé approximation of a fourth- or sixth-order exponential scheme.
 *
 * The returned approximation is numerator(w)/denominator(w), where
 * w=(ih-eps_t*lambda)/(ih+eps_t*lambda). The diagonal Padé coefficients are
 * generated from their general formula, not selected from degree-specific
 * implementations.
 *
 * @param[in] D Number of samples.
 * @param[in] q Potential samples.
 * @param[in] r Auxiliary potential samples.
 * @param[in] eps_t Sampling step size.
 * @param[in] method_order Order of the exponential scheme (4 or 6).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @param[in] h Positive scale of the linear fractional map.
 * @param[out] numerator Numerator matrix polynomial buffer.
 * @param[out] numerator_degree Degree of each numerator entry.
 * @param[out] numerator_exponent Binary normalization exponent.
 * @param[out] denominator Scalar denominator polynomial buffer.
 * @param[out] denominator_degree Degree of the denominator.
 * @param[out] denominator_exponent Binary normalization exponent.
 * @return FNFT_SUCCESS or an FNFT error code.
 */
FNFT_INT fnft__akns_fscatter_pade(const FNFT_UINT D,
        FNFT_COMPLEX const * const q, FNFT_COMPLEX const * const r,
        const FNFT_REAL eps_t, const FNFT_UINT method_order,
        const FNFT_UINT pade_degree, const FNFT_REAL h,
        FNFT_COMPLEX * const numerator,
        FNFT_UINT * const numerator_degree,
        FNFT_INT * const numerator_exponent,
        FNFT_COMPLEX * const denominator,
        FNFT_UINT * const denominator_degree,
        FNFT_INT * const denominator_exponent);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_fscatter_pade_numel(...) fnft__akns_fscatter_pade_numel(__VA_ARGS__)
#define akns_fscatter_pade_den_numel(...) fnft__akns_fscatter_pade_den_numel(__VA_ARGS__)
#define akns_fscatter_pade(...) fnft__akns_fscatter_pade(__VA_ARGS__)
#endif

#endif
