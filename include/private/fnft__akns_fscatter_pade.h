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
 * Igor Chekhovskoy (NSU, FRC ICT) 2026.
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
 * @param[in] method_order Order of the exponential scheme (4, 6 or 8).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @returns Required number of complex elements, or zero for invalid options.
 */
FNFT_UINT fnft__akns_fscatter_pade_numel(FNFT_UINT D,
        FNFT_UINT method_order, FNFT_UINT pade_degree);

/**
 * @brief Number of elements required for the scalar denominator buffer.
 *
 * @param[in] D Number of samples.
 * @param[in] method_order Order of the exponential scheme (4, 6 or 8).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @returns Required number of complex elements, or zero for invalid options.
 */
FNFT_UINT fnft__akns_fscatter_pade_den_numel(FNFT_UINT D,
        FNFT_UINT method_order, FNFT_UINT pade_degree);

/**
 * @brief Computes the rational polynomial scattering matrix for a diagonal
 * Padé approximation of a fourth-, sixth- or eighth-order exponential scheme.
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
 * @param[in] method_order Order of the exponential scheme (4, 6 or 8).
 * @param[in] pade_degree Degree of the diagonal Padé approximant.
 * @param[in] h Positive scale of the linear fractional map.
 * @param[in] periodic_flag One makes derivative stencils wrap periodically;
 *   zero makes samples outside the interval zero. Other values are invalid.
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
        const FNFT_INT periodic_flag,
        FNFT_COMPLEX * const numerator,
        FNFT_UINT * const numerator_degree,
        FNFT_INT * const numerator_exponent,
        FNFT_COMPLEX * const denominator,
        FNFT_UINT * const denominator_degree,
        FNFT_INT * const denominator_exponent);

/**
 * @brief Builds an FES8 Padé transfer matrix in the Chebyshev basis.
 *
 * The dimensionless spectral variable is mapped according to
 * \f$\zeta=c+Hx\f$. The returned numerator and denominator coefficients are
 * in ascending Chebyshev order. Their product trees require
 * \f$O(KD\log^2D)\f$ operations, where \f$K=10s\f$ and \f$s\f$ is the Padé
 * degree. periodic_flag accepts zero for zero extension and one for periodic
 * wrapping; other values are invalid.
 */
FNFT_INT fnft__akns_fscatter_pade_chebyshev(const FNFT_UINT D,
        FNFT_COMPLEX const * const q, FNFT_COMPLEX const * const r,
        const FNFT_REAL eps_t, const FNFT_UINT pade_degree,
        const FNFT_REAL c, const FNFT_REAL H,
        const FNFT_INT periodic_flag, FNFT_COMPLEX * const numerator,
        FNFT_UINT * const numerator_degree,
        FNFT_INT * const numerator_exponent,
        FNFT_COMPLEX * const denominator,
        FNFT_UINT * const denominator_degree,
        FNFT_INT * const denominator_exponent);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define akns_fscatter_pade_numel(...) fnft__akns_fscatter_pade_numel(__VA_ARGS__)
#define akns_fscatter_pade_den_numel(...) fnft__akns_fscatter_pade_den_numel(__VA_ARGS__)
#define akns_fscatter_pade(...) fnft__akns_fscatter_pade(__VA_ARGS__)
#define akns_fscatter_pade_chebyshev(...) fnft__akns_fscatter_pade_chebyshev(__VA_ARGS__)
#endif

#endif
