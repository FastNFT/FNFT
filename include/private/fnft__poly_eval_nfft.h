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

#ifndef FNFT__POLY_EVAL_NFFT_H
#define FNFT__POLY_EVAL_NFFT_H

#include "fnft_numtypes.h"

/**
 * Evaluate npoly power-basis polynomials at exp(i*angles[j]) using NFFT3.
 * Each coefficient array has degree+1 entries in descending order:
 * p(z) = p[0]*z^degree + ... + p[degree].
 * Each output array has M entries. Input and output arrays must not overlap.
 * All angles and coefficients must be finite. The polynomials share one plan
 * and node precomputation, with a cutoff of 12 and oversampling of at least 2.
 * For fixed npoly, the cost is O(degree*log(degree+2) + M).
 * Returns an error if NFFT3 support is not enabled or a result is not finite.
 */
FNFT_INT fnft__poly_eval_nfft_power(
        FNFT_UINT degree,
        FNFT_UINT npoly,
        FNFT_COMPLEX const * const *polynomials,
        FNFT_UINT M,
        FNFT_REAL const *angles,
        FNFT_COMPLEX * const *values);

/**
 * As above, but evaluate p(cos(angles[j])) with Chebyshev coefficients in
 * ascending order: p(x) = sum_{k=0}^degree p[k]*T_k(x).
 * Coefficients may be complex; the positive and negative Fourier modes are
 * equal, not complex conjugates.
 */
FNFT_INT fnft__poly_eval_nfft_chebyshev(
        FNFT_UINT degree,
        FNFT_UINT npoly,
        FNFT_COMPLEX const * const *polynomials,
        FNFT_UINT M,
        FNFT_REAL const *angles,
        FNFT_COMPLEX * const *values);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_eval_nfft_power(...) fnft__poly_eval_nfft_power(__VA_ARGS__)
#define poly_eval_nfft_chebyshev(...) fnft__poly_eval_nfft_chebyshev(__VA_ARGS__)
#endif

#endif
