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
 * @file fnft__poly_specfact.h
 * @brief Spectral factorization of polynomials.
 * @ingroup poly
 */

#ifndef FNFT__POLY_SPECFACT_H
#define FNFT__POLY_SPECFACT_H

#include "fnft.h"

/**
 * @brief Spectral factorization of polynomial.
 *
 * @ingroup poly
 * Computes the spectral factor of given polynomial with no zeros inside the unit circle.
 * @param [in] deg Degree of the polynomial.
 * @param [in] poly Array of coefficients for the polynomial.
 * @param [out] result Array in which the coefficients of the result are
 *   stored.
 * @param [in] oversampling_factor Oversampling factor is used to improve conditioning
 * of the spectral factorization problem.
 * @param [in] kappa kappa =+1 for the focusing nonlinear Schroedinger equation,
 *  =-1 for the defocusing one.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 */
FNFT_INT fnft__poly_specfact(const FNFT_UINT deg,
                             FNFT_COMPLEX const * const poly,
                             FNFT_COMPLEX * const result,
                             const FNFT_UINT oversampling_factor,
                             const FNFT_INT kappa);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_specfact(...) fnft__poly_specfact(__VA_ARGS__)
#endif

#endif
