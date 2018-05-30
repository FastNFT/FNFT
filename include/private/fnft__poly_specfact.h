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

#ifndef FNFT__POLY_SPECFACT_H
#define FNFT__POLY_SPECFACT_H

#include "fnft.h"

FNFT_INT fnft__poly_specfact(const FNFT_UINT deg,
                             FNFT_COMPLEX const * const poly,
                             FNFT_COMPLEX * const result,
                             const FNFT_UINT oversampling_factor);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define poly_specfact(...) fnft__poly_specfact(__VA_ARGS__)
#endif

#endif
