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
 * @file fnft_version.h
 * @brief Provides a function to query the version of the FNFT library.
 * @ingroup fnft
 */

#ifndef FNFT_VERSION_H
#define FNFT_VERSION_H

#include "fnft.h"

/**
 * Provides the version of the FNFT library.
 *
 * The version is of the form $major.$minor.$patch$suffix, e.g.,
 * "0.2.1-al". The function simply returns the corresponding values
 * defined fnft_config.h.
 *
 * @param[out] major *major is overwritten with the value of FNFT_VERSION_MAJOR
 * @param[out] minor *minor is overwritten with the value of FNFT_VERSION_MINOR
 * @param[out] patch *patch is overwritten with the value of FNFT_VERSION_PATCH
 * @param[out] suffix suffix is overwritten with the value of FNFT_VERSION_SUFFIX
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 *
 * @ingroup fnft
 */
FNFT_INT fnft_version(FNFT_UINT * const major,
                      FNFT_UINT * const minor,
                      FNFT_UINT * const patch,
                      char suffix[FNFT_VERSION_SUFFIX_MAXLEN+1]);

#endif
