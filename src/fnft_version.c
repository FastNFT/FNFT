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
 * Sander Wahls (TU Delft) 2017-2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_version.h"
#include "fnft__errwarn.h"
#include <stdio.h>

INT fnft_version(UINT * const major,
                 UINT * const minor,
                 UINT * const patch,
                 char suffix[FNFT_VERSION_SUFFIX_MAXLEN+1])
{
    if (major == NULL || minor == NULL || patch == NULL)
        return E_INVALID_ARGUMENT(NULL pointer passed.);
    *major = FNFT_VERSION_MAJOR;
    *minor = FNFT_VERSION_MINOR;
    *patch = FNFT_VERSION_PATCH;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-zero-length"
    int rc = snprintf(suffix, FNFT_VERSION_SUFFIX_MAXLEN+1,FNFT_VERSION_SUFFIX);
#pragma GCC diagnostic pop
    if (rc < 0)
        return E_ASSERTION_FAILED; // snprintf encoding error
    if (rc > FNFT_VERSION_SUFFIX_MAXLEN)
        return E_ASSERTION_FAILED; // suffix in fnft_config.h is too long!
    return SUCCESS;
}
