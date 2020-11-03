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

int main()
{
    UINT major, minor, patch;
    char suffix[FNFT_VERSION_SUFFIX_MAXLEN+1];

    INT ret_code = fnft_version(&major, &minor, &patch, suffix);
    if (ret_code != SUCCESS) {
        ret_code = E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    printf("%u.%u.%u%s\n",
           (unsigned int)major,
           (unsigned int)minor,
           (unsigned int)patch,
           suffix);
    return EXIT_SUCCESS;
}
