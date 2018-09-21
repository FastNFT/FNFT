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

#include <stdio.h>
#include <stdarg.h>
#include "fnft_errwarn.h"
#include "fnft__errwarn.h"
#include "fnft_config.h"

INT fnft__errmsg_aux(const INT ec, const char *func, const INT line,
    const char *msg)
{
    fnft_printf_ptr_t printf_ptr = fnft_errwarn_getprintf();
    if (printf_ptr != NULL)
        printf_ptr("FNFT Error: %s\n in %s(%i)-%d.%d.%d%s\n", msg, func, line,
                   FNFT_VERSION_MAJOR, FNFT_VERSION_MINOR, FNFT_VERSION_PATCH,
                   FNFT_VERSION_SUFFIX);
    return ec;
}

void fnft__warn_aux(const char *func, const INT line, const char *msg)
{
    fnft_printf_ptr_t printf_ptr = fnft_errwarn_getprintf();
    if (printf_ptr != NULL)
        printf_ptr("FNFT Warning: %s\n in %s(%i)-%d.%d.%d%s\n", msg, func, line,
                   FNFT_VERSION_MAJOR, FNFT_VERSION_MINOR, FNFT_VERSION_PATCH,
                   FNFT_VERSION_SUFFIX);
}
