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
#include "fnft_config.h"

// Default printf function, prints to stderr
INT fnft__default_printf(const char * format, ...)
{
    INT ret_code;

    va_list args;
    va_start(args, format);
    ret_code = vfprintf(stderr, format, args);             
    va_end(args);

    return ret_code;
}

// Pointer to printf functions used for error messages and warnings. Set to
// NULL to disable those. Make thread local if possible. 
static
#ifdef HAVE__THREAD_LOCAL
_Thread_local
#else
#ifdef HAVE___THREAD
__thread
#endif
#endif
fnft_printf_ptr_t fnft__printf_ptr = fnft__default_printf;

void fnft_errwarn_setprintf(fnft_printf_ptr_t printf_ptr)
{
    fnft__printf_ptr = printf_ptr;
}

fnft_printf_ptr_t fnft_errwarn_getprintf()
{
    return fnft__printf_ptr;
}

