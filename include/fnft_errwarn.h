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

/**
 * @file fnft_errwarn.h
 * @ingroup errwarn
 *
 * Provides error codes and functions that control error messages output.
 */

#ifndef FNFT_ERRWARN_H
#define FNFT_ERRWARN_H

#include "fnft_numtypes.h"

/**
 * Function pointer to a printf-like function.
 * @ingroup errwarn
 */
typedef FNFT_INT (* fnft_printf_ptr_t) (const char *, ...);

/**
 * Value that routines return in the absence of errors.
 * @ingroup errwarn
 */
#define FNFT_SUCCESS 0

/* --- Use the following macros when determining the type of an error. --- */

/**
 * Error code for out of memory.
 * @ingroup errwarn
 */
#define FNFT_EC_NOMEM               1

/**
 * Error code for invalid arguments.
 * @ingroup errwarn
 */
#define FNFT_EC_INVALID_ARGUMENT    2

/**
 * Error code for division by zero.
 * @ingroup errwarn
 */
#define FNFT_EC_DIV_BY_ZERO         3

/**
 * Error code for test failure.
 * @ingroup errwarn
 */
#define FNFT_EC_TEST_FAILED         4

/**
 * Error code for other, unspecified errors.
 * @ingroup errwarn
 */
#define FNFT_EC_OTHER               5

/**
 * Error code for a feature that has not yet been implemented.
 * @ingroup errwarn
 */
#define FNFT_EC_NOT_YET_IMPLEMENTED 6

/**
 * Error code if a sanity check of a final or intermediate result failed.
 * @ingroup errwarn
 */
#define FNFT_EC_SANITY_CHECK_FAILED 7

/**
 * Error code if an assertion failed.
 * @ingroup errwarn
 */
#define FNFT_EC_ASSERTION_FAILED 8

/**
 * Sets the printf function that FNFT uses to print errors
 * and warnings.
 * @ingroup errwarn
 */
void fnft_errwarn_setprintf(fnft_printf_ptr_t printf_ptr);

/**
 * Returns a pointer to the printf function that FNFT uses to print
 * errors and warnings.
 * @ingroup errwarn
 */
fnft_printf_ptr_t fnft_errwarn_getprintf();

#ifdef FNFT_ENABLE_SHORT_NAMES
#define SUCCESS             FNFT_SUCCESS
#endif

#endif
