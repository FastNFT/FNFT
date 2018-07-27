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
 * @file fnft__errwarn.h
 * @ingroup private_errwarn
 *
 * Provides macros for raising and checking errors.
 */

#ifndef FNFT__ERRWARN_H
#define FNFT__ERRWARN_H

#include "fnft_errwarn.h"

/**
 * Provides a warning message to the user.
 * @ingroup private_errwarn
 */
#define FNFT__WARN(msg)              fnft__warn_aux(__func__, __LINE__, msg);

/**
 * Macro for raising an out of memory error.
 * @ingroup private_errwarn
 */
#define FNFT__E_NOMEM                FNFT__ERRMSG(FNFT_EC_NOMEM, "Out of memory.")

/**
 * Macro for raising an invalid argument error.
 * @ingroup private_errwarn
 */
#define FNFT__E_INVALID_ARGUMENT(name) FNFT__ERRMSG(FNFT_EC_INVALID_ARGUMENT, FNFT__E_INVALID_ARGUMENT_(name))

/**
 * Macro for raising an subroutine error, where ec is the error code returned
 * by the subroutine. The error code returned by this macro is -abs(ec). Since
 * all other error codes are positive, a negative error code therefore
 * indicates that a subroutine failed with the error code abs(ec).
 * @ingroup private_errwarn
 */
#define FNFT__E_SUBROUTINE(ec)       FNFT__ERRMSG(-abs((FNFT_INT) ec), "Subroutine failure.")

/**
 * Macro for raising an division by zero error.
 * @ingroup private_errwarn
 */
#define FNFT__E_DIV_BY_ZERO          FNFT__ERRMSG(FNFT_EC_DIV_BY_ZERO, "Division by zero.")

/**
 * Macro for raising a test failed error.
 * @ingroup private_errwarn
 */
#define FNFT__E_TEST_FAILED          FNFT__ERRMSG(FNFT_EC_TEST_FAILED, "Test failed.")

/**
 * Macro for raising an unspecified error, where msg is a custom error message.
 * @ingroup private_errwarn
 */
#define FNFT__E_OTHER(msg)           FNFT__ERRMSG(FNFT_EC_OTHER, msg)

/**
 * Macro for raising a not yet implemented error, where name is an identifier
 * for the feature and msg is an additional custom error message.
 * @ingroup private_errwarn
 */
#define FNFT__E_NOT_YET_IMPLEMENTED(name,msg)    FNFT__ERRMSG(FNFT_EC_NOT_YET_IMPLEMENTED, FNFT__E_NOT_YET_IMPLEMENTED_(name,msg))

/**
 * Macro for raising a sanity check failed error, where msg is a custom error message.
 * @ingroup private_errwarn
 */
#define FNFT__E_SANITY_CHECK_FAILED(msg) FNFT__ERRMSG(FNFT_EC_SANITY_CHECK_FAILED, FNFT__E_SANITY_CHECK_FAILED_(msg))

/**
 * Macro for raising an assertion check failed error.
 * @ingroup private_errwarn
 */
#define FNFT__E_ASSERTION_FAILED FNFT__ERRMSG(FNFT_EC_ASSERTION_FAILED, "Assertion failed.")

/**
 * Macro for checking if a return code indicates an error. If ret_code!=FNFT_SUCCESS, a
 * subroutine error is risen and the macro uses goto to jump to label.
 * @ingroup private_errwarn
 */
#define FNFT__CHECK_RETCODE(ret_code, label) {if (ret_code!=FNFT_SUCCESS) {ret_code=FNFT__E_SUBROUTINE(ret_code);goto label;}}

#ifdef FNFT_ENABLE_SHORT_NAMES
#define WARN(msg)           FNFT__WARN(msg)
#define E_NOMEM             FNFT__E_NOMEM
#define E_INVALID_ARGUMENT(name) FNFT__E_INVALID_ARGUMENT(name)
#define E_SUBROUTINE(ec)    FNFT__E_SUBROUTINE(ec)
#define E_DIV_BY_ZERO       FNFT__E_DIV_BY_ZERO
#define E_TEST_FAILED       FNFT__E_TEST_FAILED
#define E_OTHER(msg)        FNFT__E_OTHER(msg)
#define E_NOT_YET_IMPLEMENTED(name,msg) FNFT__E_NOT_YET_IMPLEMENTED(name,msg)
#define E_SANITY_CHECK_FAILED FNFT__E_SANITY_CHECK_FAILED(msg)
#define E_ASSERTION_FAILED  FNFT__E_ASSERTION_FAILED
#define CHECK_RETCODE(ret_code,label) FNFT__CHECK_RETCODE(ret_code,label)
#endif

/* --- Auxiliary macros and functions. Do not use directly. --- */

/**
 * Macro that prints the given error message using \link fnft__errmsg_aux 
 * \endlink and returns the error code ec. Auxiliary macro. Do not call
 * directly.
 * @ingroup private_errwarn
 */
#define FNFT__ERRMSG(ec,msg)         fnft__errmsg_aux(ec, __func__, __LINE__, msg);

/**
 * Auxiliary macro that is used to stringify the input to
 * \link FNFT__E_INVALID_ARGUMENT \endlink. Do not call directly.
 * @ingroup private_errwarn
 */
#define FNFT__E_INVALID_ARGUMENT_(name) "Invalid argument "#name"."

/**
 * Auxiliary macro used used to to stringify the input to
 * \link FNFT__E_NOT_YET_IMPLEMENTED \endlink. Do not call directly.
 * @ingroup private_errwarn
 */
#define FNFT__E_NOT_YET_IMPLEMENTED_(name,msg) "Not yet implemented ("#name"). "#msg

/**
 * Auxiliary macro used to stringify the input to
 * \link FNFT__E_NOT_YET_IMPLEMENTED \endlink. Do not call directly.
 * @ingroup private_errwarn
 */
#define FNFT__E_SANITY_CHECK_FAILED_(msg) "Sanity check failed ("#msg")."

/**
 * Auxiliary function that prints a formated error message using
 * the printf function returned by \link fnft_errwarn_getprintf \endlink.
 * It simply returns the error code ec. Do not call directly.
 * @ingroup private_errwarn
 */
FNFT_INT fnft__errmsg_aux(const FNFT_INT ec, const char *func,
    const FNFT_INT line, const char *msg);

/**
 * Auxiliary function that prints a formated warning message using
 * the printf function returned by \link fnft_errwarn_getprintf \endlink.
 * Do not call directly.
 * @ingroup private_errwarn
 */
void fnft__warn_aux(const char *func, const FNFT_INT line, const char *msg);

#endif
