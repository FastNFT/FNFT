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

#ifndef FNFT_H
#define FNFT_H

#include "fnft_config.h"
#include "fnft_errwarn.h"

/* Doxygen main page */

/**
 * \mainpage FNFT: Fast nonlinear Fourier transforms
 *
 * This is the documention for the C interface of the library FNFT.
 *
 * You can download FNFT at https://github.com/FastNFT/FNFT/
 *
 * Please read the file README.md first if you are new to FNFT.
 *
 * The C interface is separated into a public part ("fnft_" prefix) and a private
 * part ("fnft__" prefix).
 *
 * To get start with the public part of the C interface, study the modules
 * \ref fnft and \ref fnft_inverse, and try the examples in the examples
 * directory.
 */

/* Define groups for Doxygen documentation */

/**
 * \defgroup fnft Fast nonlinear Fourier transforms
 */

/**
 * \defgroup fnft_inverse Fast inverse nonlinear Fourier transforms
 */

/**
 * \defgroup errwarn Error codes
 */

/**
 * \defgroup data_types Data types
 */

/**
 * \defgroup numtype Macros for numerical operations
 *
 * This module contains macros to work with the numerical data types
 * \link FNFT_INT \endlink, \link FNFT_UINT \endlink, \link FNFT_REAL \endlink
 * and \link FNFT_COMPLEX \endlink. The purpose of these macros is to
 * make changing one of these data types as simple as possible.
 */

/**
 * \defgroup private_errwarn PRIVATE: Macros for raising errors and warnings
 *
 * This module contains macros that can be used to raise errors and warnings.\n
 * \n
 * Most FNFT functions return one of the error codes defined in the module
 * \ref errwarn. The return type should be \link FNFT_INT \endlink. The error
 * should be risen using the corresponding macro from this group. By using
 * these macros, error messages will be printed. The macro returns the
 * corresponding error code. For example:
 * \n \n
 * if (pointer == NULL) \n
 *  return FNFT__E_NOMEM; \n
 * \n
 * The function that prints the errors can be changed using
 * \link fnft_errwarn_setprintf \endlink.
 */

/**
 * \defgroup poly PRIVATE: Polynomials
 */

/**
 * \defgroup misc PRIVATE: Miscellaneous helper functions
 */

/**
 * \defgroup nse PRIVATE: Internals related to the nonlinear \
 *  Schroedinger equation
 */

/**
 * \defgroup kdv PRIVATE: Internals related to the Korteweg-de Vries \
 *  equation
 */

/**
 * \defgroup fft_wrapper PRIVATE: Wrapper for fast Fourier transform routines
 */

#endif
