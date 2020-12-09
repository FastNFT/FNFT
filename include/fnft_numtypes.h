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
 * @file fnft_numtypes.h
 * Defines the numerical data types used internally by FNFT and
 * provides support macros.
 * @ingroup numtype
 */

#ifndef FNFT_NUMTYPES_H
#define FNFT_NUMTYPES_H

#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * The floating point data type used by FNFT.
 * @ingroup numtype
 */
typedef double FNFT_REAL;

/**
 * The complex floating point data type used by FNFT.
 * @ingroup numtype
 */
#ifndef __cplusplus
typedef double complex FNFT_COMPLEX;
#else
typedef std::complex<double> FNFT_COMPLEX;
#endif

/**
 * The signed integer used by FNFT.
 * @ingroup numtype
 */
typedef int32_t FNFT_INT;

/**
 * The unsigned integer used by FNFT.
 * @ingroup numtype
 */
typedef size_t FNFT_UINT;

/**
 * Machine precision of \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_EPSILON DBL_EPSILON

/**
 * Not-a-number for \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_NAN NAN

/**
 * INFINITY for \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_INF INFINITY

/**
 * Absolute value of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_FABS(X) fabs(X)

/**
 * Square root of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_SQRT(X) sqrt(X)

/**
 * Cosine of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_COS(X) cos(X)

/**
 * Sine of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_SIN(X) sin(X)

/**
 * Sinh of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_SINH(X) sinh(X)

/**
 * Cosh of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_COSH(X) cosh(X)

/**
 * Tanh of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_TANH(X) tanh(X)

/**
 * Arc tangent of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_ATAN(X) atan(X)

/**
 * Natural logarithm of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_LOG(X) log(X)

/**
 * Base two logarithm of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_LOG2(X) log2(X)

/**
 * Power X^Y of two \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_POW(X, Y) pow(X, Y)

/**
 * Gamma function of a \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_GAMMA(X) tgamma(X)

/**
 * \link FNFT_REAL \endlink) that contains \f$ \pi \f$ up to machine
 * precision.
 * @ingroup numtype
 */
#define FNFT_PI acos(-1.0)

/**
 * Rounds a \link FNFT_REAL \endlink) to the nearest smaller or equal integer
 * (expressed as a \link FNFT_REAL \endlink).
 * @ingroup numtype
 */
#define FNFT_FLOOR(X) floor(X)

/**
 * Rounds a \link FNFT_REAL \endlink) to the nearest integer (expressed as a
 * \link FNFT_REAL \endlink).
 * @ingroup numtype
 */
#define FNFT_ROUND(X) round(X)

/**
 * Rounds a \link FNFT_REAL \endlink) to the nearest larger or equal integer
 * (expressed as a \link FNFT_REAL \endlink).
 * @ingroup numtype
 */
#define FNFT_CEIL(X) ceil(X)

/**
 * Computes sqrt(X*X+Y*Y) for two \link FNFT_REAL \endlink.
 * @ingroup numtype
 */
#define FNFT_HYPOT(X,Y) hypot(X,Y)

/**
 * Real part of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CREAL(X) creal(X)

/**
 * Imaginary part of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CIMAG(X) cimag(X)

/**
 * Absolute value of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CABS(X) cabs(X)

/**
 * Argument of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CARG(X) carg(X)

/**
 * Complex conjugate of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CONJ(X) conj(X)

/**
 * Power X^Y of two \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CPOW(X,Y) cpow(X,Y)

/**
 * Complex exponential of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CEXP(X) cexp(X)

/**
 * Complex natural logarithm of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CLOG(X) clog(X)

/**
 * Complex square root of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CSQRT(X) csqrt(X)

/**
 * Complex hyperbolic sine of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CSINH(X) csinh(X)

/**
 * Complex hyperbolic cosine of a \link FNFT_COMPLEX \endlink.
 * @ingroup numtype
 */
#define FNFT_CCOSH(X) ccosh(X)

/**
 * Complex hyperbolic sine of a \link FNFT_COMPLEX \endlink.
 */
#define FNFT_CSIN(X) csin(X)

/**
 * Complex hyperbolic cosine of a \link FNFT_COMPLEX \endlink.
 */
#define FNFT_CCOS(X) ccos(X)

/**
 * Complex arc hyperbolic tangent of a \link FNFT_COMPLEX \endlink.
 */
#define FNFT_ATANH(X) atanh(X)

#ifdef FNFT_ENABLE_SHORT_NAMES
#define REAL            FNFT_REAL
#define COMPLEX         FNFT_COMPLEX
#define INT             FNFT_INT
#define UINT            FNFT_UINT
#define CABS(X)         FNFT_CABS(X)
#define FABS(X)         FNFT_FABS(X)
#define FLOOR(X)        FNFT_FLOOR(X)
#define CEIL(X)         FNFT_CEIL(X)
#define ROUND(X)        FNFT_ROUND(X)
#define POW(X,Y)        FNFT_POW(X,Y)
#define CPOW(X,Y)       FNFT_CPOW(X,Y)
#define LOG2(X)         FNFT_LOG2(X)
#define LOG(X)          FNFT_LOG(X)
#define CLOG(X)         FNFT_CLOG(X)
#define COS(X)          FNFT_COS(X)
#define SIN(X)          FNFT_SIN(X)
#define ATAN(X)         FNFT_ATAN(X)
#define SQRT(X)         FNFT_SQRT(X)
#define EPSILON         FNFT_EPSILON
#define CSINH(X)        FNFT_CSINH(X)
#define CCOSH(X)        FNFT_CCOSH(X)
#define CCOS(X)         FNFT_CCOS(X)
#define CSIN(X)         FNFT_CSIN(X)
#define COSH(X)         FNFT_COSH(X)
#define SINH(X)         FNFT_SINH(X)
#define TANH(X)         FNFT_TANH(X)
#define HYPOT(X,Y)      FNFT_HYPOT(X,Y)
#define CREAL(X)        FNFT_CREAL(X)
#define CIMAG(X)        FNFT_CIMAG(X)
#define CONJ(X)         FNFT_CONJ(X)
#define CSQRT(X)        FNFT_CSQRT(X)
#define CEXP(X)         FNFT_CEXP(X)
#define CARG(X)         FNFT_CARG(X)
#define PI        	    FNFT_PI
#define GAMMA(X)        FNFT_GAMMA(X)
#define ATANH(X)   FNFT_ATANH(X)
#endif

#endif
