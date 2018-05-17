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
 * @file fnft__fft_wrapper_plan_t.h
 * @ingroup fft_wrapper
 */

#ifndef FNFT__FFT_WRAPPER_H
#define FNFT__FFT_WRAPPER_H

#include "fnft.h"
#include "kiss_fft.h"

/**
 * @brief Stores information needed by \link fnft__fft_wrapper_execute_plan
 * \endlink to perform a (inverse) FFT.
 * @ingroup fft_wrapper
 */
#ifdef HAVE_FFTW3
#include <fftw3.h>
typedef fftw_plan fnft__fft_wrapper_plan_t;
#else
typedef kiss_fft_cfg fnft__fft_wrapper_plan_t;
#endif

#endif

#ifdef FNFT_ENABLE_SHORT_NAMES
#ifndef FNFT__FFT_WRAPPER_PLAN_T_SHORT_NAMES
#define FNFT__FFT_WRAPPER_PLAN_T_SHORT_NAMES
#define fft_wrapper_plan_t fnft__fft_wrapper_plan_t
#endif
#endif
