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

#ifndef FNFT__FFT_WRAPPER_H
#define FNFT__FFT_WRAPPER_H

#include "fnft.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"
#include "kiss_fft.h"

#ifdef HAVE_FFTW3
#include <fftw3.h>
typedef fftw_plan fnft__fft_wrapper_plan_t;
#else
typedef kiss_fft_cfg fnft__fft_wrapper_plan_t;
#endif

#endif

static inline FNFT_UINT fnft__fft_wrapper_next_fft_length(
    FNFT_UINT desired_length)
{
    // Also seems to work well with FFTW
    return kiss_fft_next_fast_size(desired_length);    
}

static inline fnft__fft_wrapper_plan_t fnft__fft_wrapper_safe_plan_init()
{
    return NULL;
}

static inline FNFT_INT fnft__fft_wrapper_create_plan(
    fnft__fft_wrapper_plan_t * plan_ptr,
    FNFT_UINT fft_length,
    FNFT_COMPLEX * in,
    FNFT_COMPLEX * out,
    FNFT_INT is_inverse)
{
    if (plan_ptr == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan);
    if (fft_length == 0)
        return FNFT__E_INVALID_ARGUMENT(fft_length);
    if (is_inverse > 1)
        return FNFT__E_INVALID_ARGUMENT(is_inverse);

#ifdef HAVE_FFTW3
    if (is_inverse == 0)
        *plan_ptr = fftw_plan_dft_1d(fft_length, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    else
        *plan_ptr = fftw_plan_dft_1d(fft_length, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
    (void)in;
    (void)out;
    *plan_ptr = kiss_fft_alloc(fft_length, is_inverse, NULL, NULL);
#endif

    if (*plan_ptr == NULL)
        return FNFT__E_NOMEM;
    return FNFT_SUCCESS;
}

static inline FNFT_INT fnft__fft_wrapper_execute_plan(
    fnft__fft_wrapper_plan_t plan, FNFT_COMPLEX *in, FNFT_COMPLEX *out)
{
    if (plan == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan);

#ifdef HAVE_FFTW3
    fftw_execute_dft((fftw_plan)plan, (fftw_complex *)in, (fftw_complex *)out);
#else    
    kiss_fft((kiss_fft_cfg)plan, (kiss_fft_cpx *)in, (kiss_fft_cpx *)out);
#endif

    return FNFT_SUCCESS;
}

static inline FNFT_INT fnft__fft_wrapper_destroy_plan(
    fnft__fft_wrapper_plan_t * plan_ptr)
{
    if (plan_ptr == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan_ptr);
#ifdef HAVE_FFTW3    
    fftw_destroy_plan(*plan_ptr);
#else
    KISS_FFT_FREE(*plan_ptr);
#endif    
    *plan_ptr = fnft__fft_wrapper_safe_plan_init();
    return FNFT_SUCCESS;
}

static inline void * fnft__fft_wrapper_malloc(FNFT_UINT size)
{
#ifdef HAVE_FFTW3
    return fftw_malloc(size);
#else
    return KISS_FFT_MALLOC(size);
#endif
}

static inline void fnft__fft_wrapper_free(void * ptr)
{
#ifdef HAVE_FFTW3
    fftw_free(ptr);
#else
    KISS_FFT_FREE(ptr);
#endif
}

#ifdef FNFT_ENABLE_SHORT_NAMES
#ifndef FNFT__FFT_WRAPPER_SHORT_NAMES
#define FNFT__FFT_WRAPPER_SHORT_NAMES
#define fft_wrapper_plan_t fnft__fft_wrapper_plan_t
#define fft_wrapper_next_fft_length(...) fnft__fft_wrapper_next_fft_length(__VA_ARGS__)
#define fft_wrapper_safe_plan_init(...) fnft__fft_wrapper_safe_plan_init(__VA_ARGS__)
#define fft_wrapper_create_plan(...) fnft__fft_wrapper_create_plan(__VA_ARGS__)
#define fft_wrapper_execute_plan(...) fnft__fft_wrapper_execute_plan(__VA_ARGS__)
#define fft_wrapper_destroy_plan(...) fnft__fft_wrapper_destroy_plan(__VA_ARGS__)
#define fft_wrapper_malloc(...) fnft__fft_wrapper_malloc(__VA_ARGS__)
#define fft_wrapper_free(...) fnft__fft_wrapper_free(__VA_ARGS__)
#endif
#endif

