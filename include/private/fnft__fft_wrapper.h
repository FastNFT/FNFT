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
 * @file fnft__fft_wrapper.h
 * @brief Wraps underlying FFT library.
 * 
 * @ingroup fft_wrapper
 *
 * Wraps a FFT library (currently either KISS FFT or, if HAVE_FFTW3 is set by
 * cmake, FFTW3). The function bodies are all declared as static inline and
 * directly included in the header file for speed.
 */

#include "fnft__fft_wrapper_plan_t.h"
#include "fnft__errwarn.h"

/**
 * @brief Next valid number of samples for the FFT routines.
 * @ingroup fft_wrapper
 *
 * Returns the next larger or equal valid number of samples accepted by
 * \link fnft__fft_wrapper_create_plan \endlink.
 * @param[in] desired_length Desired number of samples
 * @return Actual number of samples to be used.
 */
static inline FNFT_UINT fnft__fft_wrapper_next_fft_length(
    FNFT_UINT desired_length)
{
    // Also seems to work well with FFTW
    return kiss_fft_next_fast_size(desired_length);    
}

/**
 * @brief Value to initialize plan variables.
 * @ingroup fft_wrapper
 *
 * To avoid having uninitialized plans, initialize new plans with the value
 * returned by this function (usually, this is simply NULL).
 */
static inline fnft__fft_wrapper_plan_t fnft__fft_wrapper_safe_plan_init()
{
    return NULL;
}

/**
 * @brief Prepares a new (inverse) fast Fourier transform (FFT).
 * @ingroup fft_wrapper
 *
 * Plans can be reused as long as the parameters of the FFT (fft_length and
 * is_inverse) do not change.
 *
 * @param[in,out] plan_ptr Pointer a \link fnft__fft_wrapper_plan_t \endlink
 *   object. Will be changed by the routine.
 * @param[in] fft_length Length of the (inverse) FFT to be computed. Must be
 *   generated using \link fnft__fft_wrapper_next_fft_length \endlink.
 * @param[in,out] in Input buffer with fft_length entries. Initialize after
 *   creating the plan as it might be overwritten. Create with \link
 *   fnft__fft_wrapper_malloc \endlink to ensure correct alignment.
 * @param[in,out] out Output buffer with fft_length enties. Initialize after
 *   creating the plan as it might be overwritten. Create with \link 
 *   fnft__fft_wrapper_malloc \endlink to ensure correct alignment.
 * @param[in] is_inverse -1 => forward FFT, 1 => inverse FFT. Note that the
 *   inverse FFT will not be normalized by the factor 1/fft_length.
 * @return FFT_SUCCESS or an error code.
 */
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
    if (is_inverse != 1 && is_inverse != -1)
        return FNFT__E_INVALID_ARGUMENT(is_inverse);

#ifdef HAVE_FFTW3
    *plan_ptr = fftw_plan_dft_1d(fft_length, in, out, is_inverse, FFTW_ESTIMATE);
#else
    (void)in;
    (void)out;
    *plan_ptr = kiss_fft_alloc(fft_length, (is_inverse+1)/2, NULL, NULL);
#endif

    if (*plan_ptr == NULL)
        return FNFT__E_NOMEM;
    return FNFT_SUCCESS;
}

/**
 * @brief Computes a fast Fourier transform (FFT).
 * @ingroup fft_wrapper
 *
 * @param[in] plan Plan object created with
 *   \link fnft__fft_wrapper_create_plan \endlink
 * @param[in] in Input buffer, not neccessarily the same that was used
 *   when creating the plan. The length however has to be the same. Create
 *   with \link fnft__fft_wrapper_malloc \endlink to ensure correct alignment.
 * @param[out] out Output buffer, not neccessarily the same that was used
 *   when creating the plan. The length however has to be the same. Create
 *   with \link fnft__fft_wrapper_malloc \endlink to ensure correct alignment.
 * @return FFT_SUCCESS or an error code.
 */
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

/**
 * @brief Destroys a FFT plan when it is no longer needed.
 * @ingroup fft_wrapper
 *
 * Frees any memory that was allocated when the plan was created, and sets
 * the value of the plan to \link fnft__fft_wrapper_safe_plan_init \endlink
 * to avoid errors when a plan is destroyed several times.
 *
 * @param[in] plan_ptr Pointer to a plan object created with
 *   \link fnft__fft_wrapper_create_plan \endlink
 * @return FFT_SUCCESS or an error code.
 */
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

/**
 * @brief Memory allocation for the FFT wrapper.
 * @ingroup fft_wrapper
 *
 * Use this routine to allocate memory for the input and output buffers. It
 * ensures proper alignment.
 * @param[in] size Number of points to be allocated.
 * @return A pointer to the allocated buffer, or NULL if that failed.
 */
static inline void * fnft__fft_wrapper_malloc(FNFT_UINT size)
{
#ifdef HAVE_FFTW3
    return fftw_malloc(size);
#else
    return KISS_FFT_MALLOC(size);
#endif
}

/**
 * @brief Memory deallocation for the FFT wrapper.
 * @ingroup fft_wrapper
 *
 * Use this routine to deallocate memory that was allocated using
 * \link fnft__fft_wrapper_malloc \endlink.
 *
 * @param[in] ptr Pointer to memory block that is to be freed (or NULL).
 */
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
#define fft_wrapper_next_fft_length(...) fnft__fft_wrapper_next_fft_length(__VA_ARGS__)
#define fft_wrapper_safe_plan_init(...) fnft__fft_wrapper_safe_plan_init(__VA_ARGS__)
#define fft_wrapper_create_plan(...) fnft__fft_wrapper_create_plan(__VA_ARGS__)
#define fft_wrapper_execute_plan(...) fnft__fft_wrapper_execute_plan(__VA_ARGS__)
#define fft_wrapper_destroy_plan(...) fnft__fft_wrapper_destroy_plan(__VA_ARGS__)
#define fft_wrapper_malloc(...) fnft__fft_wrapper_malloc(__VA_ARGS__)
#define fft_wrapper_free(...) fnft__fft_wrapper_free(__VA_ARGS__)
#endif
#endif

