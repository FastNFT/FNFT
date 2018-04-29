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
#include "kiss_fft.h"

typedef struct {
    kiss_fft_cfg cfg;
    FNFT_UINT fft_length;
    FNFT_COMPLEX * in;
    FNFT_COMPLEX * out;
    FNFT_UINT work_length;
    FNFT_COMPLEX * work;   
} fnft__fft_wrapper_plan_t;

#endif

static inline FNFT_UINT fnft__fft_wrapper_next_fft_length(
    FNFT_UINT desired_length)
{
    return kiss_fft_next_fast_size(desired_length);    
}

static inline FNFT_UINT fnft__fft_wrapper_work_size(FNFT_UINT fft_length)
{
    return 0;
}

static inline FNFT_INT fnft__fft_wrapper_create_plan(
    fnft__fft_wrapper_plan_t * plan,
    FNFT_UINT fft_length,
    FNFT_COMPLEX * in,
    FNFT_COMPLEX * out,
    FNFT_INT is_inverse,    
    void * work)
{
    if (plan == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan);
    if (fft_length == 0)
        return FNFT__E_INVALID_ARGUMENT(fft_length);
    if (in == NULL)
        return FNFT__E_INVALID_ARGUMENT(in);
    if (out == NULL)
        return FNFT__E_INVALID_ARGUMENT(out);
    if (is_inverse > 1)
        return FNFT__E_INVALID_ARGUMENT(is_inverse);
    plan->cfg = kiss_fft_alloc(fft_length, is_inverse, NULL, NULL);
    if (plan->cfg == NULL)
        return FNFT__E_NOMEM;
    plan->in = in;
    plan->out = out;
    return FNFT_SUCCESS;
}

static inline FNFT_INT fnft__fft_wrapper_execute_plan(
    fnft__fft_wrapper_plan_t * plan)
{
    if (plan == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan);
    kiss_fft(plan->cfg, (kiss_fft_cpx *)plan->in, (kiss_fft_cpx *)plan->out);
    return FNFT_SUCCESS;
}

static inline FNFT_INT fnft__fft_wrapper_destroy_plan(
    fnft__fft_wrapper_plan_t * plan)
{
    if (plan == NULL)
        return FNFT__E_INVALID_ARGUMENT(plan);
    KISS_FFT_FREE(plan->cfg);
    plan->cfg = NULL;
    return FNFT_SUCCESS;
}

static inline void * fnft__fft_wrapper_malloc(FNFT_UINT size)
{
    return KISS_FFT_MALLOC(size);
}

static inline void fnft__fft_wrapper_free(void * ptr)
{
    KISS_FFT_FREE(ptr);
}

static inline INT fnft__fft_wrapper_single_fft(UINT len, COMPLEX * in,
    COMPLEX * out, FNFT_INT is_inverse, void * work)
{
    FNFT_INT ret_code = FNFT_SUCCESS;
    fnft__fft_wrapper_plan_t plan;
    ret_code = fnft__fft_wrapper_create_plan(&plan, len, in, out, is_inverse,
        work);
    FNFT__CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__fft_wrapper_execute_plan(&plan);
    FNFT__CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__fft_wrapper_destroy_plan(&plan);
    FNFT__CHECK_RETCODE(ret_code, leave_fun);
leave_fun:
    return ret_code;
}

#ifdef FNFT_ENABLE_SHORT_NAMES
#ifndef FNFT__FFT_WRAPPER_SHORT_NAMES
#define FNFT__FFT_WRAPPER_SHORT_NAMES
#define fft_wrapper_plan_t fnft__fft_wrapper_plan_t
#define fft_wrapper_next_fft_length(...) fnft__fft_wrapper_next_fft_length(__VA_ARGS__)
#define fft_wrapper_work_size(...) fnft__fft_wrapper_work_size(__VA_ARGS__)
#define fft_wrapper_create_plan(...) fnft__fft_wrapper_create_plan(__VA_ARGS__)
#define fft_wrapper_execute_plan(...) fnft__fft_wrapper_execute_plan(__VA_ARGS__)
#define fft_wrapper_destroy_plan(...) fnft__fft_wrapper_destroy_plan(__VA_ARGS__)
#define fft_wrapper_malloc(...) fnft__fft_wrapper_malloc(__VA_ARGS__)
#define fft_wrapper_free(...) fnft__fft_wrapper_free(__VA_ARGS__)
#define fft_wrapper_single_fft(...) fnft__fft_wrapper_single_fft(__VA_ARGS__)
#endif
#endif

