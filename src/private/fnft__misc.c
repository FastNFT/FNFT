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
* Sander Wahls (TU Delft) 2017.
* Peter J Prins (TU Delft) 2018-2020.
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include <stdio.h>
#include "fnft__misc.h"
#include "fnft__fft_wrapper.h"

void misc_print_buf(const UINT len, COMPLEX const * const buf,
                    char const * const varname)
{
    UINT i;
    printf("%s = [", varname);
    for (i = 0; i < len; i++) {
        printf("%1.12e+%1.12ej", CREAL(buf[i]), CIMAG(buf[i]));
        if (i != len-1)
            printf(", ");
    }
    printf("];\n");
}

REAL misc_rel_err(const UINT len, COMPLEX const * const vec_numer,
    COMPLEX const * const vec_exact)
{
    UINT i;
    double n = 0.0, d = 0.0;
    for (i=0; i<len; i++) {
        n += CABS(vec_numer[i] - vec_exact[i]);
        d += CABS(vec_exact[i]);
    }
    return n/d;
}

REAL misc_hausdorff_dist(const UINT lenA,
    COMPLEX const * const vecA, const UINT lenB,
    COMPLEX const * const vecB)
{
    UINT i, j;
    double tmp, dist, max_dist = -1.0;

    for (i=0; i<lenA; i++) {
        dist = INFINITY;
        for (j=0; j<lenB; j++) {
            tmp = CABS(vecA[i] - vecB[j]);
            if (tmp < dist)
                dist = tmp;
        }
        if (dist > max_dist)
            max_dist = dist;
    }

    for (j=0; j<lenB; j++) {
        dist = INFINITY;
        for (i=0; i<lenA; i++) {
            tmp = CABS(vecA[i] - vecB[j]);
            if (tmp < dist)
                dist = tmp;
        }
        if (dist > max_dist)
            max_dist = dist;
    }

    return max_dist;
}

COMPLEX misc_sech(COMPLEX Z)
{
    return 2.0 / (CEXP(Z) + CEXP(-Z));
}

REAL misc_l2norm2(const UINT N, COMPLEX const * const Z,
    const REAL a, const REAL b)
{
    // Check inputs
    if (a >= b || N==0)
        return NAN;

    // Integrate |q(t)|^2 numerically
    REAL val = 0.0;
    for (UINT i=0; i<N; i++) {
        REAL tmp = (REAL)CABS(Z[i]);
        val += tmp * tmp;
    }
    val *= (b - a)/N;

    return val;
}

INT misc_filter(UINT * const N_ptr, COMPLEX * const vals,
    COMPLEX * const rearrange_as_well,
    REAL const * const bounding_box)
{
    UINT i, N_local, N_filtered;

    if (N_ptr == NULL)
        return E_INVALID_ARGUMENT(N_ptr);
    if (vals == NULL)
        return E_INVALID_ARGUMENT(vals);
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
        || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);

    N_filtered = 0; // Will contain number of values that survived the
                    // filtering (the current no of candidates is in N)
    N_local = *N_ptr;

    for (i=0; i<N_local; i++) {

        // Check if this value is in the bounding box, proceed to next value
        // if not. Due to the formulation, NANs are not in the box.

        if (! (CREAL(vals[i]) >= bounding_box[0]) )
            continue;
        if (! (CREAL(vals[i]) <= bounding_box[1]) )
            continue;
        if (! (CIMAG(vals[i]) >= bounding_box[2]) )
            continue;
        if (! (CIMAG(vals[i]) <= bounding_box[3]) )
            continue;

        // Keep value since no reason to skip it has been found
        vals[N_filtered] = vals[i];
        if (rearrange_as_well != NULL)
            rearrange_as_well[N_filtered] = rearrange_as_well[i];
        N_filtered++;
    }
    *N_ptr = N_filtered;

    return SUCCESS;
}

INT misc_filter_inv(UINT * const N_ptr, COMPLEX * const vals,
    COMPLEX * const rearrange_as_well,
    REAL const * const bounding_box)
{
    UINT i, N_local, N_filtered;
    INT ok_flag;

    if (N_ptr == NULL)
        return E_INVALID_ARGUMENT(N_ptr);
    if (vals == NULL)
        return E_INVALID_ARGUMENT(vals);
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
        || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);

    N_filtered = 0; // Will contain number of values that survived the
                    // filtering (the current no of candidates is in N)
    N_local = *N_ptr;

    for (i=0; i<N_local; i++) {

        // Check if this value is in the bounding box, proceed to next value
        // if so. Due to the formulation, NANs are skipped.
        ok_flag = 0; // not ok
        if (!(CREAL(vals[i]) > bounding_box[0]))
            ok_flag = 1;
        if (!(CREAL(vals[i]) < bounding_box[1]))
            ok_flag = 1;
        if (!(CIMAG(vals[i]) > bounding_box[2]))
            ok_flag = 1;
        if (!(CIMAG(vals[i]) < bounding_box[3]))
            ok_flag = 1;
        if (ok_flag == 1) {
            vals[N_filtered] = vals[i];
            if (rearrange_as_well != NULL)
                rearrange_as_well[N_filtered] = rearrange_as_well[i];
            N_filtered++;
        }
    }
    *N_ptr = N_filtered;

    return SUCCESS;
}

INT misc_filter_nonreal(UINT *N_ptr, COMPLEX * const vals, const REAL tol_im)
{
    UINT i, N_local, N_filtered;

    if (N_ptr == NULL)
        return E_INVALID_ARGUMENT(N_ptr);
    if (vals == NULL)
        return E_INVALID_ARGUMENT(vals);
    if (!(tol_im >= 0))
        return E_INVALID_ARGUMENT(tol_im);

    N_local = *N_ptr;
    N_filtered = 0;
    for (i=0; i<N_local; i++) {
        if (!(FABS(CIMAG(vals[i])) > tol_im))
            continue;
        vals[N_filtered++] = vals[i];
    }
    *N_ptr = N_filtered;

    return SUCCESS;
}

INT misc_merge(UINT *N_ptr, COMPLEX * const vals, REAL tol)
{
    REAL dist = -1.0;
    UINT i, j, N, N_filtered;

    if (N_ptr == NULL)
        return E_INVALID_ARGUMENT(N_ptr);
    if (*N_ptr == 0)
        return SUCCESS;
    if (vals == NULL)
        return E_INVALID_ARGUMENT(vals);
    if (tol < 0.0)
        return E_INVALID_ARGUMENT(tol);

    N = *N_ptr;
    N_filtered = 1;
    for (i=1; i<N; i++) {
        for (j = 0; j < i; j++) {
            dist = CABS(vals[j] - vals[i]);
            if (dist < tol)
                break;
        }
        if (dist < tol)
            continue;

        // Keep bound value since it is not close to previous values
        vals[N_filtered++] = vals[i];
    }
    *N_ptr = N_filtered;

    return SUCCESS;
}

INT misc_downsample(const UINT D, COMPLEX const * const q,
    UINT * const Dsub_ptr, COMPLEX ** qsub_ptr, UINT * const first_last_index)
{
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (D <= 2)
        return E_INVALID_ARGUMENT(D);
    if (qsub_ptr == NULL)
        return E_INVALID_ARGUMENT(qsub_ptr);
    if (Dsub_ptr == NULL)
        return E_INVALID_ARGUMENT(Dsub_ptr);
    if (first_last_index == NULL)
        return E_INVALID_ARGUMENT(first_last_index);

    // Determine number of samples after downsampling, Dsub
    UINT Dsub = *Dsub_ptr; // desired Dsub
    if (Dsub < 2)
       Dsub = 2;
    if (Dsub > D)
        Dsub = D;
    const UINT nskip_per_step = (UINT) ROUND((REAL)D / Dsub);
    Dsub = (UINT) ROUND((REAL)D / nskip_per_step); // actual Dsub

    COMPLEX * const qsub = malloc(Dsub * sizeof(COMPLEX));
    if (qsub == NULL)
        return E_NOMEM;

    // Perform the downsampling
    UINT isub, i = 0;
    for (isub=0; isub<Dsub; isub++) {
        qsub[isub] = q[i];
        i += nskip_per_step;
    }

    // Original index of the first and last sample in qsub
    first_last_index[0] = 0;
    first_last_index[1] = i - nskip_per_step;

    *qsub_ptr = qsub;
    *Dsub_ptr = Dsub;
    return SUCCESS;
}

UINT misc_nextpowerof2(const UINT number)
{
    if (number == 0)
        return 0;
    UINT result = 1;
    while (result < number)
        result *= 2;
    return result;
}

INT misc_resample(const UINT D, const REAL eps_t, COMPLEX const * const q,
    const REAL delta, COMPLEX *const q_new)
{
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (D <= 2)
        return E_INVALID_ARGUMENT(D);
    if (q_new == NULL)
        return E_INVALID_ARGUMENT(q_new);
    if (eps_t == 0)
        return E_INVALID_ARGUMENT(eps_t);

    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    COMPLEX *buf0 = NULL, *buf1 = NULL;
    REAL *freq = NULL;
    INT ret_code;
    UINT i;
    // Allocate memory


    const UINT lenmem = D * sizeof(COMPLEX);
    buf0 = fft_wrapper_malloc(lenmem);
    buf1 = fft_wrapper_malloc(lenmem);
    freq = malloc(D * sizeof(REAL));
    if (buf0 == NULL || buf1 == NULL || freq == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    ret_code = fft_wrapper_create_plan(&plan_fwd, D, buf0, buf1, -1);
    CHECK_RETCODE(ret_code, release_mem);


    ret_code = fft_wrapper_create_plan(&plan_inv, D, buf0, buf1, 1);
    CHECK_RETCODE(ret_code, release_mem);

    // Continuing the signal periodically
    for (i = 0; i < D; i++)
        buf0[i] = q[i];

    ret_code = fft_wrapper_execute_plan(plan_fwd, buf0, buf1);
    CHECK_RETCODE(ret_code, release_mem);

    // Check that the signal spectrum decays sufficiently to ensure accurate interpolation
    
    // D samples of buf1 correspond to 100% bandwidth. Find the total l2-norm
    // of buf1 and compare it to l2-norm of 90% bandwidth. To prevent repeated 
    // computation, the two norms of the 5% bandwidths on both ends of the
    // spectrum are computed instead.
    UINT Dlp = D/20;
    REAL tmp = SQRT(misc_l2norm2(Dlp, buf1+D/2-1-Dlp, 0, Dlp*eps_t) + 
            misc_l2norm2(Dlp, buf1+D/2+1, 0, Dlp*eps_t))/SQRT(misc_l2norm2(D, buf1, 0, D*eps_t));
    if (tmp > SQRT(EPSILON))
        WARN("Signal does not appear to be bandlimited. Interpolation step may be inaccurate. Try to reduce the step size, or switch to a discretization that does not require interpolation");
    
    const REAL scl_factor = (REAL)D*eps_t;
    // Applying phase shift
    for (i = 0; i <D/2; i++)
        freq[i] = i/scl_factor;
    for (i = D/2; i < D; i++)
        freq[i] = ((REAL)i - (REAL)D)/scl_factor;

    for (i = 0; i < D; i++) {
        buf1[i] *= CEXP(2*I*PI*delta*freq[i]);
    }

    // Inverse FFT and truncation
    ret_code = fft_wrapper_execute_plan(plan_inv, buf1, buf0);
    CHECK_RETCODE(ret_code, release_mem);
    for (i = 0; i < D; i++)
        q_new[i] = buf0[i]/D;


release_mem:  
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    fft_wrapper_free(buf0);
    fft_wrapper_free(buf1);
    free(freq);
    return ret_code;
}
