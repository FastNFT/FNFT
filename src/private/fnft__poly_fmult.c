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

#include "fnft__errwarn.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "kiss_fft.h"
#include "_kiss_fft_guts.h"

static INT poly_fmult2_len(UINT deg)
{
    return kiss_fft_next_fast_size(2*(deg + 1) - 1);
}

static UINT poly_fmult2_lenmen(UINT deg)
{
    return sizeof(kiss_fft_cpx)*(4*poly_fmult2_len(deg) - 1);
}

static INT poly_fmult2(const UINT deg, COMPLEX *p1, \
    COMPLEX *p2, COMPLEX *result, void *mem, \
    kiss_fft_cfg cfg_fft, kiss_fft_cfg cfg_ifft, INT add_flag)
{
    UINT i, len;
    kiss_fft_cpx *buf0, *buf1, *buf2;

    // Prepare buffers
    len = poly_fmult2_len(deg);
    buf0 = (kiss_fft_cpx *)mem;
    buf1 = buf0 + len;
    buf2 = buf1 + len;

    // FFT of first polynomial
    for (i = 0; i <= deg; i++) {
        buf0[i].r = creal(p1[i]);
        buf0[i].i = cimag(p1[i]);
    }
    for (i = deg + 1; i < len; i++) {
        buf0[i].r = 0;
        buf0[i].i = 0;
    }
    kiss_fft(cfg_fft, buf0, buf1);

    // FFT of second polynomial
    for (i = 0; i <= deg; i++) {
        buf0[i].r = creal(p2[i]);
        buf0[i].i = cimag(p2[i]);
    }
    kiss_fft(cfg_fft, buf0, buf2);

    // Inverse FFT of product
    for (i = 0; i < len; i++)
        C_MUL(buf0[i], buf1[i], buf2[i]);
    kiss_fft(cfg_ifft, buf0, buf1);

    // Extract result
    if (!add_flag) {
        for (i = 0; i < 2*deg + 1; i++) {
            result[i] = buf1[i].r/len + I*buf1[i].i/len;
        }
    } else {
        for (i = 0; i < 2*deg + 1; i++) {
            result[i] += buf1[i].r/len + I*buf1[i].i/len;
        }
    }

    // No error
    return SUCCESS;
}

static INT poly_rescale(const UINT d, COMPLEX * const p)
{
    UINT i;
    INT a;
    REAL scl;
    REAL cur_abs;
    REAL max_abs = 0.0;

    // Find max of absolute values of coefficients
    max_abs = 0.0;
    for (i=0; i<=d; i++) {
        cur_abs = CABS( p[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
    }

    // Return if polynomial is identical to zero
    if (max_abs == 0.0)
        return 0;

    // Otherwise, rescale
    a = FLOOR( LOG2(max_abs) );
    scl = POW( 2.0, -a );
    for (i=0; i<=d; i++)
        p[i] *= scl;

    return a;
}

INT fnft__poly_fmult(UINT * const d, UINT n, COMPLEX * const p, 
    INT * const W_ptr)
{
    UINT i, deg, lenmem, len, memneeded, memneeded_buf, mod_n, deg_excess = 0;
    void *mem, *mem_fft, *mem_ifft;
    COMPLEX *p1, *p2, *result;
    kiss_fft_cfg cfg_fft = NULL, cfg_ifft = NULL;
    INT W = 0;
    INT ret_code;

    // Allocate memory for for calls to poly_fmult2
    deg = *d;   
    lenmem = poly_fmult2_lenmen(deg * n);
    mem = malloc(lenmem); // contains the memory for the actual data
    // The line below find max number of bytes needed for an (I)FFT config
    kiss_fft_alloc(poly_fmult2_lenmen(deg * n/2), 0, NULL, &memneeded);
    mem_fft = malloc(memneeded); // memory for FFT configs
    mem_ifft = malloc(memneeded); // memory for IFFT configs
    if (mem == NULL || mem_fft == NULL || mem_ifft == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Main loop, n is the current number of polynomials, deg is their degree
    while (n >= 2) {

        // Create FFT and IFFT config (computes twiddle factors, so reuse)
        len = poly_fmult2_len(deg);
        memneeded_buf = memneeded;
        cfg_fft = kiss_fft_alloc(len, 0, mem_fft, &memneeded_buf);
        memneeded_buf = memneeded;
        cfg_ifft = kiss_fft_alloc(len, 1, mem_ifft, &memneeded_buf);
        if (cfg_fft == NULL || cfg_ifft == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }   

        // Pointers to current pair of polynomials and their product
        p1 = p;
        p2 = p + (deg + 1);
        result = p;
       
        // Multiply all pairs of polynomials, normalize if desired
        mod_n = n % 2;
        for (i=0; i<n-mod_n; i+=2) {

            ret_code = poly_fmult2(deg, p1, p2, result, mem, cfg_fft,
                cfg_ifft, 0);
            CHECK_RETCODE(ret_code, release_mem);

            if (W_ptr != NULL)
                W += poly_rescale(2*deg, result);

            p1 += 2*deg + 2;
            p2 += 2*deg + 2;
            result += 2*deg + 1;
        }

        // If the n number of polynomials was odd, take care of the left over
        // polynomial -> simply keep it, but extend to next degree 2*deg+1
        if (mod_n != 0) {

            INT diff = (n - mod_n)/2;
            memmove(p1 + (deg - diff), p1, (deg+1)*sizeof(COMPLEX));
            for (i=0; i<deg; i++)
                result[i] = 0;
            deg_excess += deg; // to keep track of leading zero coefficients
                               // that we have introduced with the zero padding
        }

        // Double degrees and (more or less) half the number of polynomials
        deg *= 2;
        n = (n - mod_n)/2 + mod_n;
    }

    // Remove leading coefficients that we know are are zero, if any
    if (deg_excess > 0) {
        memmove(p, p + deg_excess, (deg-deg_excess+1)*sizeof(COMPLEX));
        deg -= deg_excess;
    }

    // Set degree of final result, free memory and return w/o error
    *d = deg;
    if (W_ptr != NULL)
        *W_ptr = W;
release_mem:  
    free(mem);
    free(mem_fft);
    free(mem_ifft);
    return ret_code;
}

static INT poly_rescale2x2(const UINT d,
    COMPLEX * const p11,
    COMPLEX * const p12,
    COMPLEX * const p21,
    COMPLEX * const p22)
{
    UINT i;
    INT a;
    REAL scl;
    REAL cur_abs;
    REAL max_abs = 0.0;

    // Find max of absolute values of coefficients
    max_abs = 0.0;
    for (i=0; i<=d; i++) {
        cur_abs = CABS( p11[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p12[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p21[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
        cur_abs = CABS( p22[i] );
        if (cur_abs > max_abs)
            max_abs = cur_abs;
    }

    // Return if polynomials are all identical to zero
    if (max_abs == 0.0)
        return 0;

    // Otherwise, rescale
    a = FLOOR( LOG2(max_abs) );
    scl = POW( 2.0, -a );
    for (i=0; i<=d; i++) {
        p11[i] *= scl;
        p12[i] *= scl;
        p21[i] *= scl;
        p22[i] *= scl;
    }

    return a;
}

/*
* length of p = m*m*n*(deg+1)
* length of result = m*m*(n/2)*(2*deg+1)
* WARNING: p is overwritten
*/
INT fnft__poly_fmult2x2(UINT * const d, UINT n, COMPLEX * const p,
    COMPLEX * const result, INT * const W_ptr)
{
    UINT i, deg, lenmem, len, memneeded, memneeded_buf;
    void *mem, *mem_fft, *mem_ifft;
    UINT o1, o2, or; // pointer offsets
    COMPLEX *p11, *p12, *p21, *p22;
    COMPLEX *r11, *r12, *r21, *r22;
    kiss_fft_cfg cfg_fft = NULL, cfg_ifft = NULL;
    INT W = 0;
    INT ret_code;

    // Setup pointers to the individual polynomials in p
    deg = *d;
    p11 = p;
    p12 = p11 + n*(deg+1);
    p21 = p12 + n*(deg+1);
    p22 = p21 + n*(deg+1);
   
    // Allocate memory for for calls to poly_fmult2
    lenmem = poly_fmult2_lenmen(deg * n);
    mem = malloc(lenmem); // memory for actual data
    // Find max number of bytes needed for an (I)FFT configuration
    kiss_fft_alloc(poly_fmult2_len(*d * n/2), 0, NULL, &memneeded);
    mem_fft = malloc(memneeded); // memory for the FFT configs
    mem_ifft = malloc(memneeded); // memory for the IFFT configs
    if (mem == NULL || mem_fft == NULL || mem_ifft == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Main loop, n is the current number of polynomials, deg is their degree
    while (n >= 2) {

        // Create FFT and IFFT config (computes twiddle factors, so reuse)
        len = poly_fmult2_len(deg);
        memneeded_buf = memneeded;
        cfg_fft = kiss_fft_alloc(len, 0, mem_fft, &memneeded_buf);
        memneeded_buf = memneeded;
        cfg_ifft = kiss_fft_alloc(len, 1, mem_ifft, &memneeded_buf);
        if (cfg_fft == NULL || cfg_ifft == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }   

        // Offsets for the current pair of polynomials and their product
        o1 = 0;
        o2 = deg + 1;
        or = 0;

        // Setup pointers to the individual polynomials in result
        r11 = result;
        r12 = r11 + (n/2)*(2*deg+1);
        r21 = r12 + (n/2)*(2*deg+1);
        r22 = r21 + (n/2)*(2*deg+1);

        // Multiply all pairs of polynomials, normalize if desired
        for (i=0; i<n; i+=2) {

            // Multiply current pair of 2x2 matrix-valued polynomials
            ret_code = poly_fmult2(deg, p11+o1, p11+o2, r11+or, mem,
                cfg_fft, cfg_ifft, 0);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p12+o1, p21+o2, r11+or, mem,
                cfg_fft, cfg_ifft, 1);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p11+o1, p12+o2, r12+or, mem,
                cfg_fft, cfg_ifft, 0);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p12+o1, p22+o2, r12+or, mem, 
                cfg_fft, cfg_ifft, 1);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p21+o1, p11+o2, r21+or, mem,
                cfg_fft, cfg_ifft, 0);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p22+o1, p21+o2, r21+or, mem,
                cfg_fft, cfg_ifft, 1);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p21+o1, p12+o2, r22+or, mem,
                cfg_fft, cfg_ifft, 0);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = poly_fmult2(deg, p22+o1, p22+o2, r22+or, mem,
                cfg_fft, cfg_ifft, 1);
            CHECK_RETCODE(ret_code, release_mem);

            // Normalize if desired
            if (W_ptr != NULL)
                W += poly_rescale2x2(2*deg, r11+or, r12+or, r21+or, r22+or);

            // Move pointers to next pair
            o1 += 2*deg + 2;
            o2 += 2*deg + 2;
            or += 2*deg + 1;
        }

        // Update degrees and number of polynomials
        deg *= 2;
        if (n%2 != 0) {
            ret_code = E_INVALID_ARGUMENT(n); // n was no power of two
            goto release_mem;
        }
        n /= 2;

        // Prepare for the next iteration
        if (n>1) {
            memcpy(p11, r11, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p12, r12, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p21, r21, n*(deg+1)*sizeof(COMPLEX));
            memcpy(p22, r22, n*(deg+1)*sizeof(COMPLEX));
        }
    }
    
    // Set degree of final result, free memory and return w/o error
    *d = deg;
    if (W_ptr != NULL)
        *W_ptr = W;
release_mem:
    free(mem);
    free(mem_fft);
    free(mem_ifft);
    return ret_code;
}

