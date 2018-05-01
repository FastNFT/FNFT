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

#include "fnft.h"
#include "fnft__errwarn.h"
#include "fnft__poly_chirpz.h"
#include "fnft__poly_fmult.h"
#include "fnft__fft_wrapper.h"

/*
 * result should be of length M
 * Z = A * W.^-(0:(M-1)); result = polyval(p, 1./Z).'
 * L.R. Rabiner, R.W. Schafer and C.M. Rader, "The Chirp z-Transform
 * Algorithm," IEEE Trans. Audio Electroacoust. 17(2), Jun. 1969.
 */
INT poly_chirpz(const UINT deg, COMPLEX const * const p,
    const COMPLEX A, const COMPLEX W, const UINT M,
    COMPLEX * const result)
{
    COMPLEX *Y, *V, *buf;
    INT ret_code = SUCCESS;
    UINT n;

    // Check inputs
    if (p == NULL)
        return E_INVALID_ARGUMENT(p);
    if (M == 0)
        return E_INVALID_ARGUMENT(M);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);

    // Allocate memory
    const UINT N = deg + 1;
    const UINT L = fft_wrapper_next_fft_length(N + M - 1);
    Y = fft_wrapper_malloc(L * sizeof(COMPLEX));
    V = fft_wrapper_malloc(L * sizeof(COMPLEX));
    buf = fft_wrapper_malloc(L * sizeof(COMPLEX));
    if (Y == NULL || V == NULL || buf == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Setup yn and compute Yr = fft(yn)
    for (n=0; n<=N-1; n++)
        buf[n] = p[deg - n] * CPOW(A, -1.0*n) * CPOW(W, 0.5*n*n);
    for (n=N; n<L; n++) 
        buf[n] = 0;
    ret_code = fft_wrapper_single_fft(L, buf, Y, 0);
    CHECK_RETCODE(ret_code, release_mem);

    // Setup vn and compute Vr = fft(vn)
    for (n=0; n<=M-1; n++)
        buf[n] = CPOW(W, -0.5*n*n);
    for (n=M; n<=L-N; n++)
        buf[n] = 0;
    for (n=L-N+1; n<L; n++)
         buf[n] = CPOW(W, -0.5*(L - n)*(L - n));
    ret_code = fft_wrapper_single_fft(L, buf, V, 0);
    CHECK_RETCODE(ret_code, release_mem);

    // Multiply V and Y
    for (n=0; n<L; n++)
        buf[n] = V[n] * Y[n];
    
    // Compute inverse FFT of the product and store it in V
    ret_code = fft_wrapper_single_fft(L, buf, V, 1);
    CHECK_RETCODE(ret_code, release_mem);

    // Form the final result
    for (n=0; n<M; n++)
        result[n] = CPOW(W, 0.5*n*n) * V[n] / L;

    // Release memory and return
release_mem:
    fft_wrapper_free(Y);
    fft_wrapper_free(V);
    fft_wrapper_free(buf);
    return ret_code;
}
