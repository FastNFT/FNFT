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
#include "kiss_fft.h"
#include "_kiss_fft_guts.h"

/*
 * result should be of length M
 * Z = A * W.^-(0:(M-1)); result = polyval(p, 1./Z).'
 * L.R. Rabiner, R.W. Schafer and C.M. Rader, "The Chirp z-Transform
 * Algorithm," IEEE Trans. Audio Electroacoust. 17(2), Jun. 1969.
 */
INT poly_chirpz(const UINT deg, COMPLEX const * const p, \
    const COMPLEX A, const COMPLEX W, const UINT M, \
    COMPLEX * const result)
{
    COMPLEX Z;
    kiss_fft_cpx *Y, *V, *buf;
    kiss_fft_cfg cfg;
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
    const UINT L = kiss_fft_next_fast_size(N + M - 1);
    cfg = kiss_fft_alloc((int)L, 0, NULL, NULL);
    Y = malloc(L * sizeof(kiss_fft_cpx));
    V = malloc(L * sizeof(kiss_fft_cpx));
    buf = malloc(L * sizeof(kiss_fft_cpx));
    if (cfg == NULL && Y == NULL && V == NULL && buf == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Setup yn and compute Yr = fft(yn)
    for (n=0; n<=N-1; n++) {
        Z = p[deg - n] * CPOW(A, -1.0*n) * CPOW(W, 0.5*n*n);
        buf[n].r = CREAL(Z);
        buf[n].i = CIMAG(Z);
    }
    for (n=N; n<L; n++) {
        buf[n].r = 0;
        buf[n].i = 0;
    }
    kiss_fft(cfg, buf, Y);

    // Setup vn and compute Vr = fft(vn)
    for (n=0; n<=M-1; n++) {
        Z = CPOW(W, -0.5*n*n);
        buf[n].r = CREAL(Z);
        buf[n].i = CIMAG(Z);
    }
    for (n=M; n<=L-N; n++) {
        buf[n].r = 0;
        buf[n].i = 0;
    }
    for (n=L-N+1; n<L; n++) {
         Z = CPOW(W, -0.5*(L - n)*(L - n));
         buf[n].r = CREAL(Z);
         buf[n].i = CIMAG(Z);
    }
    kiss_fft(cfg, buf, V);

    // Multiply V and Y
    for (n=0; n<L; n++)
        C_MUL(buf[n], V[n], Y[n]);
    
    // Compute inverse FFT of the product and store it in V
    free(cfg);   
    cfg = kiss_fft_alloc((int)L, 1, NULL, NULL);
    if (cfg == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }
    kiss_fft(cfg, buf, V);

    // Form the final result
    for (n=0; n<M; n++)
        result[n] = CPOW(W, 0.5*n*n) * (V[n].r + I*V[n].i) / L;

    // Release memory and return
release_mem:
    free(cfg);
    free(Y);
    free(V);
    free(buf);
    return ret_code;
}
