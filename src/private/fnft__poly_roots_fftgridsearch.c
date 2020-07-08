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
* Marius Brehler (TU Dortmund) 2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_roots_fftgridsearch.h"
#include "fnft__poly_chirpz.h"
#ifdef DEBUG
#include "fnft__misc.h" // for misc_filter
#include <stdio.h> // for printf
#endif

// Computation of polynomial roots on the unit circle via gridsearch.
// *M_ptr is the number of points in the grid. The array roots
// must be preallocated by the user with *M_ptr entries. Upon exit, *M_ptr
// contains the number of found roots, which are located in the beginning of
// the array roots. Returns SUCCESS or an error code.
INT poly_roots_fftgridsearch(const UINT deg,
    COMPLEX const * const p, UINT * const M_ptr,
    REAL const * const PHI, COMPLEX * const roots)  
{
    INT ret_code;
    COMPLEX A, W, c, zi, z0, yi, y0, zr;
    UINT i, j, M, nroots = 0;
    INT k;
    COMPLEX * vals;
    REAL tmp;

	// Check inputs
    if ( deg < 2 )
        return E_INVALID_ARGUMENT(deg);
    if (p == NULL)
	return E_INVALID_ARGUMENT(p);
    if (M_ptr == NULL || *M_ptr < 2)
 	return E_INVALID_ARGUMENT(M_ptr);      
    if (PHI == NULL || !(PHI[0] < PHI[1]) || PHI[0] == -INFINITY
    || PHI[1] == INFINITY)
        return E_INVALID_ARGUMENT(PHI);
    if (roots == NULL)
	return E_INVALID_ARGUMENT(roots);

    // Allocate memory
    M = *M_ptr;
    vals = malloc(3*M * sizeof(COMPLEX));
    if (vals == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Evaluate polynomial using the Chirp transform on three rings
    const REAL eps = (PHI[1] - PHI[0]) / (M - 1);
    W = CEXP(I*eps);
    for (k=-1; k<=1; k++) {

        A = (1.0 + k*eps) * CEXP(-I*PHI[0]);
        ret_code = poly_chirpz(deg, p, A, W, M, vals + (k+1)*M);
        CHECK_RETCODE(ret_code, release_mem);
    }
   
    // Approximate the roots
    for (i=1; i<M-1; i++) {

        // Minimum modulus theorem => minimum absolute value must be on the
        // boundary of a domain around the current test point unless there is a
        // root. Keep track of absolute values of moving grid of nine points.
        tmp = CABS( vals[M + i] );
        if ( tmp > CABS(vals[i - 1]) )
            continue;
        if ( tmp > CABS(vals[i]) )
            continue;
        if ( tmp > CABS(vals[i + 1]) )
            continue;
        if ( tmp > CABS(vals[M + i - 1]) )
            continue;
        if ( tmp > CABS(vals[M + i + 1]) )
            continue;
        if ( tmp > CABS(vals[2*M + i - 1]) )
            continue;
        if ( tmp > CABS(vals[2*M + i]) )
            continue;
        if ( tmp > CABS(vals[2*M + i + 1]) )
            continue;

        // Let z0 be the center point of the current grid such that
        // y0 = p(z0) = vals[M+i]. We approximate p(z) locally around z0 with
        // a linear function: p(z) ~= y0 + c(z-z0). The coefficient c
        // is found via a least squares fit w.r.t. the other vals, i.e., by
        // minimizing ||[y1-y0 ... yn-y0]-c[z1-z0 ... zn-z0]||^2.
        z0 = CEXP(I*(PHI[0] + i*eps));
        c = 0.0;
        tmp = 0.0;
        y0 = vals[M + i];
        for (j=i-1; j<i+2; j++) {
            for (k=-1; k<2; k++) {

                if (j == 0 && k == 0)
                    continue; // Skip the center point

                zi = (1 - k*eps)*CEXP(I*(PHI[0] + j*eps));
                yi = vals[(k + 1)*M + j];
                c += CONJ( zi - z0 )*( yi - y0 );
                tmp += CABS( zi - z0 )*CABS( zi - z0 );
            }
        }
        if (tmp == 0.0) { // This should never happen ...
            ret_code = E_DIV_BY_ZERO;
            goto release_mem;
        }
        c /= tmp;

        // Find the root zr of the linear approximation y0+c(z-z0)
        if (c == 0.0) { // Pathologic case: linear apprimation has no roots
                        // at all or is zero everywhere

            if (y0 != 0.0)
                continue;
            zr = z0;

        } else { // Linear approximation has a unique root

            // We solve 0=y0+c(z-z0) for z.
            zr = z0 - y0/c;
            if ( CABS(zr - z0) > eps )
                continue; // root estimate is too far from center

        }

        // Save the root
        roots[nroots++] = zr;
    }
    // Save the number of detected roots
    *M_ptr = nroots;

release_mem:
    free(vals);
    return ret_code;
}

// Computation of polynomial roots on the unit circle via gridsearch.
// The degree must be odd, and p(z) * z^(N-1), where N=deg/2+1, must be
// para-hermitian. *M_ptr is the number of points in the grid. The array roots
// must be preallocated by the user with *M_ptr entries. Upon exit, *M_ptr
// contains the number of found roots, which are located in the beginning of
// the array roots. Returns SUCCESS or an error code.
INT poly_roots_fftgridsearch_paraherm(const UINT deg,
    COMPLEX const * const p, UINT * const M_ptr,
    REAL const * const PHI, COMPLEX * const roots)  
{
    INT ret_code;
    COMPLEX A, W;
    REAL phi, phi1, phi2;
    UINT i, N, M, nroots = 0;

	// Check inputs
    if ( deg%2 == 1 || deg < 2 ) // degree must be even and >= 2
        return E_INVALID_ARGUMENT(deg);
    if (p == NULL)
	return E_INVALID_ARGUMENT(p);
    if (M_ptr == NULL || *M_ptr < 2)
 		return E_INVALID_ARGUMENT(M_ptr);      
    if (PHI == NULL || !(PHI[0] < PHI[1]) || PHI[0] == -INFINITY
    || PHI[1] == INFINITY)
        return E_INVALID_ARGUMENT(PHI);
	if (roots == NULL)
		return E_INVALID_ARGUMENT(roots);

    // Evaluate polynomial using the Chirp transform
    M = *M_ptr;
    const REAL eps = (PHI[1] - PHI[0]) / (M - 1);
    W = CEXP(I*eps);
    A = CEXP(-I*PHI[0]);
    ret_code = poly_chirpz(deg, p, A, W, M, roots);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);

    // Remove the phase factor
    N = deg/2 + 1; // remember: we checked that deg even
    for (i=0; i<M; i++) {
        phi = PHI[0] + eps*i; // angle of roots[i]
        roots[i] *= CEXP(-I*phi*(N-1));
    }

    // Approximate the roots
    for (i=1; i<M; i++) {
        // Mean value theorem => change of sign between two consecutive values
        // shows there is a root in between them
        if (CREAL(roots[i-1])*CREAL(roots[i]) <= 0.0) {
            // We approximate the polynomial linearly to estimate the location
            // of the root 
            phi1 = PHI[0] + eps*(i - 1);
            phi2 = phi1 + eps;
            if (roots[i-1] != roots[i])
                phi = phi1 - roots[i-1]*(phi2 - phi1)/(roots[i] - roots[i-1]);
            else
                phi = 0.5*(phi1 + phi2); // angle of the root estimate
            roots[nroots++] = CEXP(I*phi);
        }
    }

    *M_ptr = nroots;
    return SUCCESS;
}
