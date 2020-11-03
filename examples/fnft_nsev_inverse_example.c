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

// This example demonstrates the use of the function fnft_nsev_inverse, which
// implements the inverse nonlinear Fourier transform with respect to the
// nonlinear Schroedinger equation with vanishing boundary conditions.

// The signal is a truncated soliton pulse given in Rourke and Morris,
// Phys. Rev. A 46(7), 1992, https://doi.org/10.1103/PhysRevA.46.3631

#include <stdio.h> // for printf
#include "fnft_nsev_inverse.h"

int main()
{
    FNFT_INT ret_code;

    // The types FNFT_UINT, FNFT_INT, FNFT_REAL and FNFT_COMPLEX are defined
    // in the header fnft_numtypes.h

    /** Step 1: Set up parameters **/

    // Number of samples for the continuous spectrum.
    FNFT_UINT M = 2048;

    // Specify the desired number of samples in the time domain, should be <=M
    FNFT_UINT D = 1024;

    // Specify the location of the first and last time domain sample
    const FNFT_REAL T[2] = { -2.0, 2.0 };

    // We will use the default options when calling fnft_nsev_inverse.
    // See the header file fnft_nsev_inverse.h for possible changes.
    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();

    // Determine the location of the first and last continuous spectrum sample.
    FNFT_REAL XI[2];
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts.discretization);
    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured!\n");
        return EXIT_FAILURE;
    }

    // Focusing nonlinear Schroedinger equation
    int kappa = +1;

    /** Step 2: Allocate memory for the continuous spectrum and the time domain
        signal **/

    FNFT_COMPLEX * contspec = malloc(M * sizeof(FNFT_COMPLEX));
    FNFT_COMPLEX * q = malloc(D * sizeof(FNFT_COMPLEX));
    if (contspec == NULL || q == NULL) {
        printf("An error occured!\n");
        free(contspec);
        free(q);
        return EXIT_FAILURE;
    }

    /** Step 3: Set up the continuous spectrum **/

    const FNFT_REAL alpha = 2;
    const FNFT_REAL beta = 0.55;
    const FNFT_REAL gamma = FNFT_SQRT(FNFT_FABS(alpha)*FNFT_FABS(alpha)
                                      + beta*beta);
    const FNFT_REAL eps_xi = (XI[1] - XI[0])/(M - 1);
    for (FNFT_UINT i=0; i<M; i++) {
        const FNFT_REAL xi = XI[0] + i*eps_xi;
        contspec[i] = alpha / (xi - beta*I);
    }

    /** Step 4: Set up the discrete spectrum **/

    const FNFT_UINT K = 1;
    const FNFT_COMPLEX bound_states[1] = { I*beta };
    const FNFT_COMPLEX normconsts[1] = { -I*alpha/(gamma + beta) };

    /** Step 5: Call fnft_nsev_inverse **/

    ret_code = fnft_nsev_inverse(M, contspec, XI, K, bound_states, normconsts,
                                 D, q, T, kappa, &opts);
    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured!\n");
        free(contspec);
        free(q);
        return EXIT_FAILURE;
    }

    /** Step 6: Print some of the results **/

    FNFT_REAL eps_t = (T[1] - T[0]) / (D - 1);
    printf("Below a few of the %d computed samples are printed:\n",
           (unsigned int)D);
    for (FNFT_UINT i=0; i<D; i+=64) {
        const FNFT_REAL t = T[0] + i*eps_t;
        printf("  q(t=%f) \t= %g + %gI\n",
            (double)t,
            (double)FNFT_CREAL(q[i]),
            (double)FNFT_CIMAG(q[i])
        );
    }

    /** Step 7: Free memory and return **/

    free(contspec);
    free(q);
    return EXIT_SUCCESS;
}
