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

// This example demonstrates the use of the function fnft_nsep, which 
// implements the nonlinear Fourier transform with respect to the nonlinear
// Schroedinger equation with periodic boundary conditions. The signal is the
// same plane wave as in Sec. II of https://doi.org/10.1103/PhysRevA.37.815

#include <stdio.h> // for printf
#include "fnft_nsep.h"

int main()
{
    /** Step 1: Set up the signal **/

    // Number of time-domain samples. Increase to improve precision of results.
    FNFT_UINT D = 256;

    // Contains the samples of q(t)
    FNFT_COMPLEX q[D];

    // Location of the 1st time-domain sample and the beginning of the next
    // period. Location of last sample is T[1]. In previous version it was T[1]-eps_t
    FNFT_REAL T[2] = { 0.0, 2.0*FNFT_PI };

    // Define a simple rectangular signal
    for (FNFT_UINT i=0; i<D; i++) {
        FNFT_REAL t = T[0] + i*(T[1]-T[0])/(D-1);
        q[i] = FNFT_CEXP(2.0*I*t);
    }

    // The types FNFT_UINT, FNFT_INT, FNFT_REAL and FNFT_COMPLEX as well as
    // the macro FNFT_CEXP are defined in the header fnft_numtypes.h

    /** Step 2: Prepare calling fnft_nsep **/ 

    // Max. number of main spectrum points we expect = size of the array below
    FNFT_UINT K = D;

    // Buffer for result: Main spectrum
    FNFT_COMPLEX main_spec[K];

    // Max. number of aux. spectrum points we expect = size of the array below
    FNFT_UINT M = D;

    // Buffer for result: Auxiliary spectrum
    FNFT_COMPLEX aux_spec[M];

    // Focusing nonlinear Schroedinger equation
    int kappa = +1;

    // Default options
    fnft_nsep_opts_t opts = fnft_nsep_default_opts();

    // We only consider spectrum in the box -2<=Re(lambda)<=2,
    // -2<=Im(lambda)<=2.
    opts.filtering = fnft_nsep_filt_MANUAL;
    opts.bounding_box[0] = -2;
    opts.bounding_box[1] = 2;
    opts.bounding_box[2] = -2;
    opts.bounding_box[3] = 2;

    // See the header file fnft_nsep.h for other options

    /** Step 3: Call fnft_nsev and check for errors **/

    int ret_code = fnft_nsep(D, q, T, &K, main_spec, &M, aux_spec, NULL, kappa,
        &opts);
    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured!\n");
        return EXIT_FAILURE;
    }

    /** Step 4: Print the results. **/

    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    printf("Main spectrum:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  %g + %gI\n",
            (double)FNFT_CREAL(main_spec[i]),
            (double)FNFT_CIMAG(main_spec[i])
        );
    }

    printf("Auxiliary spectrum:\n");
    for (FNFT_UINT i=0; i<M; i++) {
        printf("  %g + %gI\n",
            (double)FNFT_CREAL(aux_spec[i]),
            (double)FNFT_CIMAG(aux_spec[i])
        );
    }

    // Note that also degenerate spectra, which manifest themselves in
    // repeated points, are printed. When comparing the results with the plot
    // in the paper mentioned above, please take into consideration that lambda
    // here corresponds to k in the paper.

    return EXIT_SUCCESS;
}
