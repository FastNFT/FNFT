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
 * Peter J. Prins (TU Delft) 2021.
 */

// This example demonstrates the use of the function fnft_kdvv, which 
// implements the nonlinear Fourier transform with respect to the
// Korteweg-de Vries equation with vanishing boundary conditions. The signal
// is a rectangular pulse.

#include <stdio.h> // for printf
#include "fnft_kdvv.h"

int main()
{
    /** Step 1: Set up the signal **/

    // Number of time-domain samples. Increase to improve precision of results.
    FNFT_UINT D = 256;

    // Contains the samples of q(t)
    FNFT_COMPLEX q[D];

    // Location of the 1st and last time-domain sample
    FNFT_REAL T[2] = { -1.0, 1.0 };

    // Define a simple rectangular signal
    for (FNFT_UINT i=0; i<D; i++) {
        // The line below would determine the location of the current sample
        //FNFT_REAL t = T[0] + i*(T[1]-T[0])/(D-1);
        q[i] = 2.0 + 0.0*I;
    }

    // The types FNFT_UINT, FNFT_INT, FNFT_REAL and FNFT_COMPLEX are defined
    // in the header fnft_numtypes.h

    /** Step 2: Prepare calling fnft_kdvv **/ 

    // Location of the 1st and last sample of the continuous spectrum
    FNFT_REAL XI[2] = {-2.0, 2.0 };

    // Number of samples of the continuous spectrum
    FNFT_UINT M = 8;

    // Buffer for result: Samples of the continuous spectrum
    FNFT_COMPLEX contspec[M];

    // Maximum number of bound states we expect = size of the two arrays below
    FNFT_UINT K = D;

    // Buffer for result: Bound states
    FNFT_COMPLEX bound_states[K];

    // Buffer for result: Norming constants
    FNFT_COMPLEX normconsts[K];

    // Default options
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();

    // Uncomment the next line to compute residues instead of norming constants
    //opts.discspec_type = fnft_kdvv_dstype_RESIDUES;

    // See the header file fnft_kdvv.h for the options

    /** Step 3: Call fnft_nsev and check for errors **/

    int ret_code = fnft_kdvv(D, q, T, M, contspec, XI, &K, bound_states,
                             normconsts, &opts);

    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured!\n");
        return EXIT_FAILURE;
    }

    /** Step 4: Print the results **/

    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    FNFT_REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    printf("Continuous spectrum:\n");
    for (FNFT_UINT i=0; i<M; i++) {
        FNFT_REAL xi = XI[0] + i*eps_xi;
        printf("  continuous_spectrum(xi=%f) \t= %g + %gI\n",
            (double)xi,
            (double)FNFT_CREAL(contspec[i]),
            (double)FNFT_CIMAG(contspec[i])
        );
    }

    printf("Discrete spectrum:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  bound state at %g + %gI with norming constant %g + %gI\n",
            (double)FNFT_CREAL(bound_states[i]),
            (double)FNFT_CIMAG(bound_states[i]),
            (double)FNFT_CREAL(normconsts[i]),
            (double)FNFT_CIMAG(normconsts[i])
        );
    }

    return EXIT_SUCCESS;
}
