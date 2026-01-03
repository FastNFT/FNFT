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
 * Sander Wahls (KIT) 2025.
 */

// This example demonstrates the use of the functions fnft_kdvp and 
// fnft_kdvp_ampmodfreq, which implement the nonlinear Fourier transform
// for the Korteweg-de Vries equation with periodic boundary conditions.
// The signal is the periodic continuation of a single soliton. Since
// the period is quite large, the amplitude of the soliton is recovered
// very well. This example was originally considered in the paper
// "The solitons of Zabusky and Kruskal revisited: Perspective in terms
// of the periodic spectral transform" by Osborne and Bergamasco,
// Physica D: Nonlin. Phen. 18(1-3),Jan. 1986, p. 26-46,
// https://doi.org/10.1016/0167-2789(86)90160-0

#include <stdio.h> // for printf
#include "fnft_kdvp.h"
#include "fnft__misc.h" // for sech

int main()
{
    /** Step 1: Set up the signal **/

    // Number of time-domain samples. Increase to improve precision of results.
    FNFT_UINT D = 256;

    // Contains the samples of q(t)
    FNFT_COMPLEX q[D];

    // Location of the 1st and last time-domain sample
    FNFT_REAL T[2] = { 0.0, 100.0 };

    // Define the soliton
    const FNFT_REAL g = 981;    // Gravity [cm/s^2]
    const FNFT_REAL h = 5;      // Water depth [cm]
    const FNFT_REAL u1 = 4;
    const FNFT_REAL P = 100;

    const FNFT_REAL c0 = FNFT_SQRT(g*h);
    const FNFT_REAL al = 3*c0/(2*h);
    const FNFT_REAL be = c0*h*h/6;
    const FNFT_REAL lam = al/(6*be);
    const FNFT_REAL P1 = FNFT_SQRT(12*be/(al*u1));

    const FNFT_REAL eps_x = (T[1]-T[0])/(D-1);
    for (FNFT_UINT i=0; i<D; i++) {
        const FNFT_REAL x = T[0] + i*eps_x;
        const FNFT_REAL tmp = fnft__misc_sech((x-P/2)/P1);
        q[i] = lam*u1*tmp*tmp;
    }

    // The types FNFT_UINT, FNFT_INT, FNFT_REAL and FNFT_COMPLEX are defined
    // in the header fnft_numtypes.h

    /** Step 2: Prepare calling fnft_kdvp **/ 

    // Interval in which the main and auxiliary spectrum will be looked for
    FNFT_REAL E[2] = {-0.03, 0.015 };

    // Maximum number expected in the main spectrum points
    FNFT_UINT K = D;

    // Buffer for result: Main spectrum
    FNFT_REAL main_spec[2*K];

    // Maximum number of points expected in the auxiliary spectrum
    FNFT_UINT M = D;

    // Buffer for result: Auxiliary spectrum
    FNFT_REAL aux_spec[M];

    // Buffer for result: Sheet indices
    FNFT_REAL sheet_indices[M];

    // Default options
    fnft_kdvp_opts_t opts = fnft_kdvp_default_opts();

    // Set grid spacing used during the localization of main and auxiliary spectra
    opts.grid_spacing = 0.001;

    // See the header file fnft_kdvp.h for more options

    /** Step 3: Call fnft_nsev and check for errors **/

    int ret_code = fnft_kdvp(D, q, T, E, &K, main_spec, &M, aux_spec,
                             sheet_indices, &opts);

    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured in fnft_kdvp!\n");
        return EXIT_FAILURE;
    }

    /** Step 4: Print the results so far **/

    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    printf("Main spectrum:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  E = %g\n", (double)main_spec[2*i]);
    }

    printf("Auxiliary spectrum and sheet indices:\n");
    for (FNFT_UINT i=0; i<M; i++) {
        printf("  E = %g,\tsheet index = %g\n",
            (double)aux_spec[i],
            (double)sheet_indices[i]
        );
    }

    /** Step 5: Compute ampltiudes, moduli and frequencies **/

    // Buffer for result: Amplitudes, moduli and wave numbers
    // (instead of frequencies since we consider a space series)
    FNFT_REAL ampmodfreq[3*K];

    ret_code = fnft_kdvp_ampmodfreq(&K, main_spec, ampmodfreq);

    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured in fnft_kdvp_ampmodfreq!\n");
        return EXIT_FAILURE;
    }

    /** Step 6: Print the next results **/

    printf("Amplitudes, moduli and wave numbers of the hyperelliptic modes:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  A = %g, m = %g, k = %g\n",
                (double)ampmodfreq[3*i]/lam,
                (double)ampmodfreq[3*i+1],
                (double)ampmodfreq[3*i+2]
        );
    }

    return EXIT_SUCCESS;
}
