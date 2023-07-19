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
 * Lianne de Vries (TU Delft student) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h> // for printf
#include "fnft_manakovv.h"
#include "fnft_manakov_discretization_t.h"
#include "fnft__misc.h"

int main()
{
    /** Step 1: Set up the signal **/

    // Number of time-domain samples. Increase to improve precision of results.
    FNFT_UINT D = 1024;

    // Contains the samples of q(t)
    FNFT_COMPLEX q1[D], q2[D];

    // Location of the 1st and last time-domain sample
    FNFT_REAL T[2] = { -2.0, 2.0 };

    // Define a simple rectangular signal
    // TODO: use rectangular signal
    for (UINT i=0; i<D; i++){
            q1[i] = 2;
            q2[i] = 0.65;
    }

    // The types FNFT_UINT, FNFT_INT, FNFT_REAL and FNFT_COMPLEX are defined
    // in the header fnft_numtypes.h


    /** Step 2: Prepare calling fnft_manakovv **/ 

    // Location of the 1st and last sample of the continuous spectrum
    FNFT_REAL XI[2] = {-1.75, 2.0 };

    // Number of samples of the continuous spectrum
    FNFT_UINT M = 8;

    // Buffer for result: Samples of the continuous spectrum
    FNFT_COMPLEX contspec[3*M];


    // Maximum number of bound states we expect = size of the two arrays below
    FNFT_UINT K = D;        

    // Buffer for result: Bound states
    FNFT_COMPLEX bound_states[K];

    // Buffer for result: Norming constants
    FNFT_COMPLEX normconsts[K];

    // Focusing nonlinear Schroedinger equation
    int kappa = +1;

    // Default options
    fnft_manakovv_opts_t opts = fnft_manakovv_default_opts();
    opts.discretization = manakov_discretization_2SPLIT3A;

    // Uncomment the next line to compute residues instead of norming constants
    //opts.discspec_type = fnft_nsev_dstype_RESIDUES;

    // See the header file fnft_nsev.h for other options

    /** Step 3: Call fnft_nsev and check for errors **/

    int ret_code = fnft_manakovv(D, q1, q2, T, M, contspec, XI, &K, bound_states,
        normconsts, kappa, &opts);
    if (ret_code != FNFT_SUCCESS) {
        printf("An error occured!\n");
        return EXIT_FAILURE;
    }

    /** Step 4: Print the results **/

    printf("Number of samples:\n  D = %u\n", (unsigned int)D);

    FNFT_REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    printf("Continuous spectrum (first entry):\n");
    for (FNFT_UINT i=0; i<M; i++) {
        FNFT_REAL xi = XI[0] + i*eps_xi;
        printf("  continuous_spectrum(xi=%f) \t= %g + %gI\n",
            (double)xi,
            (double)FNFT_CREAL(contspec[i]),
            (double)FNFT_CIMAG(contspec[i])
        );
    }

    printf("Continuous spectrum (second entry):\n");
    for (FNFT_UINT i=M; i<2*M; i++) {
        FNFT_REAL xi = XI[0] + (i-8)*eps_xi;
        printf("  continuous_spectrum(xi=%f) \t= %g + %gI\n",
            (double)xi,
            (double)FNFT_CREAL(contspec[i]),
            (double)FNFT_CIMAG(contspec[i])
        );
    }

/*    printf("Discrete spectrum:\n");
    for (FNFT_UINT i=0; i<K; i++) {
        printf("  bound state at %g + %gI with norming constant %g + %gI\n",
            (double)FNFT_CREAL(bound_states[i]),
            (double)FNFT_CIMAG(bound_states[i]),
            (double)FNFT_CREAL(normconsts[i]),
            (double)FNFT_CIMAG(normconsts[i])
        );
    }*/

    return EXIT_SUCCESS;
}
