/**
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
 * Sander Wahls (KIT) 2023.
 **/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft_kdvp.h"
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif

// This example run some tests for a soliton as given in the paper "The solitons
// of Zabusky and Kruskal revisited: Perspective in terms of the periodic
// spectral transform" by Osborne and Bergamasco, Physica D: Nonlin. Phen.
// 18(1-3),Jan. 1986, p. 26-46, https://doi.org/10.1016/0167-2789(86)90160-0

static INT run_test(const UINT D, const REAL err_bounds[3])
{
    INT ret_code = SUCCESS;

    const REAL g = 981;    // Gravity [cm/s^2]
    const REAL h = 5;      // Water depth [cm]
    const REAL u1 = 4;     // Soliton amplitude [cm]
    const REAL L = 100;    // Period [cm]
    
    // KdV parameters, see p. 27
    const REAL c0 = sqrt(g*h);
    const REAL al = 3*c0/(2*h);
    const REAL be = c0*h*h/6;
    const REAL lam = al/(6*be);

    // Generate the signal, Eq. 52
    COMPLEX q[D];
    const REAL L1 = SQRT(12*be/(al*u1));
    REAL X[2] = {0, L};
    const REAL eps_x = (X[1] - X[0])/D;
    for (UINT i=0; i<D; i++) {
        const REAL x = X[0] + i*eps_x;
        q[i] = misc_sech((x-L/2)/L1);
        q[i] *= lam*u1*q[i];
    }

    COMPLEX q_reconstructed[D];
    REAL E[2] = {-0.15, 0.001};

    // Prepare the calls to fnft_kdvp
    REAL main_spec[3*D];
    REAL aux_spec[D];
    const REAL main_spec_exact[6] = {-lam*u1/2, 1, -lam*u1/2, -1, 0, -1};

    UINT K = D;
    UINT M = D;

    fnft_kdvp_opts_t opts = fnft_kdvp_default_opts();
    opts.grid_spacing = 0.001;

    // Run some tests
ret_code = fnft_kdvp(D, q, X, E, &K, main_spec, &M, aux_spec, NULL/*sheet_indices*/, &opts);
CHECK_RETCODE(ret_code, leave_fun);

    if (K != 3) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    REAL err = misc_rel_err_real(K, main_spec, main_spec_exact);
#ifdef DEBUG
    printf("err(main_spec) = %g, err_bounds[0] = %g\n", err, err_bounds[0]);
#endif
   if (!(err <= err_bounds[0])) {
#ifdef DEBUG
        misc_print_buf_real(2*K, main_spec, "main_spec");
        misc_print_buf_real(2*K, main_spec_exact, "main_spec_exact");
#endif
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    opts.mainspec_type = kdvp_mstype_AMPLITUDES_MODULI_FREQS;
    K = D;
    M = D;
    ret_code = fnft_kdvp(D, q, X, E, &K, main_spec, &M, aux_spec, NULL/*sheet_indices*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    if (K != 1) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    const REAL A = main_spec[0];
    const REAL m = main_spec[1];

    err = FABS(A - lam*u1)/FABS(A);
#ifdef DEBUG
    printf("err(A) = %g, err_bounds[1] = %g\n", err, err_bounds[1]);
#endif
    if (!(err <= 0.003)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

#ifdef DEBUG
    printf("m = %g\n", m);
#endif
    if (!(m >= 0.99 && m<=1)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    // Reconstruct q[i] sample by sample using Eq. 13 (and run some tests)
    opts.mainspec_type = kdvp_mstype_EDGEPOINTS_AND_SIGNS;
    for (UINT i=0; i<D; i++) {

        K = D;
        M = D;
        ret_code = fnft_kdvp(D, q, X, E, &K, main_spec, &M, aux_spec, NULL/*sheet_indices*/, &opts);
        CHECK_RETCODE(ret_code, leave_fun);

        q_reconstructed[i] = 0;
        for (UINT j=0; j<K; j++)
            q_reconstructed[i] -= main_spec[2*j];
        for (UINT j=0; j<M; j++)
            q_reconstructed[i] += 2*aux_spec[j];

        const COMPLEX q0 = q[0];
        for (UINT j=0; j<D-1; j++)
            q[j] = q[j+1];
        q[D-1] = q0;
    }

    err = misc_rel_err(D, q_reconstructed, q);
#ifdef DEBUG
    printf("err(q_reconstructed) = %g, err_bounds[2] = %g\n\n", err, err_bounds[2]);
#endif
    if (!(err <= err_bounds[2])) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}

int main()
{
    REAL err_bounds[3] = {2e-4, 0.003, 0.07};
    if (run_test(128, err_bounds) != SUCCESS)
        return EXIT_FAILURE;

    err_bounds[0] = 1e-4;
    err_bounds[1] = 0.002;
    err_bounds[2] = 0.04;
    if (run_test(256, err_bounds) != SUCCESS)
        return EXIT_FAILURE;

    err_bounds[2] = 0.02;
    if (run_test(512, err_bounds) != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


