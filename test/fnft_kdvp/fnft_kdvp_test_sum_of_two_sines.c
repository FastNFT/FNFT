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
 * Sander Wahls (KIT) 2025.
 **/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft_kdvp.h"
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif

// Verifies that the amplitudes and frequencies of a sum of two sines
// with (relatively) low amplitudes are detected correctly.

static INT run_test()
{
    COMPLEX q[256];
    REAL T[2] = {0, 2};
    REAL E[2] = {-100, 800};

    REAL main_spec[100];
    REAL aux_spec[100];

    UINT D = sizeof(q)/sizeof(COMPLEX);
    UINT K = sizeof(main_spec)/sizeof(REAL);
    UINT M = sizeof(aux_spec)/sizeof(REAL);

    const REAL A1 = 2;
    const REAL A2 = 1;
    const REAL f1 = 4;
    const REAL f2 = 8;

    fnft_kdvp_opts_t opts = fnft_kdvp_default_opts();

    const REAL dt = T[1]/D;
    for (UINT i=0; i<D; i++) {
        const REAL t = i*dt;
        q[i] = A1*SIN(2*FNFT_PI*f1*t) - A2*SIN(2*FNFT_PI*f2*t);
    }

    opts.grid_spacing = 0.01;
    INT ret_code = fnft_kdvp(D, q, T, E, &K, main_spec, &M, aux_spec, NULL/*sheet_indices*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
#ifdef DEBUG
    misc_print_buf_real(2*K, main_spec, "main_spec");
    misc_print_buf_real(M, aux_spec, "aux_spec");
#endif

    REAL ampmodfreq[100];
    ret_code = fnft_kdvp_ampmodfreq(&K, main_spec, ampmodfreq);
    CHECK_RETCODE(ret_code, leave_fun);
#ifdef DEBUG
    misc_print_buf_real(3*K, ampmodfreq, "ampmodfreq");
#endif
    if (K != 2) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    REAL err = FABS(FABS(A1) - ampmodfreq[0]);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 0.004) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    err = FABS(f1 - ampmodfreq[2]);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 2e-5) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    err = FABS(FABS(A2) - ampmodfreq[3]);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 0.007) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    err = FABS(f2 - ampmodfreq[5]);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 7e-6) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}

int main()
{
    if (run_test() == SUCCESS)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}


