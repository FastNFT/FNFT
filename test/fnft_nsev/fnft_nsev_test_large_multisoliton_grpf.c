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
 * Sander Wahls (TU Delft) 2022.
 */

// Uncomment the line below to raise arithmetic exceptions e.g. on division by zero. Requires gcc.
//#define ENABLE_ARITHMETIC_EXCEPTIONS

#define FNFT_ENABLE_SHORT_NAMES
#ifdef ENABLE_ARITHMETIC_EXCEPTIONS
#define _GNU_SOURCE
#include <fenv.h>
#endif
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif
#include "fnft.h"
#include "fnft_nsev.h"
#include "fnft_nsev_inverse.h"
#include "fnft__misc.h"

/* This test utilizes Example 1 from the paper "Fast inverse nonlinear Fourier
 * transformation using exponential one-step methods: Darboux transformation"
 * by V. Vaibhav (Phys. Rev. E 96, 063302 (2017).
 */
static INT test(const REAL err_bound)
{
    const UINT D = 2048;
    REAL T[2] = {NAN};
    const UINT K_exact = 32;
    const UINT J = 4;
    const REAL th0 = PI/3;
    const REAL dth = (PI - 2*th0)/(J - 1);
    COMPLEX q[2048];
    COMPLEX bound_states_exact[32];
    COMPLEX normconsts_exact[32];
    COMPLEX bound_states[128];
    COMPLEX normconsts[128];
    UINT j, l;
    INT ret_code = SUCCESS;
    
    assert(sizeof(q)/sizeof(COMPLEX) == D);
    assert(sizeof(bound_states_exact)/sizeof(COMPLEX) == K_exact);
    assert(sizeof(normconsts_exact)/sizeof(COMPLEX) == K_exact);
    assert(sizeof(bound_states)/sizeof(COMPLEX) >= K_exact);
    assert(sizeof(normconsts) == sizeof(bound_states));
   
    for (l=1; l<=8; l++)
        for (j=1; j<=J; j++) {
            const REAL thj = th0 + (j-1)*dth;
            bound_states_exact[j+J*(l-1)-1] = l*CEXP(I*thj);
        }

    for (j=1; j<=8*J; j++) 
        normconsts_exact[j-1] = CEXP(I*PI*(j-1)/(8*J-1));

    REAL kappa = 0;
    REAL im_min = INFINITY;
    for (j=0; j<K_exact; j++) {
        const REAL im = CIMAG(bound_states_exact[j]);
        kappa += im;
        if (im < im_min)
            im_min = im;
    }
    kappa = 2*SQRT(kappa);
    for (j=0; j<K_exact; j++)
        bound_states_exact[j] /= kappa;
    const REAL L = 11*kappa/im_min;
    T[0] = -L;
    T[1] = L;
#ifdef DEBUG
    printf("T = [%g, %g]\n", T[0], T[1]);
#endif

    ret_code = fnft_nsev_inverse(0, NULL, NULL, K_exact, bound_states_exact,
        normconsts_exact, D, q, T, +1/*kappa*/, NULL);
    CHECK_RETCODE(ret_code, leave_fun);

    UINT K = sizeof(bound_states)/sizeof(COMPLEX);
    fnft_nsev_opts_t opts = fnft_nsev_default_opts();
    opts.bound_state_localization = nsev_bsloc_GRPF;
    opts.bound_state_filtering = fnft_nsev_bsfilt_MANUAL;
    opts.bounding_box[0] = -2;
    opts.bounding_box[1] = 2;
    opts.bounding_box[2] = 0;
    opts.bounding_box[3] = 5;
    opts.discretization = nse_discretization_ES4;
    opts.tol = err_bound;
    ret_code = fnft_nsev(D, q, T, 0/*M*/, NULL/*contspec*/, NULL/*XI*/, &K,
        bound_states, normconsts, +1/*kappa*/, &opts);

    const REAL err = misc_hausdorff_dist(K_exact, bound_states_exact, K, bound_states);
#ifdef DEBUG
    misc_print_buf(K_exact, bound_states_exact, "bound_states_exact");
    misc_print_buf(K, bound_states, "bound_states");
    printf("err = %g, err_bound = %g\n", err, err_bound);
#endif

    if (!(err <= err_bound)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}

int main() {
#ifdef ENABLE_ARITHMETIC_EXCEPTIONS
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    INT ret_code = test(0.01);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    ret_code = test(0.001);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
