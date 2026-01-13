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
* Sander Wahls (KIT) 2025-2026.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_kdvv.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

// This is the first example in the paper "Reliable computation of the eigenvalues
// of the discrete KdV spectrum" by Prins and Wahls, Appl. Math. Comput. 422,
// Nov. 2022, https://doi.org/10.1016/j.amc.2022.127361

INT run_test(const UINT D, const REAL err_bnd, const fnft_kdv_discretization_t discr)
{
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    opts.bound_state_localization = kdvv_bsloc_ACCOUNTING;
    opts.discretization = discr; 

    COMPLEX q[D];
    const REAL T[2] = {-10, 10};
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    for (UINT n=0; n<D; n++) {
        const REAL t = T[0] + n*eps_t;
        const COMPLEX tmp = misc_sech(2*t);
        q[n] = 99*tmp*tmp;
    }

    UINT K = D;
    COMPLEX bound_states[K];
    INT ret_code = fnft_kdvv(D, q, T, 0/*M*/, NULL/*contspec*/, NULL/*XI*/,
                             &K, bound_states, NULL/*normconsts_or_residues*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    COMPLEX bound_states_exact[5] = {1*I, 3*I, 5*I, 7*I, 9*I};
    const UINT K_exact = sizeof(bound_states_exact)/sizeof(COMPLEX);

    const REAL err = misc_hausdorff_dist(K, bound_states, K_exact, bound_states_exact);
#ifdef DEBUG
    misc_print_buf(K, bound_states, "bound_states");
    misc_print_buf(K_exact, bound_states_exact, "bound_states_exact");
    printf("err = %g\n", err);
#endif
    if (err > err_bnd) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    if (K != K_exact) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}



int main() {
    UINT D = 256;
    REAL err_bnd = 10*0.008;
    fnft_kdv_discretization_t discr = kdv_discretization_BO;

    INT ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    // check for quadratic convergence

    D *= 2;
    err_bnd /= 4;
    ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    D *= 2;
    err_bnd /= 4;
    ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    D = 256;
    err_bnd = 2e-5;
    discr = kdv_discretization_CF4_2;

    ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    // check for fourth order convergence

    D *= 2;
    err_bnd /= 16;
    ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    D *= 2;
    err_bnd /= 16;
    ret_code = run_test(D, err_bnd, discr);
    if (ret_code != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

