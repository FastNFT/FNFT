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
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_kdvv.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

// This is the sixth example in the paper "Reliable computation of the eigenvalues
// of the discrete KdV spectrum" by Prins and Wahls, Appl. Math. Comput. 422,
// Nov. 2022, https://doi.org/10.1016/j.amc.2022.127361
//
// Since there is no closed-form for the bound states, we simply check that their
// number is correct and that grid search finds the same values.

int main()
{
    fnft_kdvv_opts_t opts = fnft_kdvv_default_opts();
    opts.bound_state_localization = kdvv_bsloc_ACCOUNTING;

    COMPLEX q[2] = {10, 10};
    const UINT D = sizeof(q)/sizeof(COMPLEX);
    const REAL T[2] = {-5, 5};

    UINT K = 30;
    COMPLEX bound_states[K];
    INT ret_code = fnft_kdvv(D, q, T, 0/*M*/, NULL/*contspec*/, NULL/*XI*/,
                             &K, bound_states, NULL/*normconsts_or_residues*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    COMPLEX bound_states_grid_search[100];
    UINT K_grid_search = sizeof(bound_states_grid_search)/sizeof(COMPLEX);

    opts.bound_state_localization = kdvv_bsloc_GRIDSEARCH_AND_REFINE;
    opts.grid_spacing = 0.001;
    ret_code = fnft_kdvv(D, q, T, 0/*M*/, NULL/*contspec*/, NULL/*XI*/,
                         &K_grid_search, bound_states_grid_search, NULL/*normconsts_or_residues*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    const REAL err = misc_hausdorff_dist(K, bound_states, K_grid_search, bound_states_grid_search);
#ifdef DEBUG
    misc_print_buf(K, bound_states, "bound_states");
    misc_print_buf(K_grid_search, bound_states_grid_search, "bound_states_grid_search");
    printf("err = %g\n", err);
#endif
    if (err > 100*EPSILON) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (K != 21) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

