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

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft__errwarn.h"
#include "fnft__global_root_pole_finding_algorithm.h"
#include "fnft__misc.h"
#ifdef DEBUG
#include <stdio.h>
#endif

// We recreate the test case
// https://github.com/PioKow/GRPF/blob/master/example0_simple_rational_function/

INT grpf_fun(UINT K, COMPLEX * zs, COMPLEX * dest, void * opts)
{
    (void) opts;
    UINT k;
    for (k=0; k<K; k++) {
        const COMPLEX z = zs[k];
        if (z == I)
            return E_DIV_BY_ZERO;
        dest[k] = (z-1)*(z-I)*(z-I)*(z+1)*(z+1)*(z+1)/(z+I);
    }
    return SUCCESS;
}

INT test(const REAL Tol)
{
    const COMPLEX roots_true[3] = {1, I, -1};
    const UINT K_true = sizeof(roots_true)/sizeof(COMPLEX);
    COMPLEX roots[5];
    UINT K = sizeof(roots)/sizeof(COMPLEX);

    const UINT NodesMax = 500000;
    const REAL bounding_box[4] = {-2, 2, -2, 2};

    INT ret_code = global_root_pole_finding_algorithm(&K, roots, &grpf_fun, 
            NULL, NodesMax, bounding_box, Tol);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL err = misc_hausdorff_dist(K_true, roots_true, K, roots);
#ifdef DEBUG
    misc_print_buf(K, roots, "roots");
    printf("error = %g, error_bound = %g\n", err, Tol);
#endif
    if (!(err <= Tol)) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}

int main()
{
    INT ret_code = test(1e-2);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = test(1e-6);
    CHECK_RETCODE(ret_code, leave_fun);
 
    ret_code = test(1e-12);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
