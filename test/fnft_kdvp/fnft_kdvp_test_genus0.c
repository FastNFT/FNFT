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
 * Sander Wahls (KIT) 2023, 2025.
 **/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft_kdvp.h"
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif


/**
 * The Floquet determinant for a constant signal is provided in Eq. 3.5 of the
 * paper "Functiontheoretic properties of the discriminant of Hill's equation"
 * by Hochstadt, Math. Zeitschr. 82, 237-242 (1963),
 * https://doi.org/10.1007/BF01111426
 */

static INT run_test()
{
    COMPLEX q[16];
    REAL T[2] = {0, 2};
    REAL E[2] = {-100, 100};

    REAL main_spec[100];
    REAL aux_spec[100];

    UINT D = sizeof(q)/sizeof(COMPLEX);
    UINT K = sizeof(main_spec)/sizeof(REAL);
    UINT M = sizeof(aux_spec)/sizeof(REAL);
    
    fnft_kdvp_opts_t opts = fnft_kdvp_default_opts();

    for (UINT i=0; i<D; i++)
        q[i] = 1;

    opts.grid_spacing = 0.01;
    INT ret_code = fnft_kdvp(D, q, T, E, &K, main_spec, &M, aux_spec, NULL/*sheet_indices*/, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

#ifdef DEBUG
    misc_print_buf_real(K, main_spec, "main_spec");
    misc_print_buf_real(M, aux_spec, "aux_spec");
#endif
    if (K != 1) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (M != 0) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    REAL err = FABS(-1 - main_spec[0]);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 10*EPSILON) {
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


