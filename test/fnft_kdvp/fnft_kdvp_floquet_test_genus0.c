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

static REAL scale(REAL x)
{
    const REAL absx = FABS(x);
    if (absx <= 1) {
        return x;
    } else {
        INT sign = x >= 0 ? 1 : -1;
        return sign * (1 + LOG(absx));
    }
}

static INT run_test()
{
    COMPLEX q[16];
    REAL T[2] = {0, 2};
    REAL E[2] = {-100, 100};

    REAL floq_det[100];
    REAL floq_det_exact[100];
    REAL al21[100];

    UINT D = sizeof(q)/sizeof(COMPLEX);
    UINT L = sizeof(floq_det)/sizeof(REAL);
    
    for (UINT i=0; i<D; i++)
        q[i] = 1;

    INT ret_code = fnft_kdvp_floquet(D, q, T, E, L, floq_det, al21, NULL/*opts*/);
    CHECK_RETCODE(ret_code, leave_fun);

    const REAL eps_E = (E[1] - E[0])/(L - 1);
    for (UINT i=0; i<L; i++) {
        const REAL Ei = E[0] + i*eps_E;
        floq_det_exact[i] = scale( CREAL(CCOS(-2*CSQRT(Ei - (-1)))) );
        floq_det[i] = scale( floq_det[i] );
    }
    REAL err = fnft__misc_rel_err_real(L, floq_det, floq_det_exact);
#ifdef DEBUG
    misc_print_buf(D, q, "q");
    misc_print_buf_real(K, floq_det, "floq_det");
    misc_print_buf_real(K, floq_det_exact, "floq_det_exact");
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


