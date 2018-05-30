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
 * Sander Wahls (TU Delft) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_specfact.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"
#ifdef DEBUG
#include <stdio.h>
#endif

static INT poly_specfact_test()
{
    const COMPLEX poly[4] = { 4.0, 3.0, 2.0, 1.0 };
    COMPLEX result[4];
    const COMPLEX result_exact[4] = { 1.0, 2.0, 3.0, 4.0 };
    const UINT deg = 3;
    const UINT oversampling_factor = 32;

    INT ret_code = poly_specfact(deg, poly, result, oversampling_factor);
    if (ret_code != SUCCESS)
        ret_code = E_SUBROUTINE(ret_code);

    REAL err = misc_rel_err(deg+1, result, result_exact);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 10*FNFT_EPSILON)
        ret_code = E_TEST_FAILED;

    return ret_code;
}

int main(void)
{
    INT ret_code = poly_specfact_test();
    if (ret_code != SUCCESS) {
        ret_code = E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
