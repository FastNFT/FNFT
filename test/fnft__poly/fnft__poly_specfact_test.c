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
    /* Matlab code to generate this test case:

       roots_poly = [0.1+0.2j 0.2-0.3j 2+1j -2+2j exp(0.23j)];
       roots_specfact = roots_poly; idx = abs(roots_poly) < 1;
       roots_specfact(idx) = 1./conj(roots_specfact(idx));

       poly = 1; for r = roots_poly, poly = conv([1 r], poly); end
       result = 1; for r = roots_specfact, result = conv([1 r], result); end
       result = result/result(end); % make lowest coeff real
       result = result/abs(polyval(result,1))*abs(polyval(poly,1)); % match gain

       format long g; poly = poly.', result_exact = result.'
    */
    const COMPLEX poly[6] = {
        1 +                     0*I,
        1.27366639500537 +      3.12797752353519*I,
        -5.98903489975043 +      5.80202580257614*I,
        -7.76541973341761 +      2.99213552719788*I,
        -2.41536385774943 +      1.13047624544538*I,
        -0.509630949856206 -    0.0166221222670568*I };
    COMPLEX result[6];
    const COMPLEX result_exact[6] = {
        -0.0805797283830336 +   0.00262818829548505*I,
        -0.376517479322989 -     0.384616524437207*I,
        -0.341922071026159 -      1.56418586948815*I,
        2.39294448606361 -      3.94255780695439*I,
        9.0040568729159 -      3.65545082150954*I,
        6.32455532033676 +                     0*I };
    const UINT oversampling_factor = 32;
    INT ret_code = SUCCESS;

    const UINT deg = sizeof(poly)/sizeof(poly[0]) - 1;
    if (sizeof(result) != (deg+1)*sizeof(COMPLEX)
    || sizeof(result_exact) != (deg+1)*sizeof(COMPLEX)) {
        ret_code = E_ASSERTION_FAILED;
        goto leave_fun;
    }

    ret_code = poly_specfact(deg, poly, result, oversampling_factor, 0);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL err = misc_rel_err(deg+1, result, result_exact);
#ifdef DEBUG
    printf("err = %g\n", err);
#endif
    if (err > 0.0055) { // difficult case due to a zero on the unit circle
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
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
