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
* Sander Wahls (TU Delft) 2017.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include "fnft__akns_discretization_t.h"

INT poly_eval_test()
{
    const UINT deg = 3;
    const UINT nz = 2;
    COMPLEX p[4] = { 1.0-2.0*I, 0.3+0.4*I, -2.0-2.0*I, -3.0+4.0*I };
    COMPLEX z[2] = { 0.5+0.25*I, -1.0+2.0*I };
    COMPLEX result_exact[2];
    UINT i, j;
    INT ret_code;

    for (i=0; i<nz; i++) {
        result_exact[i] = 0.0;
        for (j=0; j<=deg; j++)
            result_exact[i] += p[j] * CPOW(z[i], deg-j);
    }

    ret_code = poly_eval(deg, p, nz, z);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if ( misc_rel_err(nz, z, result_exact) > 100*DBL_EPSILON )
        return E_TEST_FAILED;
    
    return SUCCESS;
}

INT poly_evalderiv_test()
{
    const UINT deg = 3;
    const UINT nz = 2;
    COMPLEX p[4] = { 1.0-2.0*I, 0.3+0.4*I, -2.0-2.0*I, -3.0+4.0*I };
    COMPLEX z[2] = { 0.5+0.25*I, -1.0+2.0*I };
    COMPLEX deriv[2];
    COMPLEX result_exact[2], deriv_exact[2];
    UINT i, j;
    INT ret_code;

    for (i=0; i<nz; i++) {
        result_exact[i] = 0.0;
        for (j=0; j<=deg; j++)
            result_exact[i] += p[j] * CPOW(z[i], deg-j);

        deriv_exact[i] = 0.0;
        for (j=1; j<=deg; j++)
            deriv_exact[i] += (deg - (j-1))*p[j-1] * CPOW(z[i], deg-j);
    }

    ret_code = poly_evalderiv(deg, p, nz, z, deriv);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if ( misc_rel_err(nz, z, result_exact) > 100*DBL_EPSILON )
        return E_TEST_FAILED;
    if ( misc_rel_err(nz, deriv, deriv_exact) > 100*DBL_EPSILON )
        return E_TEST_FAILED;
    return SUCCESS;
}

int main()
{
    if ( poly_eval_test() != SUCCESS )
        return EXIT_FAILURE;
    if ( poly_evalderiv_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
