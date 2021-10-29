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

#include "fnft.h"
#include "fnft__poly_eval.h"
#include "fnft__errwarn.h"

INT poly_eval(const UINT deg, COMPLEX const * const p, const UINT nz, \
    COMPLEX * const z)
{
    UINT i, j;
    COMPLEX tmp;    

    // Check inputs
    if (p == NULL)
        return E_INVALID_ARGUMENT(p);
    if (z == NULL)
        return E_INVALID_ARGUMENT(z);

    for (i=0; i<nz; i++) {
        // Evaluate p(z) at z=z[i]
        if (FNFT_CABS(z[i]) <= 1.0) { // Horner's method
            tmp = p[0];
            for (j=1; j<=deg; j++)
                tmp = p[j] + tmp*z[i];
            z[i] = tmp;
        } else { // Apply Horner's method to reversed polynomial for stability
            tmp = p[deg];
            for (j=1; j<=deg; j++)
                tmp = p[deg - j] + tmp/z[i];
            z[i] = tmp * FNFT_CPOW(z[i], deg);
        }
    }

    return SUCCESS;
}

INT poly_evalderiv(const UINT deg, COMPLEX const * const p, \
    const UINT nz, COMPLEX * const z, \
    COMPLEX * const deriv)
{
    UINT i, j;
    COMPLEX tmp, tmp2;

    // Check inputs
    if (p == NULL)
        return E_INVALID_ARGUMENT(p);
    if (z == NULL)
        return E_INVALID_ARGUMENT(z);

    for (i=0; i<nz; i++) {
        // Evaluate p(z) and p'(z) at z=z[i]
        if (CABS(z[i]) <= 1.0) { // Horner's method
            tmp = p[0];
            tmp2 = 0.0;
            for (j=1; j<=deg; j++) {
                tmp2 = tmp + tmp2*z[i];
                tmp = p[j] + tmp*z[i];
            }
            z[i] = tmp;
            deriv[i] = tmp2;
        } else { // Apply Horner's method to reversed polynomial for stabiltiy
            tmp = p[deg];
            tmp2 = 0.0;
            for (j=1; j<=deg; j++) {
                tmp2 = tmp + tmp2/z[i];
                tmp = p[deg - j] + tmp/z[i];
            }
            deriv[i] = CPOW(z[i], deg - 1)*(deg*tmp - tmp2/z[i]);
            z[i] = tmp * CPOW(z[i], deg);
        }
    }
    return SUCCESS;
}
