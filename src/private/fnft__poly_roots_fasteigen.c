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
* Sander Wahls (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"

// Interface to the EISCOR root finding routine
extern INT z_poly_roots_modified_(INT *N, double complex const * const coeffs,
    double complex * const roots, double *threshold, INT *info);

// Fast computation of polynomial roots. See the header file for details.
INT poly_roots_fasteigen(const UINT deg,
    COMPLEX const * const p, COMPLEX * const roots)
{
    INT int_deg, info;
    double threshold = 1e8;
    // This threshold was used in the original routine. Set to INFINITY to
    // enforce QR. Set to 0 to enforce QZ.

	// Check inputs
	if (p == NULL)
		return E_INVALID_ARGUMENT(p);
	if (roots == NULL)
		return E_INVALID_ARGUMENT(roots);

    // Call Fortran root finding routine
    int_deg = (int)deg;
    z_poly_roots_modified_(&int_deg, p, roots, &threshold, &info);
    
    if (info == 0)
        return SUCCESS;
    else
        return E_SUBROUTINE(FNFT_EC_OTHER);
}
