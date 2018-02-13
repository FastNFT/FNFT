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

#include "fnft__poly_roots_fasteigen.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT poly_roots_fasteigen_test()
{
    const UINT deg = 3;
    COMPLEX p[4] = { 1.0-2.0*I, 0.3+0.4*I, -2.0-2.0*I, -3.0+4.0*I };
	COMPLEX roots[3];
    COMPLEX roots_exact[3] = { \
         -0.767344914566607 -      1.47758771489852*I, \
           1.26733490516498 +     0.334189743482641*I, \
         -0.399989990598367 +     0.943397971415882*I };
    INT ret_code;   
 
	ret_code = poly_roots_fasteigen(deg, p, roots);
    if (ret_code != SUCCESS)
		return E_SUBROUTINE(ret_code);
	if ( misc_hausdorff_dist(deg, roots, deg, roots_exact)
    > 100*EPSILON )
		return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if ( poly_roots_fasteigen_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
