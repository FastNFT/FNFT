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
#include <stdio.h>
#include "fnft__poly_roots_fftgridsearch.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT poly_roots_fftgridsearch_test_deg_odd(UINT M, REAL eb)
{
    const UINT deg = 3;
    COMPLEX p[4] = {
                          1 +                     0*I,
          -1.58378511059697 +      2.52620897565978*I,
           0.46009879991909 +      1.20456190643706*I,
           3.23268341994846 +      1.59679613801836*I };
    COMPLEX * roots;
    UINT nroots;
    const UINT nroots_exact = 2;
    COMPLEX roots_exact[2] = {
          0.504846104599857 +     0.863209366648873*I,
         -0.921060994002885 -      0.38941834230865*I };
    /* Matlab code for generating the "exact" roots:
    rr = [exp(0.4j) exp(-2.1j) -2+3j];
    p = 1; for r=rr, p=conv(p, [1 r]); end
    format long g; p = p.', rts = roots(p);
    result_exact=rts(abs(1-abs(rts))<100*eps)
    */
    REAL PHI[2] = { 0, 2.0*PI };
    REAL err;
    INT ret_code;   
 
    roots = malloc(M * sizeof(COMPLEX));
    if (roots == NULL)
        return E_NOMEM;

    nroots = M;
	ret_code = poly_roots_fftgridsearch(deg, p, &nroots, PHI, roots);
    if (ret_code != SUCCESS) {
		ret_code = E_SUBROUTINE(ret_code);
        goto release_mem;
    }
	err = misc_hausdorff_dist(nroots, roots, nroots_exact, roots_exact);
    if (!(err <= eb)) {
		ret_code = E_TEST_FAILED;
        goto release_mem;
    }

release_mem:
    free(roots);
    return ret_code;
}

int main()
{
    UINT M = 128;
    REAL eb = 0.0003;
    INT ret_code;

    ret_code = poly_roots_fftgridsearch_test_deg_odd(M, eb);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_roots_fftgridsearch_test_deg_odd(2*M, eb/4);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
