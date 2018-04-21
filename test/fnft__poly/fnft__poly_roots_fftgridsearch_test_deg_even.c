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

#include "fnft__poly_roots_fftgridsearch.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT poly_roots_fftgridsearch_test(UINT M, REAL eb1, REAL eb2)
{
    const UINT deg = 6;
    COMPLEX p[7] = { 4-5*I, 3-4*I, 2-3*I, 2, 2+3*I, 3+4*I, 4+5*I };
	COMPLEX * roots;
    UINT nroots;
    COMPLEX roots_exact[6] = {
      3.992603696776205e-01 + 9.168375849652399e-01*I,
     -3.932716698145485e-01 + 9.194223152182454e-01*I,
      9.475398776668446e-01 - 3.196375763753407e-01*I,
     -9.853253052543915e-01 + 1.706869732151292e-01*I,
     -7.486910771535749e-01 - 6.629190531208312e-01*I,
      4.384974884943283e-17 - 1.000000000000000e+00*I };
    /* Matlab code for generating the "exact" roots:
    N = 4; c = (1:N) + 1j*(2:N+1); p = zeros(1,2*N-1);
    p(1:N) = conj(c(end:-1:1)); p(N:end) = p(N:end) + c;
    format long e; p.', roots(p)
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

	err = misc_hausdorff_dist(nroots, roots, deg, roots_exact);
    if (!(err <= eb1)) {
		ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    // Look for roots only in a small arc around the first root. (The angle of
    // the first root is approx 1.16.)
    PHI[0] = 1;
    PHI[1] = 1.5;
    nroots = M;
	ret_code = poly_roots_fftgridsearch(deg, p, &nroots, PHI, roots);
    if (ret_code != SUCCESS) {
		ret_code = E_SUBROUTINE(ret_code);
        goto release_mem;
    }

    err = misc_hausdorff_dist(nroots, roots, 1, roots_exact);
    if (!(err <= eb2)) {
		ret_code = E_TEST_FAILED;
        goto release_mem;
    }


release_mem:
    free(roots);
    return ret_code;
}

INT poly_roots_fftgridsearch_paraherm_test(UINT M, REAL eb1,
REAL eb2)
{
    const UINT deg = 6;
    COMPLEX p[7] = { 4-5*I, 3-4*I, 2-3*I, 2, 2+3*I, 3+4*I, 4+5*I };
	COMPLEX * roots;
    UINT nroots;
    COMPLEX roots_exact[6] = {
      3.992603696776205e-01 + 9.168375849652399e-01*I,
     -3.932716698145485e-01 + 9.194223152182454e-01*I,
      9.475398776668446e-01 - 3.196375763753407e-01*I,
     -9.853253052543915e-01 + 1.706869732151292e-01*I,
     -7.486910771535749e-01 - 6.629190531208312e-01*I,
      4.384974884943283e-17 - 1.000000000000000e+00*I };
    /* Matlab code for generating the "exact" roots:
    N = 4; c = (1:N) + 1j*(2:N+1); p = zeros(1,2*N-1);
    p(1:N) = conj(c(end:-1:1)); p(N:end) = p(N:end) + c;
    format long e; p.', roots(p)
    */
    REAL PHI[2] = { 0, 2.0*PI };
    REAL err;
    INT ret_code;   
 
    roots = malloc(M * sizeof(COMPLEX));
    if (roots == NULL)
        return E_NOMEM;

    nroots = M;
	ret_code = poly_roots_fftgridsearch_paraherm(deg, p, &nroots, PHI, roots);
    if (ret_code != SUCCESS) {
		ret_code = E_SUBROUTINE(ret_code);
        goto release_mem;
    }

	err = misc_hausdorff_dist(nroots, roots, deg, roots_exact);
    if (!(err <= eb1)) {
		ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    // Look for roots only in a small arc around the first root. (The angle of
    // the first root is approx 1.16.)
    PHI[0] = 1;
    PHI[1] = 1.5;
    nroots = M;
	ret_code = poly_roots_fftgridsearch_paraherm(deg, p, &nroots, PHI, roots);
    if (ret_code != SUCCESS) {
		ret_code = E_SUBROUTINE(ret_code);
        goto release_mem;
    }

    err = misc_hausdorff_dist(nroots, roots, 1, roots_exact);
    if (!(err <= eb2)) {
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
    REAL eb1 = 6.2e-4, eb2 = 3.6e-6;
    INT ret_code;

    ret_code = poly_roots_fftgridsearch_paraherm_test(M, eb1, eb2);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_roots_fftgridsearch_paraherm_test(2*M, eb1/4, eb2/4);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    eb1 = 0.0014;
    eb2 = 4.8e-6;
    ret_code = poly_roots_fftgridsearch_test(M, eb1, eb2);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    // Note: Second error reduces only 3.6-fold in this example ... 
    ret_code = poly_roots_fftgridsearch_test(2*M, eb1/4, eb2/3.6);
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
