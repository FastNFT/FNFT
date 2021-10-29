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

#include "fnft__poly_chirpz.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT poly_chirpz_test()
{
    const UINT deg = 3;
    COMPLEX p[4] = {1+2*I, -3-0.5*I, 0.3, -0.4*I};
    COMPLEX result[6];
    COMPLEX result_exactM3[3] = { \
          -1.84195946931039 +      1.37868493949555*I, \
          -3.23124399963083 -     0.277190872296927*I, \
           -2.9642090953629 -      2.91478167338596*I };
    COMPLEX result_exactM6[6] = { \
          -1.84195946931039 +      1.37868493949555*I, \
          -3.23124399963083 -     0.277190872296927*I, \
           -2.9642090953629 -      2.91478167338596*I, \
         -0.560346883257616 -      4.87438215021699*I, \
           2.92615008874117 -      4.55045428561155*I, \
           5.42577886268045 -      1.63749297097422*I };
    /* Matlab code to generate the exact results:
    p = [1+2*1j, -3-0.5*1j, 0.3, -0.4*1j];
    A = 0.95; W = exp(0.3j); format long g;
    M = 3; Z = A * W.^-(0:(M-1));
    result_exact3 = polyval(p, 1./Z).'
    M = 6; Z = A * W.^-(0:(M-1));
    result_exact6 = polyval(p, 1./Z).'
    */
    COMPLEX A, W;
    UINT M;
    INT ret_code;
   
    A = 0.95;
    W = CEXP(0.3*I);

    M = 3;
    ret_code = poly_chirpz(deg, p, A, W, M, result);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(M, result, result_exactM3) > 100*EPSILON)
        return E_TEST_FAILED;

    M = 6;
    ret_code = poly_chirpz(deg, p, A, W, M, result);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(M, result, result_exactM6) > 100*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (poly_chirpz_test() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

