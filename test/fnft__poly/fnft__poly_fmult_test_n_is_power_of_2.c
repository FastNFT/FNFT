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

#include "fnft__poly_fmult.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static INT poly_fmult_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 2, n = 8, i;
    COMPLEX p[poly_fmult_numel(deg, n)]; // since n is already a power of two
    /* Matlab code to generate the "exact" result:
        i=0:23; p=sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
        r=p(1:3); for (k=4:3:24), r=conv(r,p(k:(k+2))); end
        format long g; r.'
    */
    COMPLEX result_exact[17] = {
           10690.1093069785 - I*7285.45565246217,
           22605.6811426431 - I*19510.065631389, 
           57389.5004107115 - I*20763.2575587549,
           53847.0938223071 - I*18963.5560762862,
           88648.7268642591 - I*2309.98190237675,
           66515.8698354361 - I*3878.11603829108,
           45160.5175520559 + I*12083.5337878339,
           64021.1527510823 + I*34487.976136233, 
          -26614.7840829275 + I*18912.0673174305,
           624.669240812309 + I*18294.2454730283,
           7755.19752834515 + I*18447.5313133858,
          -33812.2918248361 + I*1257.46809193287,
          -711.103319555697 - I*2023.66741077336,
           9781.54461527385 + I*407.75394285625, 
           39.3498537600159 + I*420.245708764338,
          -508.134987301475 + I*156.212629759735,
           4.42978955469815 - I*0.0813363246183392 };
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0]) - 1;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
    }
    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult(&deg, n, p, W_ptr);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (deg != deg_exact)
        return E_TEST_FAILED;
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<deg+1; i++)
            p[i] *= scl;
    }
    if (misc_rel_err(deg+1, p, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;
    return SUCCESS;
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult_test_n_is_power_of_2(0); // 1x1 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_n_is_power_of_2(1); // 1x1 without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
