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

static INT poly_fmult_test_n_is_no_power_of_2(INT normalize_flag)
{
    UINT deg = 2, n = 11, i;
    COMPLEX p[poly_fmult_numel(deg, n)];
    /* Matlab code to generate the "exact" result:
        deg=2; n=11; N=(deg+1)*n; i=0:(N-1);
        p=sqrt(i+1).*(cos(i) + 1j*sin(-2*i)); r=p(1:(deg+1));
        for (k=(deg+2):(deg+1):N), r=conv(r,p(k:(k+deg))); end
        format long g; r.'
    */
    COMPLEX result_exact[23] = {
          -319308.705234574 -      166784.182995777*I,
          -2743910.31810172 -      794705.736952064*I,
           -9075447.9140096 -      882018.812864758*I,
          -18191941.8165097 -      122049.064791405*I,
          -30115074.8590672 +      760936.799294374*I,
          -44814146.2207504 +      1917904.17878258*I,
          -59509427.3900375 +       1151280.7345219*I,
          -67870895.8105611 +       6758914.6822708*I,
          -65004270.8979727 +      2378651.90695769*I,
           -62350519.934101 -      4351181.14834246*I,
           -49802387.816833 +      6324801.48028804*I,
          -25632453.5014797 -      15474491.9397379*I,
          -19770246.1054521 -      11669704.6069006*I,
          -5178072.79922669 +      454130.609010228*I,
           7883730.23930349 -      23609117.6699819*I,
          -1012854.76439042 -      6787648.75896478*I,
           2010055.29299395 +      1306386.13720826*I,
           4346918.37729613 -      8535073.39464626*I,
          -2082296.12229051 -      649844.924066411*I,
          -1785325.93258564 +      2098730.17691877*I,
            50207.974399346 -      84587.5774434785*I,
           56394.1550202597 -      147664.217837454*I,
          -801.466849402119 +      1048.08952721767*I };
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
    if (!(misc_rel_err(deg+1, p, result_exact) <= 100*EPSILON))
        return E_TEST_FAILED;
    return SUCCESS;
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult_test_n_is_no_power_of_2(0); // without normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    ret_code = poly_fmult_test_n_is_no_power_of_2(1); // with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
