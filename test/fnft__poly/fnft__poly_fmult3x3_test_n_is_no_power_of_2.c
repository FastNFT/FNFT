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
#include <stdio.h>

static INT poly_fmult3x3_test_n_is_power_of_2(INT normalize_flag)
{
    UINT deg = 1, n = 5;
    UINT i;
    INT W, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    /* Matlab code to generate the exact result:
    % Define the polynomial coefficients:
% p11 = [p11_1(1) p11_1(0) p11_2(1) ... p11_4(0)]
% where the subscript indicates the number of the matrix and (.) the
% coefficient. So p11_4(1) is the coefficient in front of z^1 of the upper
% left element of matrix 4. In C, we have the same, but put all pij in one
% long vector: p = [p11 p12 ... p33] and have pointers *pij indicating the
% memory address where each pij starts.
i = 0:7;
    p11 = sqrt(i+1).*(cos(i) + 1j*sin(-2*i));
    p12 = sqrt(i+1).*(cos(i+0.1) + 1j*sin(-2*i+0.1));
    p13 = sqrt(i+1).*(cos(i+0.2) + 1j*sin(-2*i+0.2));
    p21 = sqrt(i+1).*(cos(i+0.3) + 1j*sin(-2*i+0.3));
    p22 = sqrt(i+1).*(cos(i+0.4) + 1j*sin(-2*i+0.4));
    p23 = sqrt(i+1).*(cos(i+0.5) + 1j*sin(-2*i+0.5));
    p31 = sqrt(i+1).*(cos(i+0.6) + 1j*sin(-2*i+0.6));
    p32 = sqrt(i+1).*(cos(i+0.7) + 1j*sin(-2*i+0.7));
    p33 = sqrt(i+1).*(cos(i+0.8) + 1j*sin(-2*i+0.8));
    
% First multiply matrices 1 and 2:
    ind1 = 1:2; ind2 = 3:4;
    s11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    s12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    s13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    s21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    s22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    s23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    s31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    s32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    s33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    ind1 = 5:6; ind2 = 7:8;
    t11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    t12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    t13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    t21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    t22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    t23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    t31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    t32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    t33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    r11 = conv(s11,t11) + conv(s12,t21) + conv(s13,t31);
    r12 = conv(s11,t12) + conv(s12,t22) + conv(s13,t32);
    r13 = conv(s11,t13) + conv(s12,t23) + conv(s13,t33);
    r21 = conv(s21,t11) + conv(s22,t21) + conv(s23,t31);
    r22 = conv(s21,t12) + conv(s22,t22) + conv(s23,t32);
    r23 = conv(s21,t13) + conv(s22,t23) + conv(s23,t33);
    r31 = conv(s31,t11) + conv(s32,t21) + conv(s33,t31);
    r32 = conv(s31,t12) + conv(s32,t22) + conv(s33,t32);
    r33 = conv(s31,t13) + conv(s32,t23) + conv(s33,t33);

    format long g; result_exact = [r11 r12 r13 r21 r22 r23 r31 r32 r33]'
    */

    const UINT memsize = poly_fmult3x3_numel(deg, n);
    COMPLEX p[memsize], result[memsize];
    COMPLEX result_exact[54] = {
        -19.964853888108 -      I*989.577413598393,
          -4028.12279685482 -    I*  3715.67752882116,
          -8870.76111133752 +    I*  671.209686728029,
          -2306.15588988156 +    I*  3424.04889656632,
          -823.761315041304 -    I*   3784.0022643467,
          -4179.15486304363 +    I*  548.743992726927,
           195.961714933656 -    I*  1207.11620709836,
          -4142.66264753167 -    I*  4423.15361221925,
          -9407.35648728203 +    I*  513.501498474698,
          -2033.90304034706 +    I*  3600.67695074125,
          -921.546282469662 -    I*  4248.05228844997,
          -4273.04414595308 +    I*  625.012189927769,
           409.930299076136 -    I*  1412.59389447658,
           -4215.8103824166 -    I*  5086.43500682424,
          -9849.95666686423 +    I*  350.662572989634,
          -1741.32810395236 +    I*  3741.32823104993,
          -1010.12346406628 -    I*  4669.65717830643,
          -4324.23858423676 +    I*  695.035471928413,
            273.36696835462 -    I*  951.279429772063,
          -2985.97395042892 -    I*  3902.11832044015,
          -7309.52841104911 +    I*  1284.11184519067,
          -1313.45101088836 +    I*  4317.77394136328,
           -354.04669043496 -    I*  4099.52127253097,
          -3733.17279668002 +    I*  1616.95822840883,
           543.936607733713 -    I*  1095.29151276106,
          -3123.59806606296 -    I*  4483.87132743425,
          -7893.61732904301 +    I*  1167.64401776559,
          -944.536732180249 +    I*  4557.86184553525,
          -545.371795406077 -    I*   4545.2396403835,
          -3800.09541876213 +    I*  1712.62683597514,
            809.07141232987 -    I*  1228.35980500979,
          -3230.01222234515 -    I*  5020.82297429544,
           -8398.8358319681 +    I*  1039.50947728679,
           -566.18495466652 +    I*  4752.40910077545,
          -731.247725673446 -    I*  4945.54347620579,
          -3829.04874356451 +    I*  1791.18344231555,
           542.279733469735 -    I*   828.00648763331,
          -1677.09694399174 -    I*  3739.99450398284,
          -5095.35730741353 +    I*  1782.30811693009,
           -203.41946487957 +    I*  4825.80509939374,
           147.293870587954 -    I*  4048.84225484424,
          -2953.71752271539 +    I*  2540.73440125477,
           843.323263344778 -    I*  885.627689642094,
          -1825.51177221258 -    I*  4144.05817106477,
          -5674.76484397596 +    I*  1717.48437448669,
           229.202229204552 +    I*  5107.90651612515,
          -120.480870113078 -    I*  4436.41427210705,
          -2987.69548545195 +    I*   2647.2576273978,
           1135.94058593809 -    I*  934.399992525566,
          -1955.68669023943 -    I*  4506.71577874493,
          -6197.47200604524 +    I*  1635.50009589827,
           659.533810378649 +    I*  5338.97141939686,
          -387.051805785622 -    I*  4779.65910444655,
           -2991.8213824987 +    I*  2727.33033039489};
    const UINT deg_exact = sizeof(result_exact)/sizeof(result_exact[0])/9 - 1;

    for (i=0; i<(deg+1)*n; i++) {
        p[i] = SQRT(i+1.0)*(COS(i) + I*SIN(-2.0*i));
        p[i+n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.1) + I*SIN(-2.0*i+0.1));
        p[i+2*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.2) + I*SIN(-2.0*i+0.2));
        p[i+3*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.3) + I*SIN(-2.0*i+0.3));
        p[i+4*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.4) + I*SIN(-2.0*i+0.4));
        p[i+5*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.5) + I*SIN(-2.0*i+0.5));
        p[i+6*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.6) + I*SIN(-2.0*i+0.6));
        p[i+7*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.7) + I*SIN(-2.0*i+0.7));
        p[i+8*n*(deg+1)] = SQRT(i+1.0)*(COS(i+0.8) + I*SIN(-2.0*i+0.8));
    }

    if (normalize_flag)
        W_ptr = &W;
    ret_code = poly_fmult3x3(&deg, n, p, result, W_ptr);
        
    if (ret_code != SUCCESS){
        return E_SUBROUTINE(ret_code);
    }
    if (deg != deg_exact){
        return E_TEST_FAILED;
    }
    if (normalize_flag) {
        if (W == 0)
            return E_TEST_FAILED;
        scl = POW(2.0, W);
        for (i=0; i<9*(deg+1); i++)
            result[i] *= scl;
    }
    if (!(misc_rel_err(9*(deg+1), result, result_exact) <= 100*EPSILON)){
        REAL err;
        err = misc_rel_err(9*(deg+1), result, result_exact);
        return E_TEST_FAILED;
    }

    return SUCCESS; 
}

INT main(void)
{
    INT ret_code;

    ret_code = poly_fmult3x3_test_n_is_power_of_2(0); // 3x3 w/o normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    
    ret_code = poly_fmult3x3_test_n_is_power_of_2(1); // 3x3 with normalization
    if (ret_code != SUCCESS) {
        E_SUBROUTINE(ret_code);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS; // Replaced SUCCES FAILURE
}
