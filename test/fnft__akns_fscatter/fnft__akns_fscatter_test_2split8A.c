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
 * Shrinivas Chimmalgi (TU Delft) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__akns_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static INT akns_fscatter_test_2split8A()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT8A;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q[8], r[8];
    COMPLEX result[20];
    const REAL err_bnd = 250*EPSILON;
    // The following MATLAB code has been used to compute the values below.
// eps_t = 0.13;
// kappa = 1; D=8; q = (0.41*cos(1:D)+0.59j*sin(0.28*(1:D)))*50;
// r = (0.33*sin(1:D) + 0.85j*cos(0.43*(1:D)))*25;
// zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
// result_exact = zeros(1,20);
// for i = 1:1:5
//     z = zlist(i);
//     S = eye(2);
//     for n=1:D
//         % Eq. 24 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (1024/315)*[1,0;0,z]^3*(expm([0,q(n);r(n),0]*eps_t/4)*[1,0;0,z]^6)^3*expm([0,q(n);r(n),0]*eps_t/4)*[1,0;0,z]^3 ...
//             -(729/280)*[1,0;0,z]^4*(expm([0,q(n);r(n),0]*eps_t/3)*[1,0;0,z]^8)^2*expm([0,q(n);r(n),0]*eps_t/3)*[1,0;0,z]^4 ...
//             +(16/45)*[1,0;0,z]^6*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^12*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^6 ...
//             -(1/360)*[1,0;0,z]^12*expm([0,q(n);r(n),0]*eps_t)*[1,0;0,z]^12;
//         S = U*S;
//     end
//     result_exact(i) = S(1,1);
//     result_exact(5+i) = S(1,2);
//     result_exact(10+i) = S(2,1);
//     result_exact(15+i) = S(2,2);
// end
// fprintf('result_exact \n')
// for i=1:1:20
//     fprintf('%.16e + %.16e*I,\n',real(result_exact(i)),imag(result_exact(i)))
// end
    COMPLEX result_exact[20] = {
        -2.5364652588246862e+05 + -5.8469865460008185e+05*I,
        6.5248496461951542e+09 + 7.0756507272697372e+09*I,
        -9.7338641425550270e+09 + -3.7347957102244296e+09*I,
        -2.2049427843689697e+10 + -2.6180438223164539e+10*I,
        6.5157962546570420e+07 + -2.4295972449897211e+07*I,
        -1.1310351652489761e+05 + -2.5628499486967403e+05*I,
        -2.8304770746896200e+09 + 7.4948296064692706e+07*I,
        -6.2491650296852016e+09 + -1.1872461292865297e+10*I,
        -9.7714982659138412e+09 + -5.1009751989001297e+10*I,
        -4.4595609617289320e+07 + -1.6704471323865645e+07*I,
        -5.9351529942122707e+05 + 1.7179924202664773e+05*I,
        -5.9623057491959763e+09 + 7.2454865762560492e+09*I,
        -1.8008589622840471e+09 + 1.0007851338901169e+10*I,
        -2.9111712281588634e+10 + 1.2828863288476505e+10*I,
        -1.0068598503289243e+07 + 6.2935929141706645e+07*I,
        -2.6038691672786948e+05 + 7.7154271127689281e+04*I,
        -4.5054308678994119e+08 + -2.7242498190607200e+09*I,
        -1.0232506807521317e+10 + 8.1562538994686470e+09*I,
        -4.8051047602493446e+10 + -4.6088618832487392e+09*I,
        3.3471303864763573e+07 + -2.8012785370724306e+07*I};
        
        i = akns_fscatter_numel(D, akns_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto leave_fun;
        }
        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<D; i++){
            q[i] = (0.41*COS(i+1) + 0.59*I*SIN(0.28*(i+1)))*50;
            r[i] = (0.33*SIN(i+1) + 0.85*I*COS(0.43*(i+1)))*25;
        }
        
        
        // without normalization
        ret_code = akns_fscatter(D, q, r, eps_t, transfer_matrix, &deg, NULL, akns_discretization);
        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }
        
        for (i=0; i<4; i++){
            for (j=0; j<nz; j++){
                result[i*nz+j] = z[j];
            }
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
#ifdef DEBUG
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact), err_bnd);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        // with normalization
        W_ptr = &W;
        ret_code = akns_fscatter(D, q, r, eps_t, transfer_matrix, &deg, W_ptr, akns_discretization);
        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }
        if (W == 0){
            return E_TEST_FAILED;
            goto leave_fun;
        }
        scl = POW(2.0, W);
        for (i=0; i<4*(deg+1); i++)
            transfer_matrix[i] *= scl;
        
        for (i=0; i<4; i++){
            for (j=0; j<nz; j++)
                result[i*nz+j] = z[j];
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
        
#ifdef DEBUG
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact), err_bnd);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    if (akns_fscatter_test_2split8A() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
