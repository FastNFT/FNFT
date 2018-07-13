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

static INT akns_fscatter_test_2split7A()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT7A;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q[8], r[8];
    COMPLEX result[20];
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
//         % Eq. 23 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (117649/46080)*[1,0;0,z]^15*(expm([0,q(n);r(n),0]*eps_t*2/7)*[1,0;0,z]^30)^3*expm([0,q(n);r(n),0]*eps_t/7) ...
//             -(15625/9216)*[1,0;0,z]^21*(expm([0,q(n);r(n),0]*eps_t*2/5)*[1,0;0,z]^42)^2*expm([0,q(n);r(n),0]*eps_t/5) ...
//             +(729/5120)*[1,0;0,z]^35*expm([0,q(n);r(n),0]*eps_t*2/3)*[1,0;0,z]^70*expm([0,q(n);r(n),0]*eps_t/3) ...
//             -(1/9216)*[1,0;0,z]^105*expm([0,q(n);r(n),0]*eps_t);
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
        -2.5364652588246911e+05 + -5.8469865460008197e+05*I,
        5.5684993140680064e+05 + 2.8265658008869342e+05*I,
        1.6392424056182402e+00 + 1.3470061979758112e+00*I,
        -5.6939399865803745e+04 + -2.1185324915824656e+05*I,
        -3.8380249663240457e+08 + -4.2844617857751119e+08*I,
        -1.1310351652489658e+05 + -2.5628499486967357e+05*I,
        2.1387292082607779e+06 + 1.2058715540076580e+06*I,
        -9.9578322635933603e-01 + -4.8035319037339548e+00*I,
        -7.0076924368406559e+04 + -6.9845422537580365e+04*I,
        -2.0946844733366787e+08 + -4.4544899434339178e+08*I,
        -5.9351529942122707e+05 + 1.7179924202664781e+05*I,
        -5.5186571565985063e+05 + 1.7082946209650242e+05*I,
        -1.7192897350216003e+00 + -4.2563092702919718e+00*I,
        -2.1232210113851487e+05 + -1.7381562687929122e+03*I,
        -1.6582700400487608e+08 + -4.8044348256275415e+08*I,
        -2.6038691672786878e+05 + 7.7154271127688175e+04*I,
        -2.0021710776243692e+06 + 2.9893454268236761e+05*I,
        -1.0565001435081031e+01 + 5.2428685944452162e+00*I,
        -8.3281131732680835e+04 + 4.7275648712395399e+04*I,
        -1.7754433935030937e+07 + -4.3378024805734760e+08*I};
        
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
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),250*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 250*EPSILON)
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
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),250*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 250*EPSILON)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    if (akns_fscatter_test_2split7A() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}