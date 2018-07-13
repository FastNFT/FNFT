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

static INT akns_fscatter_test_2split6A()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT6A;
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
//         % Eq. 22 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (81/40)*[1,0;0,z]^2*(expm([0,q(n);r(n),0]*eps_t/3)*[1,0;0,z]^4)^2*expm([0,q(n);r(n),0]*eps_t/3)*[1,0;0,z]^2 ...
//             -(16/15)*[1,0;0,z]^3*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^6*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^3 ...
//             +(1/24)*[1,0;0,z]^6*expm([0,q(n);r(n),0]*eps_t)*[1,0;0,z]^6;
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
        -2.5364652588246652e+05 + -5.8469865460007684e+05*I,
        -3.1940532587536265e+04 + 7.8406021841001384e+03*I,
        1.4818388796144698e+05 + 4.4849902361451183e+05*I,
        7.5467806196995417e+05 + 1.3810431197497118e+06*I,
        4.4168769434060800e+02 + 5.1610689383516421e+02*I,
        -1.1310351652489563e+05 + -2.5628499486967092e+05*I,
        -1.5739550685316608e+04 + -2.9854102436549605e+04*I,
        2.6720200360064390e+06 + -2.6079321973812673e+06*I,
        6.6458909731201595e+05 + -1.2167496530698312e+06*I,
        2.7819281185956402e+02 + -3.8706553425849904e+02*I,
        -5.9351529942122172e+05 + 1.7179924202664616e+05*I,
        -2.1693362642996952e+04 + -3.1699171833128632e+04*I,
        1.5267038841447749e+05 + 4.6316131591804413e+05*I,
        8.8672955425362033e+05 + -1.1501341645130366e+06*I,
        -5.3369316795189525e+02 + -3.4138150111192851e+02*I,
        -2.6038691672786610e+05 + 7.7154271127687694e+04*I,
        2.6173176898736769e+04 + -2.9472056965449494e+04*I,
        9.0680009906209656e+05 + -1.8436217405413520e+06*I,
        -1.3345570589554703e+06 + -2.2007659382939557e+04*I,
        -1.4319334639165652e+02 + 4.2061378257570385e+02*I};
        
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
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),255*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 255*EPSILON)
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
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),255*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 255*EPSILON)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    if (akns_fscatter_test_2split6A() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
