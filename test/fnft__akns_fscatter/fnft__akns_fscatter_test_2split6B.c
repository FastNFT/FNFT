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

static INT akns_fscatter_test_2split6B()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT6B;
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
//         % Dual of Eq. 22 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (81/40)*expm([0,q(n);r(n),0]*eps_t/6)*([1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t/3))^2*[1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t/6) ...
//             -(16/15)*expm([0,q(n);r(n),0]*eps_t/4)*[1,0;0,z]^3*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^3*expm([0,q(n);r(n),0]*eps_t/4) ...
//             +(1/24)*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^6*expm([0,q(n);r(n),0]*eps_t/2);
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
        -2.5364652588246582e+05 + -5.8469865460007498e+05*I,
        -2.7099648014637768e+02 + 1.4927489637477302e+03*I,
        2.7272793690706906e+06 + -5.9529212157737939e+06*I,
        -8.1963691658376362e+06 + -1.5153773332728422e+06*I,
        -1.0774716172053819e+04 + -1.3636165736406330e+04*I,
        -1.1310351652489512e+05 + -2.5628499486967054e+05*I,
        -2.1181952267983247e+03 + 6.5055781576754143e+03*I,
        -8.7029920880795777e+04 + -2.2806404945649626e+06*I,
        -1.6720450241802111e+07 + -5.5286143279438373e+06*I,
        -6.5491095407010135e+01 + 1.9870041510296114e+03*I,
        -5.9351529942121962e+05 + 1.7179924202664560e+05*I,
        7.4176295952072076e+02 + 1.9968902105092286e+03*I,
        -5.8677075229820963e+06 + -1.9456364068709419e+06*I,
        -3.8762137646675911e+06 + 7.3043660505557228e+06*I,
        -5.0129085258888617e+03 + 1.2582863347005448e+04*I,
        -2.6038691672786558e+05 + 7.7154271127687156e+04*I,
        2.1008666009300291e+03 + 9.3738428646725788e+03*I,
        -2.1323116082749330e+06 + 3.1678920051892614e+05*I,
        -1.0214905907203307e+07 + 1.4175887935035056e+07*I,
        1.3670949074215628e+03 + -7.2922553981908118e+02*I};
        
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
    if (akns_fscatter_test_2split6B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}