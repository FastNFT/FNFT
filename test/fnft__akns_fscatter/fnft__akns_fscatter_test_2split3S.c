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

static INT akns_fscatter_test_2split3S()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT3S;
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
//         % Eq. 14.4 in S. Brustein and A. Mirin, Third Order Difference
//         % Methods for Hyperbolic Equations, Journal of Computational
//         % Physics, 5, 547-571, 1970.
//         U = (2/3)*[1,0;0,z]*expm([0,q(n);r(n),0]*eps_t)*[1,0;0,z]...
//             +(2/3)*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t/2)...
//             -(1/6)*[1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t)...
//             -(1/6)*expm([0,q(n);r(n),0]*eps_t)*[1,0;0,z]^2;
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
        -2.5364652588246707e+05 + -5.8469865460007964e+05*I,
        -7.9981350034971489e+04 + -1.7337645683718630e+05*I,
        4.2935389670501365e+02 + 2.9500273425209325e+03*I,
        4.8797621855183916e+03 + 4.0336698076848978e+02*I,
        2.1855133883050995e+05 + -1.8636825746137061e+05*I,
        -1.1310351652489601e+05 + -2.5628499486967269e+05*I,
        2.3735323609251605e+04 + -2.9808549682597150e+04*I,
        -8.8762882538016333e+02 + 5.3405784115865899e+01*I,
        -2.8712010132545101e+03 + -1.0483687446436270e+03*I,
        -1.5408223888561050e+04 + -5.9244514897569890e+04*I,
        -5.9351529942122439e+05 + 1.7179924202664633e+05*I,
        -1.9815596709593511e+05 + 3.5162462471491817e+03*I,
        -1.9100469862363084e+03 + 1.7393343960369830e+03*I,
        -4.2975806638588610e+03 + 2.8139163715884179e+03*I,
        -2.0162341054384012e+05 + -1.8027830730184040e+05*I,
        -2.6038691672786785e+05 + 7.7154271127687796e+04*I,
        -1.8395972246269564e+04 + -3.5012665249175712e+04*I,
        -4.6894773669419806e+02 + -6.1150182595382080e+02*I,
        3.0519535368991951e+03 + -9.8419663437310237e+02*I,
        -5.3323981102272301e+04 + 2.1893810345041584e+04*I};
        
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
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),100*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 100*EPSILON)
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
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),100*EPSILON);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > 100*EPSILON)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    if (akns_fscatter_test_2split3S() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
