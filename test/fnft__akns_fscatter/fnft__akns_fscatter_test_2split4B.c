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

static INT akns_fscatter_test_2split4B()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT4B;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q[8], r[8];
    COMPLEX result[20];
    const REAL err_bnd = 100*EPSILON;
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
//         % Dual of Eq. 20 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (4/3)*expm([0,q(n);r(n),0]*eps_t/4)*[1,0;0,z]*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]*expm([0,q(n);r(n),0]*eps_t/4)...
//             -(1/3)*expm([0,q(n);r(n),0]*eps_t/2)*[1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t/2);
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
        -2.5364652588246705e+05 + -5.8469865460007975e+05*I,
        -1.9269176615336383e+05 + -5.6824314081065403e+05*I,
        3.5272061057158676e+03 + -5.1243841508914575e+03*I,
        -4.9167324137236783e+03 + 1.2146712453521726e+04*I,
        2.2285153754727158e+05 + -3.2081153373712313e+05*I,
        -1.1310351652489594e+05 + -2.5628499486967301e+05*I,
        -6.6724032849576033e+04 + -6.0562863077309029e+05*I,
        2.2224874874736379e+04 + 2.0791975431497322e+04*I,
        3.5102605995107024e+03 + -4.3974085422505304e+03*I,
        5.3562344034025460e+04 + -1.1753258831555711e+05*I,
        -5.9351529942122451e+05 + 1.7179924202664624e+05*I,
        -6.3412474766768841e+05 + -1.3593081265779503e+04*I,
        -2.6992785806193310e+03 + -7.0224208237957191e+03*I,
        1.0113908880406820e+04 + 1.0458629254715470e+02*I,
        -2.9980917206972407e+05 + -1.8742370452712273e+05*I,
        -2.6038691672786820e+05 + 7.7154271127687680e+04*I,
        -6.2581539865337207e+05 + -1.5224200159962551e+05*I,
        3.5937018679523862e+04 + -7.9498750502803323e+03*I,
        -4.1497239586833502e+03 + -1.2809074693008786e+03*I,
        -1.0860779405192050e+05 + -4.3276066382267825e+04*I};
        
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
    if (akns_fscatter_test_2split4B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
