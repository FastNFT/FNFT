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

static INT akns_fscatter_test_2split3A()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT3A;
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
//         % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (9/8)*[1,0;0,z]*expm([0,q(n);r(n),0]*eps_t*2/3)*[1,0;0,z]^2*expm([0,q(n);r(n),0]*eps_t/3)...
//             -(1/8)*[1,0;0,z]^3*expm([0,q(n);r(n),0]*eps_t);
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
        -2.5364652588246533e+05 + -5.8469865460007777e+05*I,
        1.7636438284456890e+04 + 1.0608756090546973e+05*I,
        1.4294201885996270e+04 + -7.4864114204988437e+03*I,
        -1.6328182495286962e+03 + -2.2010327957618201e+03*I,
        -2.4640183481245738e+05 + 1.5288077292310388e+05*I,
        -1.1310351652489492e+05 + -2.5628499486967159e+05*I,
        -8.9579651615119496e+04 + 1.7332126889865543e+05*I,
        2.9492916905778534e+03 + -1.9456543842367112e+03*I,
        2.0111360642352120e+03 + -7.7357216465656256e+03*I,
        -2.2815521817782745e+04 + 5.6777725415047214e+04*I,
        -5.9351529942122241e+05 + 1.7179924202664488e+05*I,
        9.4933453818192793e+04 + 5.3061986298330201e+04*I,
        1.4015700305286677e+04 + 4.2936467305829392e+03*I,
        2.8156557144783246e+03 + 8.2386022079673103e+02*I,
        2.1259683641302871e+05 + 1.7087881129191990e+05*I,
        -2.6038691672786663e+05 + 7.7154271127686894e+04*I,
        8.0334651215980208e+04 + 1.8020674298012193e+05*I,
        3.1477767063484312e+03 + 6.2687394337069372e+02*I,
        3.2815730827537477e+03 + 7.9022907955611372e+03*I,
        5.7502108420400960e+04 + 2.5117150379776690e+03*I};
        
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
    if (akns_fscatter_test_2split3A() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
