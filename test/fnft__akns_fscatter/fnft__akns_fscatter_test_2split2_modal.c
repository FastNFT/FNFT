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

static INT akns_fscatter_test_2split2_MODAL()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT2_MODAL;
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
//         % Eq. 25 in S. Wahls and V. Vaibhav, Fast Inverse Nonlinear Fourier Transforms for
//         % Continuous Spectra of Zakharov-Shabat Type, unpublished, https://arxiv.org/pdf/1607.01305v2.pdf
//         scl = 1/sqrt(1-eps_t*eps_t*q(n)*r(n));
//         U = [scl,scl*eps_t*q(n)*z;scl*eps_t*r(n),scl*z];
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
        2.4362083607185463e+00 + -3.9170503382218982e+00*I,
        -1.3775888434186105e+00 + 3.3565651922854451e+00*I,
        6.5584600528490977e-01 + -8.7176911462774553e-01*I,
        2.7515455576697201e-01 + 2.1257183073182084e-01*I,
        -3.4051346976221817e+00 + 1.3645300930588482e+00*I,
        3.6315325472324100e+00 + -1.3862968593776615e+00*I,
        -4.0552537495933514e+00 + -2.7334772090754500e-01*I,
        -1.3756179033993112e+00 + 1.1296566364231189e+00*I,
        -1.4416480197228005e+00 + 1.4332286545692863e-01*I,
        -3.0602833232433859e+00 + 4.7825952796141658e-02*I,
        -4.1382755348287787e+00 + -2.4363543653018587e+00*I,
        4.0273393046082075e+00 + 3.8384046079866480e-01*I,
        -1.4775976472447645e+00 + -1.9000091794310772e-01*I,
        -1.6349973835942599e-01 + 4.2350987018717462e-01*I,
        7.9869755208222437e-01 + 3.6323875551340112e+00*I,
        -1.4201678648804372e+00 + -3.5603238837569884e+00*I,
        9.1586007875177333e-01 + 4.1605872295350750e+00*I,
        2.0909580620036201e+00 + -2.8757367304714149e-02*I,
        -3.3166463294382809e-01 + 1.0995367474961899e+00*I,
        -4.4261059849256179e-01 + 2.7966366782242136e+00*I};
        
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
    if (akns_fscatter_test_2split2_MODAL() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}