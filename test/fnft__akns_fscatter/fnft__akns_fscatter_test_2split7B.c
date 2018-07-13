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

static INT akns_fscatter_test_2split7B()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT7B;
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
//         % Dual of Eq. 23 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = (117649/46080)*expm([0,q(n);r(n),0]*eps_t/7)*([1,0;0,z]^30*expm([0,q(n);r(n),0]*eps_t*2/7))^3*[1,0;0,z]^15 ...
//             -(15625/9216)*expm([0,q(n);r(n),0]*eps_t/5)*([1,0;0,z]^42*expm([0,q(n);r(n),0]*eps_t*2/5))^2*[1,0;0,z]^21 ...
//             +(729/5120)*expm([0,q(n);r(n),0]*eps_t/3)*[1,0;0,z]^70*expm([0,q(n);r(n),0]*eps_t*2/3)*[1,0;0,z]^35 ...
//             -(1/9216)*expm([0,q(n);r(n),0]*eps_t)*[1,0;0,z]^105;
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
        -2.5364652588246754e+05 + -5.8469865460008138e+05*I,
        7.5094380515345838e+04 + 2.1450732536718221e+06*I,
        1.0672637631622646e+01 + 7.8122470160806756e+00*I,
        -1.5501478895904984e+03 + -2.0501365898495686e+05*I,
        -3.4276265116125607e+08 + -4.6089893041513109e+08*I,
        -1.1310351652489572e+05 + -2.5628499486967319e+05*I,
        1.0356957614662722e+06 + -1.6735228923633788e+06*I,
        -5.5319562089589427e+00 + -1.2977887537436871e+01*I,
        -1.3040237982477980e+04 + -8.9029264674476930e+04*I,
        4.9224377882745564e+08 + 2.0857479839317107e+08*I,
        -5.9351529942122602e+05 + 1.7179924202664645e+05*I,
        1.5786492788206984e+06 + -9.7518345117490971e+05*I,
        5.2946172000061242e-02 + -1.0476370744956368e+01*I,
        -1.9571690075310669e+05 + -2.0955501768994513e+04*I,
        -2.6031246951505399e+08 + 4.1769822096180433e+08*I,
        -2.6038691672786820e+05 + 7.7154271127687418e+04*I,
        -1.7435356474756282e+06 + -1.4956288193373394e+05*I,
        -6.9095988546833320e+00 + 8.7244949762009618e+00*I,
        -8.6345345832495979e+04 + 2.6958880659339256e+03*I,
        1.2577631212938935e+07 + -4.5716557846451843e+08*I};
        
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
    if (akns_fscatter_test_2split7B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}