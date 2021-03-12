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
 * Lianne de Vries (TU Delft) 2021.
 */

 // This file has only 1 timestep. Use to debug.

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"


static INT manakov_fscatter_test_2split3A_1D()
{
    UINT i, j, D = 1, deg, nz = 1;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT3A;
    const REAL eps_t = 0.13;
    COMPLEX z = CEXP(I*PI/4);
    COMPLEX q[16];
    COMPLEX result[9];
    const REAL err_bnd = 100*EPSILON;
    UINT kappa = 1;
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
//         % Eq. 17 in P. J. Prins and S. Wahls, Higher order exponential
//         % splittings for the fast non-linear Fourier transform of the KdV
//         % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
//         U = [1,0;0,z]*expm([0,q(n);r(n),0]*eps_t);
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

    COMPLEX result_exact[9] = {
-0.532021427865377 - 0.532021427865377*I,
  0.087949118638512 + 0.578425618436629*I,
 -0.727999719102720 + 0.214750495886504*I,
 -0.087949118638512 + 0.578425618436628*I,
 -0.641848873830694 + 0.641848873830694*I,
  0.095352160547856 + 0.072402928578374*I,
  0.727999719102720 + 0.214750495886504*I,
 -0.072402928578375 - 0.095352160547856*I,
 -0.597279335221230 + 0.597279335221230*I
    };
        
        i = manakov_fscatter_numel(D, akns_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto leave_fun;
        }

        printf("i = %d\n",i);
        
        transfer_matrix = malloc(i*sizeof(COMPLEX));        // replaced i by 7*9, because we do not use the trick to reduce degree, we need coefficients (and memory) up to z^6 coefficients
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<2*D; i++){
            q[i] = (0.41*COS(i+1) + 0.59*I*SIN(0.28*(i+1)))*50;
        }
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q, kappa, eps_t, transfer_matrix, &deg, NULL, akns_discretization);  // with kappa =1

        printf("transfer matrix = \n");
        for (j=0; j<7*9; j++){
            printf("%f + i%f\n", creal(transfer_matrix[j]), cimag(transfer_matrix[j]));
        }



        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }
        // TODO: fill in result
  /*      for (j=0; j<9; j++){        // filling in values of result

        ret_code = poly_eval(deg, transfer_matrix+(j-1), nz, &z);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
*/
        for (i=0; i<9; i++){    // replaced 4 by 9
                result[i*nz] = z;
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
        }
        }

        for (j=0; j<9; j++){
            printf("%f + i%f/n", creal(result[j]), cimag(result[j]));
        }

#ifdef DEBUG
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),err_bnd);
#endif
 /*       if (misc_rel_err(4*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
   */     
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{
    printf("Goes into main loop\n");

    if (manakov_fscatter_test_2split3A_1D() != SUCCESS)
        return EXIT_FAILURE;

    printf("Finishes main loop\n");
    printf("Compile test 1\n");
    
    return EXIT_SUCCESS;
}
