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
 * Lianne de Vries (TU Delft student) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"


static INT manakov_fscatter_test_2split3B()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_2SPLIT3B;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[8], q2[8];
    COMPLEX result[45];
    const REAL err_bnd = 100*EPSILON;
    UINT kappa = 1;
    /* matlab code to compute results:
%% Computing values for the test file
eps_t = 0.13;
kappa = 1;
D=8;
q = (0.41*cos(1:2*D)+0.59j*sin(0.28*(1:2*D)))*50;
r = -conj(q);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_3 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q(2*n-1), q(2*n); r(2*n-1), 0, 0; r(2*n), 0, 0];
        U = (9/8)*expm(B*eps_t/3)*eA_3^2*expm(B*eps_t*2/3)*eA_3 - ...
            (1/8)*expm(B*eps_t)*eA_3^3;
        S = U*S;
    end
    S = S*z^(D*3);
    result_exact(i) = S(1,1);
    result_exact(5+i) = S(1,2);
    result_exact(10+i) = S(1,3);
    result_exact(15+i) = S(2,1);
    result_exact(20+i) = S(2,2);
    result_exact(25+i) = S(2,3);
    result_exact(30+i) = S(3,1);
    result_exact(35+i) = S(3,2);
    result_exact(40+i) = S(3,3);
end
*/

    COMPLEX result_exact[45] = {    //(for D = 8)
        0.717079769494588 + 0.603849188423541*I,
 -0.669517227404168 - 1.126212701993989*I,
  0.107341438735568 - 1.422498517878053*I,
 -0.597737345404608 + 0.492495081969725*I,
 -0.077978335772884 - 0.888630188632644*I,
 -0.019426985252253 - 0.323809345097831*I,
 -0.949261504910700 - 0.497448887020238*I,
  0.062792417014885 + 0.176722176898955*I,
  1.224929978004748 - 0.301158036920894*I,
  1.079207028505957 - 0.722598764932581*I,
  0.117897313765894 + 0.045089749388624*I,
  0.577910492565164 + 0.199416759535317*I,
  0.653522152298728 - 0.430722704570378*I,
  0.053732162582479 + 0.121281245502546*I,
  0.258675766394777 - 0.438409426387958*I,
 -0.035861335498323 + 0.195304943612685*I,
 -0.576426953733082 - 0.276023805001874*I,
 -0.802743460977440 - 0.100946215413665*I,
  1.335571166981514 - 0.101984603650326*I,
  0.913042805984224 + 0.539247607532264*I,
 -0.354959894220306 + 0.563607365810547*I,
  0.006371392883519 - 0.354242388105807*I,
 -0.389046311115062 + 0.221849970356909*I,
  0.315364018852212 + 0.503048662661443*I,
 -0.010679233117618 - 0.377052706775228*I,
  0.718384875296164 - 0.029040063044732*I,
 -1.314042299506552 - 0.319358058388561*I,
  0.456818114245353 + 0.791482282050150*I,
  0.296332923344227 + 0.608355157214008*I,
 -0.931539880498571 - 0.838557136430182*I,
  0.043902236727246 + 0.282498316895673*I,
  0.680049142075385 - 0.971144231475341*I,
  0.001270817449440 + 0.081998130859728*I,
  0.581756869218355 - 0.476729319683220*I,
  1.019367412221558 - 0.959114325402739*I,
 -0.545709285372266 + 0.391563132622336*I,
 -0.390833633338943 + 1.176985480133516*I,
 -1.278575765063722 + 0.549706594921188*I,
  0.302449371117210 + 0.232873133362699*I,
 -0.743705443037592 - 0.359665672647954*I,
 -0.658806653082993 - 0.181991146610252*I,
 -1.292184403231214 + 0.056331639648755*I,
 -0.330903683696019 - 0.282628497671965*I,
 -1.477406394017100 - 0.126378941088928*I,
  0.626337813546028 - 0.614180907014008*I
    };
        
        i = manakov_fscatter_numel(D, manakov_discretization);
        if (i == 0) {
            ret_code = E_INVALID_ARGUMENT(manakov_discretization);
            goto leave_fun;
        }
        

        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<D; i++){
            q1[i] = (0.41*COS(2*i+1) + 0.59*I*SIN(0.28*(2*i+1)))*50;
            q2[i] = (0.41*COS(2*(i+1)) + 0.59*I*SIN(0.28*(2*(i+1))))*50;
        }
        
        // without normalization 
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, manakov_discretization);  // with kappa =1
    
        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }

        
        for (i=0; i<9; i++){
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
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }

#endif
        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        
        // with normalization
        W_ptr = &W;
        ret_code = manakov_fscatter(D, q1, q2, 1, eps_t, transfer_matrix, &deg, W_ptr, manakov_discretization); // with kappa = 1

        if (ret_code != SUCCESS){
            return E_SUBROUTINE(ret_code);
            goto leave_fun;
        }
        if (W == 0){
            return E_TEST_FAILED;
            goto leave_fun;
        }
        scl = POW(2.0, W);
        for (i=0; i<9*(deg+1); i++)
            transfer_matrix[i] *= scl;
        
        for (i=0; i<9; i++){
            for (j=0; j<nz; j++)
                result[i*nz+j] = z[j];
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
        
#ifdef DEBUG
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
#endif
        if (misc_rel_err(9*nz, result, result_exact) > 100*err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{

    if (manakov_fscatter_test_2split3B() != SUCCESS)
        return EXIT_FAILURE;
    
    return EXIT_SUCCESS;
}
