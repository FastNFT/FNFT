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


static INT manakov_fscatter_test_2split3A()
{
    UINT i, j, D = 8, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    manakov_discretization_t manakov_discretization = manakov_discretization_2SPLIT3A;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q1[8], q2[8];
    COMPLEX result[45];
    const REAL err_bnd = 100*EPSILON;
    INT kappa = 1;
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
        U = (9/8)*eA_3*expm(B*eps_t*2/3)*eA_3^2*expm(B*eps_t/3) - (1/8)*eA_3^3*expm(B*eps_t);
        S = U*S;
    end
    S=S*z^(D*3);
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
        0.717079769494587 + I*0.603849188423540,
 -0.669517227404166 - I*1.126212701993991,
 -1.205426270741966 - I*0.245452300373313,
  0.769324774846953 + I*0.458993691653882,
  0.366662225901886 - I*0.770892826233883,
 -0.019426985252253 - I*0.323809345097833,
  0.497448887020238 - I*0.949261504910699,
 -0.261436847501865 + I*0.690268022871104,
  0.405707233863645 + I*0.418097072974667,
 -0.773822833028319 - I*1.154624507704179,
  0.117897313765893 + I*0.045089749388625,
 -0.199416759535318 + I*0.577910492565164,
  0.149103076456717 + I*0.503216619744914,
 -0.046905113011129 - I*0.994090689472043,
 -0.286722955970471 - I*0.671651740246069,
 -0.035861335498322 + I*0.195304943612686,
 -0.276023805001875 + I*0.576426953733081,
 -0.361265241519421 + I*0.521300650126837,
 -0.564216734514218 - I*1.188404890187716,
  0.237984451838333 + I*1.156420490694242,
 -0.354959894220305 + I*0.563607365810547,
  0.006371392883521 - I*0.354242388105808,
 -0.035220527488110 - I*1.078493629107347,
  0.502424512067295 + I*0.515360830795772,
 -0.267789113789004 - I*0.145934634563106,
  0.718384875296165 - I*0.029040063044732,
 -1.314042299506553 - I*0.319358058388560,
  0.382835690520380 + I*0.630645653681684,
  0.627309252316094 - I*0.153204191181139,
 -0.929044002718357 - I*0.620229959736022,
  0.043902236727244 + I*0.282498316895673,
 -0.971144231475345 - I*0.680049142075383,
 -0.764328966546298 + I*0.398253528370027,
  0.142408341594341 - I*0.645985906800806,
  0.945456896799084 + I*0.942301939995414,
 -0.545709285372266 + I*0.391563132622335,
 -0.390833633338943 + I*1.176985480133516,
 -0.440343978096879 - I*0.450650866020793,
 -0.970994183180071 - I*0.003891979641344,
 -0.739198781192066 - I*0.057247192380030,
 -0.658806653082993 - I*0.181991146610252,
 -1.292184403231215 + I*0.056331639648754,
 -0.810528680854822 - I*0.511282172237863,
  0.049872508514115 - I*1.244805644515133,
  0.877255914293173 + I*0.195774534543759,
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
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, NULL, manakov_discretization);
    
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

#define DEBUG
#ifdef DEBUG
        printf("error without normalization = %2.1e < %2.1e\n",misc_rel_err(9*nz, result, result_exact),err_bnd);
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }
#endif  // DEBUG

        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        // with normalization
        W_ptr = &W;
        ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, W_ptr, manakov_discretization);

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
        printf("result and exact result\n");
        for (UINT j = 0; j<45; j++){
                printf("%f + i%f,    %f + i%f\n", creal(result[j]), cimag(result[j]), creal(result_exact[j]), cimag(result_exact[j]));
            }
#endif
        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{

    if (manakov_fscatter_test_2split3A() != SUCCESS)
        return EXIT_FAILURE;
        
    return EXIT_SUCCESS;
}
