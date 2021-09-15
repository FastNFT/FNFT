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
    INT kappa = -1;
    /* matlab code to compute results:
    %% Computing values for the test file
eps_t = 0.13;
kappa = -1;
D=8;
q = (0.1*cos(1:2*D)+0.08j*sin(0.28*(1:2*D)))*50;
r = -kappa*conj(q);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_3 = [z^-1 0 0; 0 z 0; 0 0 z];
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
          3.696021961216117 - 3.163512799264804*I,
  0.518535430451489 + 0.498626683962884*I,
  2.968428782523723 + 1.056506587397136*I,
 -1.264546299775730 - 2.218475127487348*I,
  0.682898305426271 + 0.466594908071217*I,
  2.074654567631780 + 2.627592817550341*I,
  0.337615808337937 - 0.361534720361276*I,
 -1.713647029153448 + 1.271388414710695*I,
  0.029801567670084 - 1.612513368120410*I,
  0.467369255888984 - 0.600384289047545*I,
  0.797406971584739 + 3.289999711357828*I,
  0.320755889946965 + 0.406236444083529*I,
 -2.258316091301153 + 0.648357813303253*I,
 -0.265406426558306 - 1.993993472354354*I,
 -0.383556764490316 - 0.191381975452150*I,
  0.459242940020080 - 3.370536549379147*I,
  0.184540988911419 + 0.354501396419772*I,
  2.345727486139551 + 0.491070214076699*I,
  0.838558959951546 + 1.645816000481939*I,
 -0.446603735335630 + 0.535965001027631*I,
  1.937938557501514 + 0.259333678472660*I,
  0.769345010979206 - 0.265110704648237*I,
 -1.517008756252161 + 1.516273002289303*I,
  0.455738709731700 + 1.076859829608434*I,
  0.652681714196180 + 0.700141300817626*I,
  2.384769317977041 + 1.749696030530174*I,
  0.078943928414987 + 0.114006247883830*I,
 -1.275506906701426 + 0.521996378119739*I,
 -0.269569604092390 + 1.605663173089122*I,
 -0.093426306581564 - 0.740009196088632*I,
  0.342034153362398 + 3.313605815147460*I,
  0.537436327679516 - 0.349859730280854*I,
  2.054399959518573 + 0.202584595978677*I,
 -0.555402926406675 + 1.739592717107572*I,
 -0.387452130341798 - 0.055674358889240*I,
 -2.894077762481048 + 0.099446895924724*I,
 -0.571076724870293 - 0.167996778690739*I,
 -0.600722070915458 + 0.774297827297973*I,
 -1.364515114884358 + 0.599230604603874*I,
  0.153734118132853 + 0.376283459374063*I,
 -1.826099212463049 - 0.613805136904898*I,
  0.829612720660173 - 0.304358844431496*I,
 -1.819734105619276 + 1.014594076633236*I,
 -0.633017098029578 + 1.302521585464421*I,
  0.270347001491963 + 0.574881744398677*I
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
            q1[i] = (0.1*COS(2*i+1) + 0.08*I*SIN(0.28*(2*i+1)))*50;
            q2[i] = (0.1*COS(2*(i+1)) + 0.08*I*SIN(0.28*(2*(i+1))))*50;
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
