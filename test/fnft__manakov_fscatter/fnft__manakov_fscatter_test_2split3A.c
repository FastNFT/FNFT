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

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"


static INT manakov_fscatter_test_2split3A()
{
    UINT i, j, D = 2, deg, nz = 5;
    INT W = 0, *W_ptr = NULL;
    REAL scl;
    INT ret_code;
    COMPLEX *transfer_matrix = NULL;
    akns_discretization_t akns_discretization = akns_discretization_2SPLIT3A;
    const REAL eps_t = 0.13;
    COMPLEX z[5] = {1.0+0.0*I, CEXP(I*PI/4), CEXP(I*9*PI/14), CEXP(I*4*PI/3), CEXP(I*-PI/5)};
    COMPLEX q[16];
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
        U = (9/8)*eA_3*expm(B*eps_t*2/3)*eA_3^2*expm(B*eps_t/3) - (1/8)*eA_3^3*expm(B*eps_t);
        S = U*S;
    end
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
 -0.669517227404162 - I*1.126212701993993,
  0.507530878324161 - I*1.120585442302792,
  0.769324774846946 + I*0.458993691653892,
  0.156483462408368 + I*0.839184046222803,
 -0.019426985252253 - I*0.323809345097833,
  0.497448887020241 - I*0.949261504910698,
 -0.614786390906891 - I*0.408481165031849,
  0.405707233863639 + I*0.418097072974671,
  1.304707080119249 + I*0.479269199713260,
  0.117897313765893 + I*0.045089749388625,
 -0.199416759535320 + I*0.577910492565163,
 -0.523778484292012 + I*0.033388518821856,
 -0.046905113011117 - I*0.994090689472042,
  0.626750731650744 + I*0.374846147147430,
 -0.035861335498322 + I*0.195304943612686,
 -0.276023805001876 + I*0.576426953733080,
 -0.427841475498329 - I*0.468207875196650,
 -0.564216734514202 - I*1.188404890187721,
 -0.872260375813117 - I*0.795680078549576,
 -0.354959894220305 + I*0.563607365810547,
  0.006371392883522 - I*0.354242388105808,
  1.059290846798113 + I*0.205649934284978,
  0.502424512067288 + I*0.515360830795777,
  0.302424169958795 - I*0.039338892380196,
  0.718384875296165 - I*0.029040063044732,
 -1.314042299506551 - I*0.319358058388564,
 -0.700023005876833 + I*0.232905340614981,
  0.627309252316095 - I*0.153204191181131,
  1.116174410084062 - I*0.044301785681687,
  0.043902236727244 + I*0.282498316895673,
 -0.971144231475342 - I*0.680049142075386,
 -0.218189285447087 - I*0.833785390659370,
  0.142408341594348 - I*0.645985906800804,
 -1.318761880495349 - I*0.206612662672064,
 -0.545709285372266 + I*0.391563132622335,
 -0.390833633338946 + I*1.176985480133515,
  0.537337861200740 - I*0.329024383612660,
 -0.970994183180068 - I*0.003891979641355,
  0.631673431621760 - I*0.388176190581570,
 -0.658806653082993 - I*0.181991146610252,
 -1.292184403231214 + I*0.056331639648751,
  0.678822859777843 - I*0.676436048107697,
  0.049872508514131 - I*1.244805644515130,
 -0.824788327258352 + I*0.357253163396130
    };
        
        i = manakov_fscatter_numel(D, akns_discretization);
        if (i == 0) { // size D>=2, this means unknown discretization
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto leave_fun;
        }
        

        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<2*D; i++){
            q[i] = (0.41*COS(i+1) + 0.59*I*SIN(0.28*(i+1)))*50;
        }
        
        
        // without normalization 
        ret_code = manakov_fscatter(D, q, kappa, eps_t, transfer_matrix, &deg, NULL, akns_discretization);  // with kappa =1
        misc_print_buf(deg+1,transfer_matrix,"p1");

/*        for (i = 0; i<9; i++){
                printf("transfer_matrix in manakov_fscatter_test (p%d)=\n",i);
            for (UINT j = 0; j<25; j++){
                printf("%f + i%f\n", creal(transfer_matrix[j+(i*25)]), cimag(transfer_matrix[j+(i*25)]));
            }
        }*/
    
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

        /*Note: For z[0], z[1] and z[3] result and result_exact match, but the values for z[2] and z[4] are different.
        If we change D, some of the other values with z=z[i] match but not all.
        The values for z an p passed to poly_eval are okay, and when evaluating using horners method in matlab
        result and result_exact do agree, so there is probably a mistake in the matlab test file
        */
        if (misc_rel_err(9*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        
        // with normalization
        W_ptr = &W;
        ret_code = manakov_fscatter(D, q, 1, eps_t, transfer_matrix, &deg, W_ptr, akns_discretization); // with kappa = 1

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
        
        for (i=0; i<9; i++){    // replaced 4 by 9
            for (j=0; j<nz; j++)
                result[i*nz+j] = z[j];
            
            ret_code = poly_eval(deg, transfer_matrix+i*(deg+1), nz, result+i*nz);
            if (ret_code != SUCCESS){
                return E_SUBROUTINE(ret_code);
                goto leave_fun;
            }
        }
        
#ifdef DEBUG
        printf("error with normalization = %2.1e < %2.1e\n",misc_rel_err(4*nz, result, result_exact),err_bnd);
#endif
        if (misc_rel_err(4*nz, result, result_exact) > err_bnd)
            return E_TEST_FAILED;
        
        
        leave_fun:
            free(transfer_matrix);
            
            return SUCCESS;
}

INT main()
{

    if (manakov_fscatter_test_2split3A() != SUCCESS)
        return EXIT_FAILURE;
    
    printf("SUCCES");
    return EXIT_SUCCESS;
}
