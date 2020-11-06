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
* Sander Wahls (TU Delft) 2017-2018.
* Shrinivas Chimmalgi (TU Delft) 2017-2018,2020.
* Peter J Prins (TU Delft) 2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdv_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT kdv_scatter_matrix_test_bo()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8], r[8];
    COMPLEX result[2*8];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[16] = { \

         -0.624077442652352 -     0.972094744721483*I, \
        -0.0537871905813229 +      0.31788855476614*I, \
         -0.436609430226878 +    0.0502846841920214*I, \
         -0.368100364626828 +     0.800103013699221*I, \
         -0.902896255692831 +     0.605429144799636*I, \
        -0.0347608179364008 -    0.0966843075729412*I, \
          0.471300816414904 -    0.0241075859152959*I, \
         -0.900366766550114 -     0.419194866897238*I, \
          0.696728302069386 -      1.67157070981052*I, \
          0.016448265631049 +     0.376956465467117*I, \
         -0.887223992518918 +     0.244638039704334*I, \
           0.35816545902093 +     0.385053888196248*I, \
          -1.59133969487813 -     0.785147047382414*I, \
        -0.0169491883108973 -    0.0508169741358181*I, \
          0.369249525725664 +     0.117437867995581*I, \
         -0.493860821628215 +      0.35954197998881*I };
            
    /* Matlab code to generate result_exact:
    eps_t = 0.13;
    D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D)); r = -ones(1,D);
    result_exact = []; lam = [2 , 1+0.5i];
    for l = lam
        S=eye(4);
        for n=D:-1:1
            ks=(q(n)*r(n)-l^2);
            k=sqrt(ks);        
            U=[cosh(k*eps_t)-1i*l*sinh(k*eps_t)/k,q(n)*sinh(k*eps_t)/k;...
                r(n)*sinh(k*eps_t)/k,cosh(k*eps_t)+1i*l*sinh(k*eps_t)/k];
            Udash=[1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t+1i+1i*l^2/ks)*sinh(k*eps_t)/k,...
                -q(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks;...
                -r(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks,...
                -1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t-1i-1i*l^2/ks)*sinh(k*eps_t)/k];
            T=[U,zeros(2);Udash,U];
            S=S*T;
        end
        result_exact = [result_exact,[S(1,1) S(1,2) S(2,1) S(2,2) S(3,1) S(3,2) S(4,1) S(4,2)]];
    end
    format long g; result_exact.'
    */
    for (i=0; i<D; i++) {
        q[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));
        r[i] = -1.0;
    }

    INT kappa=1; // unused
    ret_code = kdv_scatter_matrix(D, q, r, eps_t, kappa, 2, lam, result, kdv_discretization_BO, 1);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(16, result, result_exact) > 10*EPSILON)
        return E_TEST_FAILED;
    
    return SUCCESS;
}

INT main()
{
    if (kdv_scatter_matrix_test_bo() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
