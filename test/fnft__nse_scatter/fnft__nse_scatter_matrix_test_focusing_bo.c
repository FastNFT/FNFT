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
* Shrinivas Chimmalgi (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT nse_scatter_matrix_test_focusing_bo()
{
    UINT i, D = 8;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8];
    COMPLEX * r = NULL;
    COMPLEX result[2*8];
    COMPLEX lam[2] = {2, 1+0.5*I};
    COMPLEX result_exact[16] = { \
            
         -0.498549252810195 -     0.807291550799286*I, \
        -0.0696038174888031 +     0.308032957849424*I, \
         0.0696038174888031 +     0.308032957849424*I, \
         -0.498549252810195 +     0.807291550799286*I, \
         -0.863696939697574 +     0.499208026880278*I, \
        -0.0277240013123306 -    0.0958298043252835*I, \
         0.0277240013123306 -    0.0958298043252835*I, \
         -0.863696939697574 -     0.499208026880278*I, \
          0.783403403735771 -      1.37034786314873*I, \
       -0.00603468101196356 +     0.369501136447118*I, \
         0.0546061165237757 +     0.419105385345959*I, \
          0.255920041779507 +     0.470188584111973*I, \
          -1.45563260652517 -     0.850404464483267*I, \
        -0.0118133743801815 -    0.0497326883518552*I, \
         0.0927410716579911 -    0.0565939759073677*I, \
         -0.506785105670114 +     0.299370915096908*I };
    /* Matlab code to generate result_exact:
    eps_t = 0.13;
    kappa = +1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    result_exact = []; lam = [2 , 1+0.5i];
    for l = lam
        S=eye(4);
        for n=D:-1:1
            ks=(-kappa*abs(q(n))^2-l^2);
            k=sqrt(ks);        
            U=[cosh(k*eps_t)-1i*l*sinh(k*eps_t)/k,q(n)*sinh(k*eps_t)/k;...
                -kappa*conj(q(n))*sinh(k*eps_t)/k,cosh(k*eps_t)+1i*l*sinh(k*eps_t)/k];
            Udash=[1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t+1i+1i*l^2/ks)*sinh(k*eps_t)/k,...
                -q(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks;...
                kappa*conj(q(n))*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks,...
                -1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t-1i-1i*l^2/ks)*sinh(k*eps_t)/k];
            T=[U,zeros(2);Udash,U];
            S=S*T;
        end
        result_exact = [result_exact,[S(1,1) S(1,2) S(2,1) S(2,2) S(3,1) S(3,2) S(4,1) S(4,2)]];
    end
    format long g; result_exact.'
    */
    for (i=0; i<D; i++)
        q[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));

    ret_code = nse_scatter_matrix(D, q, r, eps_t, +1, 2, lam, result, nse_discretization_BO);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(16, result, result_exact) > 10*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (nse_scatter_matrix_test_focusing_bo() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
