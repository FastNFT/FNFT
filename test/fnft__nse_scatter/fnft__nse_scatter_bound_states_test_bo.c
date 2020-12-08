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
* Sander Wahls (TU Delft) 2017-2018, 2020.
* Shrinivas Chimmalgi (TU Delft) 2017-2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#ifdef DEBUG
#include <stdio.h> // for printf
#endif

INT nse_scatter_bound_states_test_bo()
{
    UINT i, D = 256;
    INT ret_code;
    REAL eps_t, T[2] = {-16,16};
    REAL errs[3], error_bounds[3];
    COMPLEX q[256], r[256];
    COMPLEX a_vals[3], aprime_vals[3], b_vals[3];
    COMPLEX bound_states[3] = {0.5*I, 1.5*I, 2.5*I};
    eps_t = (T[1] - T[0])/(D - 1);

    COMPLEX a_vals_exact[3] = {
         1.863969116250583e-05 +   0.0*I,
        -2.098244573838071e-05 +   0.0*I,
         4.682264057040961e-05 +   0.0*I  };
    COMPLEX aprime_vals_exact[3] = {
         0.0 -    0.333127126071093*I,
         0.0 +    0.0416102597562244*I,
         0.0 -    0.0333984387408527*I };
    COMPLEX b_vals_exact[3] = {
        -0.999850812020245 + 0.0*I,
          1.00269448509783 + 0.0*I,
        -0.998504589821944 + 0.0*I };
        /* Matlab code to generate result_exact:
	D=256; T = [-16,16]; bound_states =[0.5i ,1.5i, 2.5i];
	eps_t = (T(2)-T(1))/(D-1);
	t = T(1):eps_t:T(2); q = 3*sech(t);
	Tend=T(2)+eps_t/2;
	T1 = T(1)-eps_t/2;
	a=[];b=[];aprime=[];
	i0 = 1;
	i1 = D;
	norm_left = 0.0;
	norm_right = 0.0;
	while (i0 < i1)
    	    if (norm_left < norm_right)
        	i0=i0+1;
        	norm_left=norm_left+ eps_t*abs(q(i0));
    	    else
        	i1=i1-1;
        	norm_right=norm_right+ eps_t*abs(q(i0));
    	    end
	end
	split_point=i0;

	for j=1:1:numel(bound_states)
    		l=bound_states(j);
   		 SR=eye(4);
    		for n=length(q):-1:split_point+1
    	    		ks=(-abs(q(n))^2-l^2);
       			 k=sqrt(ks);
        		U=[cosh(k*eps_t)-1i*l*sinh(k*eps_t)/k,q(n)*sinh(k*eps_t)/k;...
            		-conj(q(n))*sinh(k*eps_t)/k,cosh(k*eps_t)+1i*l*sinh(k*eps_t)/k];
        		Udash=[1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t+1i+1i*l^2/ks)*sinh(k*eps_t)/k,...
            		-q(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks;...
           		 conj(q(n))*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks,...
            		-1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t-1i-1i*l^2/ks)*sinh(k*eps_t)/k];
       		 	T=[U,zeros(2);Udash,U];
        		SR=SR*T;
    		end
   	 	SL=eye(4);
    		for n=split_point:-1:1
        		ks=(-abs(q(n))^2-l^2);
        		k=sqrt(ks);
        		U=[cosh(k*eps_t)-1i*l*sinh(k*eps_t)/k,q(n)*sinh(k*eps_t)/k;...
           		 -conj(q(n))*sinh(k*eps_t)/k,cosh(k*eps_t)+1i*l*sinh(k*eps_t)/k];
        		Udash=[1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t+1i+1i*l^2/ks)*sinh(k*eps_t)/k,...
           		 -q(n)*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks;...
           		 conj(q(n))*l*(eps_t*cosh(k*eps_t)-sinh(k*eps_t)/k)/ks,...
           		 -1i*eps_t*l^2*cosh(k*eps_t)/ks-(l*eps_t-1i-1i*l^2/ks)*sinh(k*eps_t)/k];
        		T=[U,zeros(2);Udash,U];
        		SL=SL*T;
    		end
    		TM = SR*SL;
    		a(j)=(TM(1,1))*exp(1i*l*Tend)*exp(-1i*l*T1);
    		aprime(j)= (TM(3,1)+1i*Tend*(TM(1,1)+TM(3,3)))*exp(1i*l*Tend)*exp(-1i*l*T1);
    		b(j)=(SL(2,1)/SR(1,1))*exp(-1i*l*T1)*exp(-1i*l*Tend);
	end
	format long g; a_vals_exact = a.'
	aprime_vals_exact = aprime.'
	bvals_exact = b.'
    */
    for (i=0; i<D; i++) {
        q[i] = 3.0*misc_sech(T[0] + i*eps_t);
        r[i] = -CONJ(q[i]);
    }
    ret_code = nse_scatter_bound_states(D, q, r, T, 3,
        bound_states, a_vals, aprime_vals, b_vals, nse_discretization_BO, 0);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    error_bounds[0] = 100*EPSILON;
    error_bounds[1] = 100*EPSILON;
    error_bounds[2] = 100*EPSILON;
    errs[0] = misc_rel_err(3, a_vals, a_vals_exact);
    errs[1] = misc_rel_err(3, aprime_vals, aprime_vals_exact);
    errs[2] = misc_rel_err(3, b_vals, b_vals_exact);

#ifdef DEBUG
    printf(
        "fnft__nse_scatter_bound_states_test_bo: %2.1e<=%2.1e %2.1e<=%2.1e %2.1e<=%2.1e\n",
        errs[0], error_bounds[0],
        errs[1], error_bounds[1],
        errs[2], error_bounds[2]
    );
#endif

    if (!(errs[0] <= error_bounds[0]))
        ret_code = E_TEST_FAILED;
    if (!(errs[1] <= error_bounds[1]))
        ret_code = E_TEST_FAILED;
    if (!(errs[2] <= error_bounds[2]))
        ret_code = E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (nse_scatter_bound_states_test_bo() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
