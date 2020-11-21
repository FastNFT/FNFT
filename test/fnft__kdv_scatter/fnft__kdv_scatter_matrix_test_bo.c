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
         -0.544605303960817 -     0.958647947076152*I, \
         0.0794721386915349 +    0.0134467976453307*I, \
        -0.0246337974340963 +    0.0123132398678595*I, \
         -0.447572503318363 +      0.78665621605389*I, \
         -0.966803401931833 +     0.607395950461071*I, \
        -0.0639071462390028 +   0.00196680566143484*I, \
         0.0622976106588603 -    0.0146064121524493*I, \
         -0.836459620311112 -     0.421161672558673*I, \
          0.844221235130024 -      1.75354130915637*I, \
          0.147492933060637 -    0.0819705993458431*I, \
        -0.0881078629988439 +     0.119509172610443*I, \
          0.210672525960293 +     0.467024487542091*I, \
          -1.69348275358045 -     0.643630324529962*I, \
          -0.10214305870232 +     0.141516722852452*I, \
          0.106220590826741 -     0.150214497063317*I, \
         -0.391717762925894 +     0.218025257136358*I };

    /* Matlab code to generate result_exact:
    eps_t = 0.13;
    D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D)); r = -ones(1,D);
    result_exact = []; lam = [2 , 1+0.5i]; U=nan(2); Udash=nan(2);
    sincdash = @(x) (cos(x)-sinc(x./pi))./x;
    for l = lam
        S=eye(4);
        for n=D:-1:1
            Delta = eps_t * sqrt(l^2+q(n));

            U(1,1) =  (q(n) * eps_t / (2i*l) - 1i * l * eps_t) * sinc(Delta/pi) + cos(Delta);
            U(1,2) =   q(n) * eps_t / (2i*l) * sinc(Delta/pi);
            U(2,1) =  -U(1,2);
            U(2,2) = -(q(n) * eps_t / (2i*l) - 1i * l * eps_t) * sinc(Delta/pi) + cos(Delta);

            Udash(1,1) = -1i * eps_t * (1+q(n)/(-2*l^2)) * sinc(Delta/pi) ...
                         + -1i * (eps_t.^3/Delta) * (l.^2 + q(n)/2) * sincdash(Delta) ...
                         - eps_t.^2 * l * sinc(Delta/pi);
            Udash(1,2) = q(n).*eps_t.^3./(2i*Delta^2) * ...
                                           ( cos(Delta) - (2+q(n)./l.^2) .* sinc(Delta/pi) );
            Udash(2,1) = - Udash(1,2);
            Udash(2,2) =  1i * eps_t * (1+q(n)/(-2*l^2)) * sinc(Delta/pi) ...
                         +  1i * (eps_t.^3/Delta) * (l.^2 + q(n)/2) * sincdash(Delta) ...
                         - eps_t.^2 * l * sinc(Delta/pi);
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
    ret_code = kdv_scatter_matrix(D, q, r, eps_t, kappa, 2, lam, result, kdv_discretization_BO_VANILLA, 1);

    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);

    #ifdef DEBUG
        misc_print_buf(16, result, "result");
        misc_print_buf(16, result_exact, "resxct");
    #endif

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
