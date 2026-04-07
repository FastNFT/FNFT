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
* Fabian Fischer (Hiwi KIT) 2026
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft__kdvv_inverse_testcases.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

/* This is a testcase with analytic reference */


/* This testcase is based on the example discussed on p. 74-78 in "Solitons - an introduction" by 
* P. G. Drazin, R. S. Johnson, Cambridge University Press (1989)
* 
* - the eigenvalues are k_1 = 1 and k_2 = 2
* - the eigenfunctions are: 
*    > psi_1 ~ sqrt(6)*exp(-t) as t->+infinity and psi_1 ~ sqrt(6)*exp(t) as t->-infinity
*       ==> normconst_1 = sqrt(6)

* the values for the normconsts are 
* 
* Note: The usage of x and t is inverted in the example towards the usage in this library
* - the resulting signal has to be multiplied by -1
why COMPLEX normconsts[K] = { CEXP(64.0*x), -CEXP(8.0*x) } instead of
COMPLEX normconsts[K] = { 2.0*SQRT(3.0)*CEXP(32.0*x), -SQRT(6.0)*CEXP(4.0*x) };????
*/


#define K 2
#define DEBUG

REAL analytic_signal(REAL t, REAL x){
    return 12.0*(3.0 + 4.0 * COSH(2.0*t-8.0*x) + COSH(4.0*t-64.0*x))/POW((3.0*COSH(t-28.0*x) + COSH(3.0*t-36.0*x)), 2);
}


INT main()
{
    INT ret_code = SUCCESS;

    REAL const x = 0.1;

    UINT D = 256;
    UINT M = 10;
    REAL XI[2] = {-2.0, 2.0};
    REAL T[2] = {-10.0, 10.0};

    COMPLEX * contspec = NULL;
    contspec = malloc(10*sizeof(COMPLEX));
    CHECK_NOMEM(contspec, ret_code, leave_fun);

    COMPLEX * q = NULL;
    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);

    COMPLEX bound_states[K] = { I*2.0, I*1.0 };
    COMPLEX normconsts[K] = {   CEXP(8.0*I*CPOW(bound_states[0], 3)*x), 
                                -CEXP(8.0*I*CPOW(bound_states[1], 3)*x) };
    
    ret_code = fnft_kdvv_inverse(M, contspec, XI, K, bound_states, normconsts, D, q, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);

    COMPLEX * t_grid = NULL;
    t_grid = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t_grid, ret_code, leave_fun);

    const COMPLEX eps_t = (T[1] - T[0])/(D - 1);

    for (UINT n=0; n<D; n++) {
        t_grid[n]= T[0] + n*eps_t;
    }

    COMPLEX * analytic_q_a = NULL;
    analytic_q_a = malloc(D*sizeof(COMPLEX));

    for (UINT i=0; i < D; i++) {
        analytic_q_a[i] = analytic_signal(t_grid[i], x);
    }

    #ifdef DEBUG
        misc_print_buf(D, t_grid, "t_grid");
        misc_print_buf(D, q, "output_inverse");
        misc_print_buf(D, analytic_q_a, "analytic_q");
    #endif

    REAL analytic_q = 0.0;
    REAL error = 0.0;
    REAL max_quadratic_error = 0.0;
    REAL quadratic_error_sum = 0.0;
    

    for (UINT i=0; i < D; i++) {
        error = POW((analytic_q_a[i] - q[i]), 2);
        quadratic_error_sum += error;
        if (error > max_quadratic_error) {max_quadratic_error = error;}
    }

    #ifdef DEBUG
        printf("resulting max quadratic error: %e \n", max_quadratic_error);
        printf("resulting quadratic error sum: %e \n", quadratic_error_sum);
    #endif

    UINT is_quadratic_error_sum_in_tolerance = quadratic_error_sum < 1e-27;
    UINT is_max_quadratic_error_in_tolerance = max_quadratic_error < 1e-28;

    
    if (is_max_quadratic_error_in_tolerance &&
        is_quadratic_error_sum_in_tolerance) 
    {
        ret_code = SUCCESS;
    } 
    else {
        ret_code = FNFT_EC_TEST_FAILED;
    }

leave_fun:
    free(contspec);
    free(q);
    free(analytic_q_a);

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
