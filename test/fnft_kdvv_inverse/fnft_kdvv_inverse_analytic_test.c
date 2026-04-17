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


/* This is a testcase with analytic reference
*
* This test case is based on the example discussed on p. 74-78 in [1]. Similar to this example
* is the example discussed in section 5.2 (p. 13) in [2]. Note that we use different notation and 
* the definition of the norming constants in [2]. Thus we have to do the following adaptions of the
* formulas given in [1]:
*
*   - substitute x with t and t with x
*   - the analytic signal needs to be multiplied by -1
*   - norming constants:
*       b(2i, x) = exp(8i*(2i)^3*x) --> b(2i, 0) = 1
*       b(1i, x) = -exp(8i*(1i)^3*x) --> b(1i, 0) = -1
*    
* [1] P. G. Drazin, R. S. Johnson (1989). Solitons - an introduction. Cambridge University Press
* [2] Prins, P. J., & Wahls, S. (2021). An accurate O(N 2) floating point algorithm for the Crum 
* transform of the KdV equation. Communications in Nonlinear Science and Numerical Simulation,
* 102, Article 105782. https://doi.org/10.1016/j.cnsns.2021.105782
*/


#define K 2
#define DEBUG
#define QUADRATIC_ERROR_SUM_TOLERANCE 1e-27
#define MAX_QUADRATIC_ERROR 1e-28

REAL analytic_signal(REAL t, REAL x){
    return 12.0*(3.0 + 4.0 * COSH(2.0*t-8.0*x) + COSH(4.0*t-64.0*x))/POW((3.0*COSH(t-28.0*x) + COSH(3.0*t-36.0*x)), 2);
}


INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX * contspec = NULL;
    COMPLEX * q = NULL;
    COMPLEX * t_grid = NULL;
    COMPLEX * analytic_q_a = NULL;

    REAL const x = 0.5;

    UINT D = 256;
    UINT M = 10;
    REAL XI[2] = {-2.0, 2.0};
    REAL T[2] = {-8.0, 12.0};

    contspec = malloc(10*sizeof(COMPLEX));
    CHECK_NOMEM(contspec, ret_code, leave_fun);

    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);

    COMPLEX bound_states[K] = { I*2.0, I*1.0 };
    COMPLEX normconsts[K] = {   CEXP(8.0*I*CPOW(bound_states[0], 3)*x), 
                                -CEXP(8.0*I*CPOW(bound_states[1], 3)*x) };
    
    ret_code = fnft_kdvv_inverse(M, contspec, XI, K, bound_states, normconsts, D, q, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);

    t_grid = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t_grid, ret_code, leave_fun);

    const COMPLEX eps_t = (T[1] - T[0])/(D - 1);

    for (UINT n=0; n<D; n++) {
        t_grid[n]= T[0] + n*eps_t;
    }

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

    UINT is_quadratic_error_sum_in_tolerance = quadratic_error_sum < QUADRATIC_ERROR_SUM_TOLERANCE;
    UINT is_max_quadratic_error_in_tolerance = max_quadratic_error < MAX_QUADRATIC_ERROR;

    
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
    free(t_grid);
    free(analytic_q_a);

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
