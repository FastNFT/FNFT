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
 * [2] Prins, P. J., & Wahls, S. (2021). An accurate O(N^2) floating point algorithm for the Crum 
 * transform of the KdV equation. Communications in Nonlinear Science and Numerical Simulation,
 * 102, Article 105782. https://doi.org/10.1016/j.cnsns.2021.105782
 */


#define K 2
#define DEBUG
#define QUADRATIC_ERROR_SUM_TOLERANCE 1e-27
#define MAX_QUADRATIC_ERROR 1e-28

static REAL analytic_signal(REAL t, REAL x){
    return 12.0*(3.0 + 4.0 * COSH(2.0*t-8.0*x) + COSH(4.0*t-64.0*x))/POW((3.0*COSH(t-28.0*x) + COSH(3.0*t-36.0*x)), 2);
}


INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX * t_grid = NULL;
    COMPLEX * q_1 = NULL;
    COMPLEX * q_2 = NULL;
    COMPLEX * analytic_q_1 = NULL;
    COMPLEX * analytic_q_2 = NULL;

    UINT D = 256;
    REAL const T[2] = {-8.0, 12.0};

    t_grid = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t_grid, ret_code, leave_fun);
    
    COMPLEX const eps_t = (T[1] - T[0])/(D - 1);
    
    for (UINT n=0; n<D; n++) {
        t_grid[n]= T[0] + n*eps_t;
    }

    // Calculating norming constants for both x-values.
    // Attention: both norming constants arrays need to be calculated here! If there are further
    // commands between the calculation of both norming constants arrays, the processor calculates the
    // 2. norming constants with rounded intermediate results. 
    
    COMPLEX const bound_states_1[K] = { I*2.0, I*1.0 };
    
    REAL const x_1 = -0.1;
    COMPLEX normconsts_1[K] = { CREAL(CEXP(8.0*I*CPOW(bound_states_1[0], 3)*x_1)), 
                                -CREAL(CEXP(8.0*I*CPOW(bound_states_1[1], 3)*x_1)) };

    REAL const x_2 = 0.5;                                
    COMPLEX normconsts_2[K] = { CREAL(CEXP(8.0*I*CPOW(bound_states_1[0], 3)*x_2)), 
                                -CREAL(CEXP(8.0*I*CPOW(bound_states_1[1], 3)*x_2)) };


    // Value #1 for x                                
    q_1 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q_1, ret_code, leave_fun);
    
    ret_code = fnft_kdvv_inverse(0, NULL, NULL, K, bound_states_1, normconsts_1, D, q_1, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);
    
    analytic_q_1 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(analytic_q_1, ret_code, leave_fun);
    
    for (UINT i=0; i < D; i++) {
        analytic_q_1[i] = analytic_signal(t_grid[i], x_1);
    }
    
    REAL error_1 = 0.0;
    REAL max_quadratic_error_1 = 0.0;
    REAL quadratic_error_sum_1 = 0.0;
        
    for (UINT i=0; i < D; i++) {
        error_1 = CABS(CPOW((analytic_q_1[i] - q_1[i]), 2));
        quadratic_error_sum_1 += error_1;
        if (error_1 > max_quadratic_error_1) {max_quadratic_error_1 = error_1;}
    }

    #ifdef DEBUG
        // misc_print_buf(K, bound_states_1, "bs1");
        // misc_print_buf(D, q_1, "q_1");
        printf("resulting max quadratic error: %e \n", max_quadratic_error_1);
        printf("resulting quadratic error sum: %e \n", quadratic_error_sum_1);
    #endif


    // Value #2 for x
    q_2 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q_2, ret_code, leave_fun);

    analytic_q_2 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(analytic_q_2, ret_code, leave_fun);
    
    ret_code = fnft_kdvv_inverse(0, NULL, NULL, K, bound_states_1, normconsts_2, D, q_2, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);
    
    for (UINT i=0; i < D; i++) {
        analytic_q_2[i] = analytic_signal(t_grid[i], x_2);
    }

    REAL error_2 = 0.0;
    REAL max_quadratic_error_2 = 0.0;
    REAL quadratic_error_sum_2 = 0.0;
        
    for (UINT i=0; i < D; i++) {
        error_2 = CABS(CPOW((analytic_q_2[i] - q_2[i]), 2));
        quadratic_error_sum_2 += error_2;
        if (error_2 > max_quadratic_error_2) {max_quadratic_error_2 = error_2;}
    }

    #ifdef DEBUG
        // misc_print_buf(D, q_2, "q_2");
        printf("resulting max quadratic error: %e \n", max_quadratic_error_2);
        printf("resulting quadratic error sum: %e \n", quadratic_error_sum_2);
    #endif

    UINT is_quadratic_error_sum_in_tolerance =  (quadratic_error_sum_1 < QUADRATIC_ERROR_SUM_TOLERANCE) &&
                                                (quadratic_error_sum_2 < QUADRATIC_ERROR_SUM_TOLERANCE);
    UINT is_max_quadratic_error_in_tolerance =  (max_quadratic_error_1 < MAX_QUADRATIC_ERROR) &&
                                                (max_quadratic_error_2 < MAX_QUADRATIC_ERROR);

    
    if (is_max_quadratic_error_in_tolerance &&
        is_quadratic_error_sum_in_tolerance) 
    {
        ret_code = SUCCESS;
    } 
    else {
        ret_code = FNFT_EC_TEST_FAILED;
    }
       

leave_fun:
    free(t_grid);
    free(q_1);
    free(q_2);
    free(analytic_q_1);
    free(analytic_q_2);

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
