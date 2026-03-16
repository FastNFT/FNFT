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

/* This is a testcase with an assymetric window */

#define K_I 8
#define DEBUG

INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX bound_states_i[K_I] = { 40.0, 35.0, 28.0, 23.0, 19.0,
                                    16.0, 12.0, 4.0};

    COMPLEX normconsts_i[K_I] = {   1e-3, -1e-17, 1e-15, -1e-20, 1e-1,
                                    -1e-5, 1e-12, -1e-14};

    for (UINT i = 0; i<K_I; i++){
        bound_states_i[i] = I*SQRT(bound_states_i[i]/2.0);
    }


    fnft_kdvv_params kdvv_parameters = {
        .D = 256,
        .T = {-18.0, 5.0},
        .K = K_I,
        .bound_states = bound_states_i,
        .normconsts = normconsts_i,
        .M = 10,
        .XI = {-2.0, 2.0},
        .contspec = malloc(10*sizeof(COMPLEX))
    };

    #ifdef DEBUG
        printf("\n Initial parameters: \n");
        kdvv_print_spectrum(kdvv_parameters.bound_states, kdvv_parameters.normconsts, 
                            kdvv_parameters.contspec, kdvv_parameters.XI, 
                            kdvv_parameters.M, kdvv_parameters.D, 
                            kdvv_parameters.K);
    #endif



    #ifdef DEBUG
        printf("\n -- 1. Iteration: -- \n");
    #endif
    
    REAL err_bnd_bound_states = 3.8e-3;
    REAL err_bnd_spurious_bound_states = 0.12;
    REAL err_bnd_normconst = 2e-1;
    REAL err_bnd_contspec = 2.5e-1;
   
    ret_code = kdvv_testcases_get_spectrum_of_inverse(  kdvv_parameters, err_bnd_bound_states, 
                                                        err_bnd_spurious_bound_states, 
                                                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);

    

    #ifdef DEBUG
        printf("\n -- 2. Iteration: -- \n");
    #endif

    // Check quadratic convergence
    kdvv_parameters.D *= 2;
    err_bnd_bound_states /= 4;
    err_bnd_spurious_bound_states /= 4;
    err_bnd_normconst /= 4;
    err_bnd_contspec /= 4;

    ret_code = kdvv_testcases_get_spectrum_of_inverse(  kdvv_parameters, err_bnd_bound_states, 
                                                        err_bnd_spurious_bound_states, 
                                                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);



    #ifdef DEBUG
        printf("\n -- 3. Iteration: -- \n");
    #endif
    
    kdvv_parameters.D *= 2;
    err_bnd_bound_states /= 4;
    err_bnd_spurious_bound_states /= 4;
    err_bnd_normconst /= 4;
    err_bnd_contspec /= 4;

    ret_code = kdvv_testcases_get_spectrum_of_inverse(  kdvv_parameters, err_bnd_bound_states, 
                                                        err_bnd_spurious_bound_states, 
                                                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);
    

leave_fun:
    

    if (ret_code != SUCCESS)
        return EXIT_FAILURE;
    else
	    return EXIT_SUCCESS;
}
