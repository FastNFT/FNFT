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

#define K_I 25

INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX bound_states_i[K_I] = { 50.0, 49.0, 48.0, 47.0, 46.0,
                                    40.0, 39.0, 38.0, 37.0, 36.0,
                                    30.0, 29.0, 28.0, 27.0, 26.0,
                                    20.0, 19.0, 18.0, 17.0, 16.0,
                                    10.0, 9.0, 8.0, 7.0, 6.0};

    COMPLEX normconsts_i[K_I] = {   1e20, -1e-7, 1e5, -1e3, 1e1,
                                    -1e0, 1e2, -1e4, 1e-6, -1e8,
                                    1e2, -1e4, 1e6, -1e8, 1e-10,
                                    -1e7, 1e-6, -1e5, 1e-9, -1e-11,
                                    1e-3, -1e3, 1e2, -1e2, 1e1};

    for (UINT i = 0; i<K_I; i++){
        bound_states_i[i] = I*SQRT(bound_states_i[i]/2.0);
    }


    fnft_kdvv_params kdvv_parameters = {
        .D = 512,
        .T = {-15.0, 15.0},
        .K = K_I,
        .bound_states = bound_states_i,
        .normconsts = normconsts_i,
        .M = 10,
        .XI = {-2.0, 2.0},
        .contspec = malloc(10*sizeof(COMPLEX))
    };

    
    REAL err_bnd_bound_states = 5e-3;
    REAL err_bnd_spurious_bound_states = 0.2;
    REAL err_bnd_normconst = 2.5e-1;
    REAL err_bnd_contspec = 2e-1;

    
    ret_code = kdvv_testcases_get_spectrum_of_inverse(  kdvv_parameters, err_bnd_bound_states, 
                                                        err_bnd_spurious_bound_states, 
                                                        err_bnd_normconst, err_bnd_contspec);

    CHECK_RETCODE(ret_code, leave_fun);

    
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
