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



INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX bound_states_i[5] = {I*SQRT(1.0/2.0), I*SQRT(2.0/2.0), I*SQRT(3.0/2.0), I*SQRT(4.0/2.0), I*SQRT(5.0/2.0)};
    COMPLEX normconsts_i[5] = {1*10, -1*0.1, 1*1, -1*1e-5, 1*1e7};

    fnft_kdvv_params kdvv_parameters = {
        .D = 256,
        .T = {-20.0, 20.0},
        .K = 5,
        .bound_states = bound_states_i,
        .normconsts = normconsts_i,
        .M = 10,
        .XI = {-2.0, 2.0},
        .contspec = malloc(10*sizeof(COMPLEX))
    };

    
    REAL err_bnd_bound_states = 1.4e-3;
    REAL err_bnd_spurious_bound_states = 1e-2;
    REAL err_bnd_normconst = 4.5e-2;
    REAL err_bnd_contspec = 1e-1;

    
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




