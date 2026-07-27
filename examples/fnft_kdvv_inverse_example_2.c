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

// This is an example for the usage of the inverse kdvv transform with many 
// bound states

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft_kdvv_inverse.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

#define K 19    // Number of bound states to add

INT main()
{
    INT ret_code = SUCCESS;

    COMPLEX * bound_states = NULL;
    COMPLEX * q = NULL;

    COMPLEX desired_solitions_height[K] = { 40.0, 39.0, 38.0, 37.0, 36.0,
                                            30.0, 29.0, 28.0, 27.0, 26.0,
                                            20.0, 19.0, 18.0, 17.0, 16.0,
                                            10.0, 9.0, 8.0, 7.0};

    COMPLEX normconsts[K] = {   1e20, -1e-7, 1e5, -1e3, 1e1,
                                -1e0, 1e2, -1e4, 1e-6, -1e8,
                                1e2, -1e4, 1e6, -1e8, 1e-10,
                                -1e7, 1e-6, -1e5, 1e-9};

    // resulting bound states out of desired solitions height                                
    bound_states = malloc(K * sizeof(COMPLEX));
    CHECK_NOMEM(bound_states, ret_code, leave_fun);
                                
    for (UINT i = 0; i<K; i++){
        bound_states[i] = I*SQRT(desired_solitions_height[i]/2.0);
    }

    // General parameters
    UINT D = 256;
    REAL T[2] = {-18.0, 12.0};
    
    // allocating memory for the output of the inverse kdvv
    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);

    // call inverse kdvv transform
    ret_code = fnft_kdvv_inverse(0, NULL, NULL, K, bound_states, normconsts, D, q, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);

    // print output of inverse kdvv into console
    misc_print_buf(D, q, "output");

leave_fun:
    free(bound_states);
    free(q);
}
