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

// This is a simple example for the usage of the inverse kdvv transform in C

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft_kdvv_inverse.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"

#define K 5     // number of bound states to add

INT main()
{
    INT ret_code = SUCCESS;
    
    COMPLEX * q = NULL;
    COMPLEX * contspec = NULL;

    // five desired bound states and norming constants
    // - bound states needs to be true imaginary positive
    // - every norming constants belong to the bound state at the same index
    // - the signs of the norming constants have to alternate regarding the order
    //   of bound states. The norming constant of the biggest bound state has to be positive.
    COMPLEX bound_states[K] = { I*5.0, I*4.0, I*3.0, I*2.0, I*1.0 };
    COMPLEX norming_constants[K] = { 1e5, -1e-2, 1e0, -1e3, 1e1 };

    // General parameters
    UINT D = 256;                   // Number of samples of the computed output
    REAL T[2] = {-10.0, 10.0};      // area for which the output should be computed
    
    // allocation of memory for the computed output
    q = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(q, ret_code, leave_fun);
    
    // Define continuous spectrum
    // continuous spectrum is out of function, but needs to be defined/allocated for using 
    // the inverse kdvv (state 04/2026)
    UINT M = 10;
    REAL XI[2] = {-2.0, 2.0};

    contspec = malloc(M * sizeof(COMPLEX));
    CHECK_NOMEM(contspec, ret_code, leave_fun);

    // call of the inverse kdvv
    ret_code = fnft_kdvv_inverse(M, contspec, XI, K, bound_states, norming_constants, D, q, T, NULL);
    CHECK_RETCODE(ret_code, leave_fun);

    // prints the results in the console
    misc_print_buf(D, q, "output");

leave_fun:
    free(q);
    free(contspec);
}
