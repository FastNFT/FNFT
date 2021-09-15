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
* Contributor:
* Lianne de Vries (TU Delft student) 2021.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__manakovv_testcases.h"
#include "fnft_manakovv.h"

// TODO: debugger was not able to find fnft__manakovv_testcases.h, even though includePath was set to right folder.
// fnft__manakovv_testcases.h was copied to workspace folder. Should be removed.

INT main(){
    printf("Goes into main loop\n");

    manakovv_testcases_t tc = manakovv_testcases_SECH_FOCUSING;
    UINT D = 16;
    REAL error_bounds[7] = {200, 200, 200, 200, 200, 200, 200};
//    REAL error_bounds[7] = {
//        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN };
    fnft_manakovv_opts_t opts = fnft_manakovv_default_opts();
    INT ret_code = manakovv_testcases_test_fnft(tc, D, error_bounds, &opts);
    printf("ret_code = %d\n",ret_code);    
  

/*
   // Check if segmentation fault is in manakovv_testcases. answer: no
   manakovv_testcases_t tc = manakovv_testcases_SECH_FOCUSING;
   UINT D = 16;
   COMPLEX * q1 = NULL;
    COMPLEX * q2 = NULL;
    COMPLEX * contspec = NULL;
    COMPLEX * bound_states = NULL;
    COMPLEX * normconsts_and_residues = NULL;
    REAL T[200], XI[200];
    COMPLEX * contspec_exact = NULL;
    COMPLEX * ab_exact = NULL;
    COMPLEX * bound_states_exact = NULL;
    COMPLEX * normconsts_exact = NULL;
    COMPLEX * residues_exact = NULL;
    UINT K, K_exact=0, M;
    INT kappa = 0;
    REAL errs[7] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN };
    INT ret_code;


   ret_code = manakovv_testcases(tc, D, &q1, &q2, T, &M, &contspec_exact, &ab_exact,
        XI, &K_exact,&bound_states_exact, &normconsts_exact, &residues_exact,
        &kappa);
*/
    printf("ret_code = %d\n",ret_code);


    return 0;
}

