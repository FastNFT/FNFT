/*
File tesing the fnft_manakov_testcases function
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakov_testcases.h"
#include <stdio.h>

INT main(){
    printf("Goes into main loop\n");

    manakov_testcases_t tc = manakov_testcases_SECH_FOCUSING;
    INT D = 16;
    REAL error_bounds[7] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN };
    fnft_manakov_opts_t opts = fnft_manakov_default_opts();
    INT ret_code = manakov_testcases_test_fnft(tc, D, error_bounds, &opts);
    printf("ret_code = %d\n",ret_code);    
  

/*
   // Check if segmentation fault is in manakov_testcases. answer: no
   manakov_testcases_t tc = manakov_testcases_SECH_FOCUSING;
   INT D = 16;
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


   ret_code = manakov_testcases(tc, D, &q1, &q2, T, &M, &contspec_exact, &ab_exact,
        XI, &K_exact,&bound_states_exact, &normconsts_exact, &residues_exact,
        &kappa);
*/
    printf("ret_code = %d\n",ret_code);


    return 0;
}

