// Just a file to see if manakov_fscatter runs

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft.h"
#include "fnft__manakov_fscatter.h"
#include "fnft__poly_eval.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

INT main(){
    printf("goes into main\n");
    INT ret_code;
    UINT deg;
    INT D = 1;      // one sample
    COMPLEX q[2] = {1, 2};
    INT kappa = 1;
    REAL eps_t = 0.1;
    COMPLEX *transfer_matrix = NULL;

    akns_discretization_t akns_discretization = akns_discretization_2SPLIT3A;
    INT i = manakov_fscatter_numel(D, akns_discretization);
    transfer_matrix = malloc(i*sizeof(COMPLEX));
    
    ret_code = manakov_fscatter(D, q, kappa, eps_t, transfer_matrix, &deg, NULL, akns_discretization);  // with kappa =1

    printf("compile test1\n");

    return 0;
}

