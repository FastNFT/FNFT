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
* Shrinivas Chimmalgi (TU Delft) 2019.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__misc.h"
#include "fnft__fft_wrapper.h"

static INT misc_resample_test()
{

    UINT i, j, D = 1256;
    COMPLEX *q = NULL;
    COMPLEX *q_s = NULL;
    COMPLEX *q_s_exact = NULL;
    REAL *t = NULL;
    REAL delta[4] = {-0.25,0.1387,0.036,-0.015};

    INT ret_code = SUCCESS;   
    t = malloc(D * sizeof(REAL));
    q = malloc(D * sizeof(COMPLEX));
    q_s = malloc(D * sizeof(COMPLEX));
    q_s_exact = malloc(D * sizeof(COMPLEX));

    if (q == NULL || q_s == NULL || q_s_exact == NULL || t == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    REAL eps_t = 64.0/(D-1);
    for (i=0; i<D; i++) {
        t[i] = -32.0 + i*eps_t;
        q[i] = 5.23345*misc_sech(t[i])*CEXP(I*5.72254*t[i]);
    }

    for (j=0; j<4; j++) {
        ret_code = misc_resample(D, eps_t, q, delta[j], q_s);
        CHECK_RETCODE(ret_code, leave_fun);
        for (i=0; i<D; i++)
            q_s_exact[i] = 5.23345*misc_sech(t[i]+delta[j])*CEXP(I*5.72254*(t[i]+delta[j]));
        REAL err = misc_rel_err(D, q_s, q_s_exact);
#ifdef DEBUG
    printf("error = %g\n", err);
#endif
        if (err > 3e-7)
           return E_TEST_FAILED;
    }

leave_fun:
    free(q);
    free(q_s);
    free(q_s_exact);
    free(t);
    return ret_code;
}

INT main()
{
    if ( misc_resample_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
