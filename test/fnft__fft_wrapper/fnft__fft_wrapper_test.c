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
* Sander Wahls (TU Delft) 2018.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__misc.h"
#include "fnft__fft_wrapper.h"

static INT fft_wrapper_test()
{
    const UINT fft_length = 4;
    UINT i;
    COMPLEX *in = NULL;
    COMPLEX *out = NULL;
    COMPLEX in_exact[4] = { 1.0-2.0*I, 0.3+0.4*I, -2.0-2.0*I, -3.0+4.0*I };
    COMPLEX out_exact[4] = { -3.7+0.4*I, -0.6-3.3*I, 1.7-8.4*I, 6.6+3.3*I };
    fft_wrapper_plan_t plan = fft_wrapper_safe_plan_init();
    INT ret_code = SUCCESS;   

    in = fft_wrapper_malloc(fft_length * sizeof(COMPLEX));
    out = fft_wrapper_malloc(fft_length * sizeof(COMPLEX));
    if (in == NULL || out == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    for (i=0; i<fft_length; i++)
        in[i] = in_exact[i];
    ret_code = fft_wrapper_create_plan(&plan, fft_length, in, out, 0);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_execute_plan(plan, in, out);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_destroy_plan(&plan);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL err = misc_rel_err(fft_length, out, out_exact);
    if (err > 100*EPSILON)
		return E_TEST_FAILED;

    for (i=0; i<fft_length; i++)
        in[i] = out_exact[i] / fft_length;
    ret_code = fft_wrapper_create_plan(&plan, fft_length, in, out, 1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_execute_plan(plan, in, out);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_destroy_plan(&plan);
    CHECK_RETCODE(ret_code, leave_fun);

    err = misc_rel_err(fft_length, out, in_exact);
    if (err > 100*EPSILON)
		return E_TEST_FAILED;

leave_fun:
    fft_wrapper_destroy_plan(&plan);
    free(in);
    free(out);
    return ret_code;
}

INT main()
{
    if ( fft_wrapper_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
