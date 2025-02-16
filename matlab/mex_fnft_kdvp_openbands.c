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
* Sander Wahls (KIT) 2025.
*/

#include <string.h>
#include "mex.h"
#ifndef SKIP_MATRIX_H
#include "matrix.h"
#else
#include <stdio.h>
#endif
#include "fnft_kdvp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FNFT_UINT k = 0, K = 0, len = 0;
    FNFT_REAL * main_spec = NULL;
    FNFT_REAL * openbands = NULL;
    double * re = NULL;
    char msg[128]; // buffer for error messages
    int ret_code;

    /* Check types and dimensions of the first input: main_spec */

    if (nrhs < 1)
        mexErrMsgTxt("At least one input expected.");
    if ( mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q should be a real row vector.");

    /* Check first input */

    len = mxGetNumberOfElements(prhs[0]);
    if ( len<2 )
        mexErrMsgTxt("Length of the first input main_spec should be at least two.");
    if ( len%2 )
        mexErrMsgTxt("Length of the first input main_spec should be even.");

    K = len/2;
    main_spec = mxGetPr(prhs[0]);
    
    /* Redirect FNFT error messages and warnings to Matlabs command window */

    fnft_errwarn_setprintf(mexPrintf);

    /* Allocate memory */

    openbands = mxMalloc(3*K * sizeof(mxREAL));
    if (openbands == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }
  
    /* Call the C routine */

    ret_code = fnft_kdvp_openbands(&K, main_spec, openbands);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_kdvp_openbands failed (error code %i).",
                ret_code);
        goto on_error;
    }

    /* Allocate memory for output and convert result */

    plhs[0] = mxCreateDoubleMatrix(1, 2*K, mxREAL);
    if (plhs[0] == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    re = mxGetPr(plhs[0]);
    for (k=0; k<3*K; k++)
        re[k] = openbands[k];

   return;

on_error:
    mexErrMsgTxt(msg);
}
