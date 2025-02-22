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
* Sander Wahls (KIT) 2023, 2025.
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
    FNFT_UINT D;
    FNFT_COMPLEX * q;
    FNFT_REAL * T;
    FNFT_REAL * E;
    FNFT_UINT i;
    FNFT_INT k;
    FNFT_UINT L;
    FNFT_REAL * DEL, * al21;
    double *re;
    char msg[128]; // buffer for error messages
    fnft_kdvp_opts_t opts;
    int ret_code;

    /* Check types and dimensions of the first three inputs: q, T, E */
    if (nrhs < 3)
        mexErrMsgTxt("At least three inputs expected.");
    if ( mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q should be a real row vector.");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 2 )
        mexErrMsgTxt("Third input E should be a double 1x2 vector.");

    D = mxGetNumberOfElements(prhs[0]);
    T = mxGetPr(prhs[1]);
    E = mxGetPr(prhs[2]);
    L = D;

    /* Check values of first four inputs */

    if ( D<2 )
        mexErrMsgTxt("Length of the first input q should be at least two.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( E[0] >= E[1] )
        mexErrMsgTxt("E(1) >= E(2).");

    /* Default options for fnft_kdvv */

    opts = fnft_kdvp_default_opts();

    /* Redirect FNFT error messages and warnings to Matlabs command window */

    fnft_errwarn_setprintf(mexPrintf);

    /* Check remaining inputs, if any */

    for (k=3; k<nrhs; k++) {

        /* Check if current input is a string as desired and convert it */
        if ( !mxIsChar(prhs[k]) ) {
            snprintf(msg, sizeof msg, "%uth input should be a string.",
                (unsigned int)(k+1));
            goto on_error;
        }
        char *str = mxArrayToString(prhs[k]);
        if ( str == NULL ) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }

        /* Try to interpret value of string input */
        if ( strcmp(str, "L") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'L' should be followed by a non-negative real scalar.");
                goto on_error;
            }
            L = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        } else if ( strcmp(str, "skip_normalization") == 0 ) {

            opts.normalization_flag = 0;

        } else if ( strcmp(str, "quiet") == 0 ) {

            fnft_errwarn_setprintf(NULL);

        } else if ( strcmp(str, "default") == 0 ) {

            // Ignore this. Useful e.g. if the user wants to pass a string variable is_quiet,
            // which is either "quiet" to turn it on or "default" to keep the default (off)

        } else {
            snprintf(msg, sizeof msg, "%uth input has invalid value.",
                (unsigned int)(k+1));
            goto on_error;
        }
    }

    /* Allocate memory */

    q = mxMalloc(D * sizeof(FNFT_COMPLEX));
    if (q == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    DEL = mxMalloc(L * sizeof(FNFT_REAL));
    al21 = mxMalloc(L * sizeof(FNFT_REAL));
    if (DEL == NULL || al21 == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }
    
    /* Convert input */

    re = mxGetPr(prhs[0]);
    for (i=0; i<D; i++)
        q[i] = (FNFT_COMPLEX) re[i];

    /* Call the C routine */

    ret_code = fnft_kdvp_floquet(D, q, T, E, L, DEL, al21, &opts);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_kdvp_floquet failed (error code %i).",
                ret_code);
        goto on_error;
    }

    /* Allocate memory for outputs and convert results */

    plhs[0] = mxCreateDoubleMatrix(1, L, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, L, mxREAL);
    if (plhs[0] == NULL || plhs[1] == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }
 
    re = mxGetPr(plhs[0]);
    for (i=0; i<L; i++)
        re[i] = DEL[i];

    re = mxGetPr(plhs[1]);
    for (i=0; i<L; i++)
        re[i] = al21[i];

    /* Free memory that is no longer needed */

    mxFree(DEL);
    mxFree(al21);
    mxFree(q);
    return;

on_error:
    mexErrMsgTxt(msg);
}
