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
* Sander Wahls (KIT) 2026.
*/

#include <string.h>
#include <stdio.h>
#include "mex.h"
#ifndef SKIP_MATRIX_H
#include "matrix.h"
#endif
#include "fnft_kdvv_inverse.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FNFT_UINT D;
    FNFT_COMPLEX * q = NULL;
    FNFT_REAL * T;
    FNFT_UINT M, K;
    FNFT_COMPLEX * contspec = NULL;
    FNFT_COMPLEX * bound_states = NULL;
    FNFT_COMPLEX * norming_constants = NULL;
    FNFT_REAL * XI = NULL;
    FNFT_UINT i;
    double *re, *im;
    char msg[128]; // buffer for error messages
    FNFT_INT ret_code;
    FNFT_UINT k;

    if (nlhs < 1)
        return;

    /* Check types and dimensions of the first seven inputs: contspec, XI,
       bound_states, norming_constants, D, T */

    if ( nrhs < 6 )
        mexErrMsgTxt("At least seven inputs expected.");
    if ( !mxIsEmpty(prhs[0]) )
        mexErrMsgTxt("First input contspec should be empty. Dealing with continuous spectrum has not yet been implemented!");
    if ( !mxIsEmpty(prhs[1]) )
        mexErrMsgTxt("Second input XI should be empty. Dealing with continuous spectrum has not yet been implemented!");
    if ( !mxIsEmpty(prhs[2]) && (!mxIsDouble(prhs[2]) || !mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != 1) )
        mexErrMsgTxt("Third input bound_states should be a complex row vector (double precision) or []. Try passing complex(double(bound_states(:)')).");
    if ( !mxIsEmpty(prhs[3]) && (!mxIsDouble(prhs[3]) || !mxIsComplex(prhs[3]) || mxGetM(prhs[3]) != 1) )
        mexErrMsgTxt("Fourth input norming_constants should be a complex row vector (double precision) or []. Try passing complex(double(norming_constants(:)')).");
    if ( mxIsComplex(prhs[4]) || !mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1 )
        mexErrMsgTxt("Fifth input D should be a real scalar (double precision).");
    if ( mxIsComplex(prhs[5]) || !mxIsDouble(prhs[5]) || mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 2 )
        mexErrMsgTxt("Sixth input T should be a real 1x2 vector (double precision).");

    M = mxGetNumberOfElements(prhs[0]);
    K = mxGetNumberOfElements(prhs[2]);
    T = mxGetPr(prhs[5]);
    D = (unsigned int)mxGetScalar(prhs[4]);
    // Dealing with continuous spectrum has not yet been implemented! XI should be empty
    // XI = mxGetPr(prhs[1]);

    /* Check values of first four inputs */

    if ( K != mxGetNumberOfElements(prhs[3]) )
        mexErrMsgTxt("bound_states and norming_constants should have the same lengths.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    // Dealing with continuous spectrum has not yet been implemented! XI should be empty
    // if ( XI[0] >= XI[1] )
    //     mexErrMsgTxt("XI(1) >= XI(2).");
    if ( D<2 )
        mexErrMsgTxt("D < 2.");

    /* Redirect FNFT error messages and warnings to Matlabs command window */

    fnft_errwarn_setprintf(mexPrintf);

    /* Check remaining inputs, if any */

    for (k=7; k<(FNFT_UINT)nrhs; k++) {

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
        if ( strcmp(str, "quiet") == 0 ) {

            fnft_errwarn_setprintf(NULL);

        } else {
            snprintf(msg, sizeof msg, "%uth input has invalid value.",
                 (unsigned int)(k+1));
            goto on_error;
        }
    }

    /* Allocate memory */

    q = mxMalloc(D * sizeof(FNFT_COMPLEX));
    if (M>0)
        contspec = mxMalloc(M * sizeof(FNFT_COMPLEX));
    if (K>0) {
        bound_states = mxMalloc(K * sizeof(FNFT_COMPLEX));
        norming_constants = mxMalloc(K * sizeof(FNFT_COMPLEX));
    }
    if ( q == NULL || (M>0 && contspec == NULL) || (K>0 && bound_states == NULL)
        || (K>0 && norming_constants == NULL) ) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    /* Convert inputs */

    re = mxGetPr(prhs[0]);
    im = mxGetPi(prhs[0]);
    for (i=0; i<M; i++)
        contspec[i] = re[i] + I*im[i];
    re = mxGetPr(prhs[2]);
    im = mxGetPi(prhs[2]);
    for (i=0; i<K; i++)
        bound_states[i] = re[i] + I*im[i];
    re = mxGetPr(prhs[3]);
    im = mxGetPi(prhs[3]);
    for (i=0; i<K; i++)
        norming_constants[i] = re[i] + I*im[i];

    /* Call the C routine */

    ret_code = fnft_kdvv_inverse(M, contspec, XI, K, bound_states,
                                 norming_constants, D, q, T, NULL);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_kdvv_inverse failed (error code %i).",
                ret_code);
        goto on_error;
    }

    /* Allocate memory for the output */

    plhs[0] = mxCreateDoubleMatrix(1, D, mxCOMPLEX);

    /* Allocate memory for outputs and convert results */

    re = mxGetPr(plhs[0]);
    im = mxGetPi(plhs[0]);
    for (i=0; i<D; i++) {
        re[i] = FNFT_CREAL(q[i]);
        im[i] = FNFT_CIMAG(q[i]);
    }

    /* Free memory that is no longer needed */

    mxFree(q);
    mxFree(contspec);
    mxFree(bound_states);
    mxFree(norming_constants);
    return;

on_error:
    mexErrMsgTxt(msg);
}
