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
* Sander Wahls (TU Delft) 2017-2018.
*/

#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "fnft_kdvv.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t D;
    double complex * q;
    double * T;
    size_t M;
    double complex * contspec;
    double * XI;
    size_t i;
    double *re, *im;
    double *csr, *csi;
    char msg[128]; // buffer for error messages
    int ret_code;

    /* To suppress unused parameter warning */
    (void) nlhs;

    /* Check types and dimensions of the first three inputs: q, T, XI */
    if (nrhs < 3)
        mexErrMsgTxt("At least three inputs expected.");
    if ( !mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q should be a complex row vector. Try passing complex(q).");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 2 )
        mexErrMsgTxt("Third input XI should be a double 1x2 vector.");

    D = mxGetNumberOfElements(prhs[0]);
    M = D;
    T = mxGetPr(prhs[1]);
    XI = mxGetPr(prhs[2]);

    /* Check values of first four inputs */
    if ( D<2 || (D & (D-1)) != 0 )
        mexErrMsgTxt("Length of the first input q should be a positive power of two.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( XI[0] >= XI[1] )
        mexErrMsgTxt("XI(1) >= XI(2).");

    /* Redirect FNFT error messages and warnings to Matlabs command window */
    fnft_errwarn_setprintf(mexPrintf);

    /* Allocate memory */
    q = mxMalloc(D * sizeof(double complex));
    contspec = mxMalloc(D * sizeof(double complex));
    if ( q == NULL || contspec == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    /* Convert input */
    re = mxGetPr(prhs[0]);
    im = mxGetPi(prhs[0]);
    for (i=0; i<D; i++)
        q[i] = re[i] + I*im[i];

    /* Call the C routine */
    ret_code = fnft_kdvv(D, q, T, M, contspec, XI, NULL, NULL, NULL, NULL);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_kdvv failed (error code %i).",
            ret_code);
        goto on_error;
    }

    /* Allocate memory for the outputs */
    plhs[0] = mxCreateDoubleMatrix(1, D, mxCOMPLEX);

    /* Convert outputs */
    csr = mxGetPr(plhs[0]);
    csi = mxGetPi(plhs[0]);
    for (i=0; i<D; i++) {
        csr[i] = creal(contspec[i]);
        csi[i] = cimag(contspec[i]);
    }

    /* Free memory that is no longer needed */
    mxFree(q);
    mxFree(contspec);
    return;

on_error:
    mexErrMsgTxt(msg);
}
