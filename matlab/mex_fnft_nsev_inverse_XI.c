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
#include "fnft_nsev_inverse.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FNFT_UINT D;
    FNFT_REAL * T;
    FNFT_UINT M;
    double * XI;
    char msg[128]; // buffer for error messages
    fnft_nsev_inverse_opts_t opts;
    FNFT_INT ret_code;
    FNFT_UINT i;

    if (nlhs < 1)
        return;

    /* Check types and dimensions of the first four inputs: D, T, M */

    if (nrhs < 3)
        mexErrMsgTxt("At least three inputs expected.");
    if ( !mxIsDouble(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 1 )
        mexErrMsgTxt("Fourth input D should be a scalar.");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 )
        mexErrMsgTxt("Third input M should be a scalar.");

    D = (int)mxGetScalar(prhs[0]);
    T = mxGetPr(prhs[1]);
    M = (int)mxGetScalar(prhs[2]);

    /* Allocate memory for the output */

    opts = fnft_nsev_inverse_default_opts();
    plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    XI = mxGetPr(plhs[0]);
    ret_code = fnft_nsev_inverse_XI(D, T, M, XI, opts.discretization);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_nsev_inverse_XI failed (error code %i).",
                 ret_code);
        goto on_error;
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(1, M, mxREAL);
        double * const xi = mxGetPr(plhs[1]);
        const double eps_xi = (XI[1] - XI[0])/(M - 1);
        for (i=0; i<M; i++)
            xi[i] = XI[0] + i*eps_xi;
    }

    return;

on_error:
    mexErrMsgTxt(msg);
}
