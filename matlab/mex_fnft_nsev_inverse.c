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
    FNFT_COMPLEX * q = NULL;
    FNFT_REAL * T;
    FNFT_UINT M;
    FNFT_COMPLEX * contspec = NULL;
    FNFT_REAL * XI;
    FNFT_INT kappa;
    FNFT_UINT i;
    double *re, *im, *csr, *csi;
    char msg[128]; // buffer for error messages
    fnft_nsev_inverse_opts_t opts;
    FNFT_INT ret_code;
    FNFT_INT k;

    if (nlhs < 1)
        return;

    /* Check types and dimensions of the first four inputs: contspec, XI, D, T,
       kappa */

    if (nrhs < 5)
        mexErrMsgTxt("At least five inputs expected.");
    if ( !mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input contspec should be a complex row vector. Try passing complex(q).");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input XI should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 )
        mexErrMsgTxt("Second input D should be a scalar.");
    if ( !mxIsDouble(prhs[3]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[3]) != 2 )
        mexErrMsgTxt("Fourth input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1 )
        mexErrMsgTxt("Fifth input kappa should be a scalar.");

    M = mxGetNumberOfElements(prhs[0]);
    T = mxGetPr(prhs[1]);
    D = (unsigned int)mxGetScalar(prhs[2]);
    XI = mxGetPr(prhs[3]);
    kappa = (int)mxGetScalar(prhs[4]);

    /* Check values of first four inputs */

    if ( M<2 )
        mexErrMsgTxt("Length of the first input contspec should be at least two.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( XI[0] >= XI[1] )
        mexErrMsgTxt("XI(1) >= XI(2).");
    if ( D<2 )
        mexErrMsgTxt("D < 2.");
    if ( kappa != +1 && kappa != -1 )
        mexErrMsgTxt("Fourth input kappa should be +1.0 or -1.0.");

    /* Default options for fnft_nsev */

    opts = fnft_nsev_inverse_default_opts();

    /* Redirect FNFT error messages and warnings to Matlabs command window */

    fnft_errwarn_setprintf(mexPrintf);

    /* Check remaining inputs, if any */

    for (k=5; k<nrhs; k++) {

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
        if ( strcmp(str, "cstype_reflection_coefficient") == 0 ) {

            opts.contspec_type =
                fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT;

        } else if ( strcmp(str, "cstype_B_of_tau") == 0 ) {

            opts.contspec_type = fnft_nsev_inverse_cstype_B_OF_TAU;

        } else if ( strcmp(str, "csmethod_tfmatrix_contains_refl_coeff") == 0 ){

            opts.contspec_inversion_method =
                fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_REFL_COEFF;

        } else if ( strcmp(str, "csmethod_tfmatrix_contains_ab_from_iter") == 0){

            opts.contspec_inversion_method =
                fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER;

        } else if ( strcmp(str, "discr_modal") == 0 ) {

            opts.discretization = fnft_nse_discretization_2SPLIT2_MODAL;

        } else if ( strcmp(str, "oversampling_factor") == 0 ) {

            /* Extract desired oversampling factor */
            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 1 ) {
                snprintf(msg, sizeof msg, "'oversampling_factor' should be followed by a positive natural number.");
                goto on_error;
            }
            opts.oversampling_factor = (FNFT_UINT)mxGetScalar(prhs[k+1]);


        } else if ( strcmp(str, "quiet") == 0 ) {

            fnft_errwarn_setprintf(NULL);

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
    contspec = mxMalloc(M * sizeof(FNFT_COMPLEX));
    if (contspec == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    /* Convert input */

    re = mxGetPr(prhs[0]);
    im = mxGetPi(prhs[0]);
    for (i=0; i<M; i++)
        contspec[i] = re[i] + I*im[i];

    /* Call the C routine */

    ret_code = fnft_nsev_inverse(M, contspec, XI, 0, NULL, NULL, D, q, T, kappa, &opts);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_nsev_inverse failed (error code %i).",
                ret_code);
        goto on_error;
    }

    /* Allocate memory for the output */

    plhs[0] = mxCreateDoubleMatrix(1, D, mxCOMPLEX);

    /* Allocate memory for outputs and convert results */

    csr = mxGetPr(plhs[0]);
    csi = mxGetPi(plhs[0]);
    for (i=0; i<D; i++) {
        csr[i] = FNFT_CREAL(q[i]);
        csi[i] = FNFT_CIMAG(q[i]);
    }

    /* Free memory that is no longer needed */

    mxFree(q);
    mxFree(contspec);
    return;

on_error:
    mexErrMsgTxt(msg);
}
