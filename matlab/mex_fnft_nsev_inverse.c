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
    FNFT_UINT q0_index = 0; // used to localize a given seed potential, if any
    FNFT_REAL * T;
    FNFT_UINT M, K;
    FNFT_COMPLEX * contspec = NULL;
    FNFT_COMPLEX * bound_states = NULL;
    FNFT_COMPLEX * normconsts_or_residues = NULL;
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

    /* Check types and dimensions of the first seven inputs: contspec, XI,
       bound_states, normconsts_or_residues, D, T, kappa */

    if ( nrhs < 7 )
        mexErrMsgTxt("At least seven inputs expected.");
    if ( !mxIsEmpty(prhs[0]) && (!mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1) )
        mexErrMsgTxt("First input contspec should be a complex row vector or []. Try passing complex(contspec).");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input XI should be a double 1x2 vector.");
    if ( !mxIsEmpty(prhs[2]) && (!mxIsComplex(prhs[2]) || mxGetM(prhs[2]) != 1) )
        mexErrMsgTxt("Third input bound_states should be a complex row vector or []. Try passing complex(bound_states).");
    if ( !mxIsEmpty(prhs[3]) && (!mxIsComplex(prhs[3]) || mxGetM(prhs[3]) != 1) )
        mexErrMsgTxt("Fourth input normconsts_or_residues should be a complex row vector or []. Try passing complex(normconsts_or_residues).");
    if ( !mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1 )
        mexErrMsgTxt("Fifth input D should be a scalar.");
    if ( !mxIsDouble(prhs[5]) || mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 2 )
        mexErrMsgTxt("Sixth input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1 )
        mexErrMsgTxt("Seventh input kappa should be a scalar.");

    M = mxGetNumberOfElements(prhs[0]);
    K = mxGetNumberOfElements(prhs[2]);
    T = mxGetPr(prhs[5]);
    D = (unsigned int)mxGetScalar(prhs[4]);
    XI = mxGetPr(prhs[1]);
    kappa = (int)mxGetScalar(prhs[6]);

    /* Check values of first four inputs */

    if ( K != mxGetNumberOfElements(prhs[3]) )
        mexErrMsgTxt("bound_states and normconsts_or_residues should have the same lengths.");
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

    for (k=7; k<nrhs; k++) {

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

        } else if ( strcmp(str, "dstype_residues") == 0 ) {

            opts.discspec_type = fnft_nsev_inverse_dstype_RESIDUES;

        } else if ( strcmp(str, "csmethod_use_seed_potential_intead") == 0 ) {

            opts.contspec_inversion_method =
                fnft_nsev_inverse_csmethod_USE_SEED_POTENTIAL_INSTEAD ;
            /* Extract initial potential */
            if ( k+1 == nrhs || !mxIsComplex(prhs[k+1])
                 || mxGetM(prhs[k+1]) != 1 || mxGetN(prhs[k+1]) != D ) {
                snprintf(msg, sizeof msg, "'csmethod_use_seed_potential_intead' should be followed by a complex 1xD vector. Try passing complex(q0).");
                goto on_error;
            }
            if ( M > 0 ) {
                snprintf(msg, sizeof msg, "Both contspec and an intial potential q0 have been provided. Pass contspec=[] when using dsmethod_addsoliton_cdt.");
                goto on_error;
            }
            q0_index = k+1;
            k++;

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
            k++;

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
    if (M>0)
        contspec = mxMalloc(M * sizeof(FNFT_COMPLEX));
    if (K>0) {
        bound_states = mxMalloc(K * sizeof(FNFT_COMPLEX));
        normconsts_or_residues = mxMalloc(K * sizeof(FNFT_COMPLEX));
    }
    if ( q == NULL || (M>0 && contspec == NULL) || (K>0 && bound_states == NULL)
        || (K>0 && normconsts_or_residues == NULL) ) {
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
        normconsts_or_residues[i] = re[i] + I*im[i];
    if (q0_index != 0) {
        re = mxGetPr(prhs[q0_index]);
        im = mxGetPi(prhs[q0_index]);
        for (i=0; i<D; i++)
            q[i] = re[i] + I*im[i];
    }

    /* Call the C routine */

    ret_code = fnft_nsev_inverse(M, contspec, XI, K, bound_states,
                                 normconsts_or_residues, D, q, T, kappa, &opts);
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
