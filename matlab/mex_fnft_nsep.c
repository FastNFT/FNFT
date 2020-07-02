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
* Sander Wahls (TU Delft) 2017-2018, 2020.
* Shrinivas Chimmalgi (TU Delft) 2020.
*/

#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "fnft_nsep.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t D;
    double complex * q;
    double * T;
    double phase_shift;
    size_t M;
    double complex * main_spec;
    size_t K;
    double complex * aux_spec = NULL;
    int * sheet_indices = NULL;
    int kappa;
    int upsampling_factor = 1;
    size_t i;
    ptrdiff_t k;
    double *re, *im;
    double *msr, *msi, *asr, *asi;
    char msg[128]; // buffer for error messages
    fnft_nsep_opts_t opts;
    int ret_code;

    /* To suppress unused parameter warning */
    (void) nlhs;

    /* Check types and dimensions of the first three inputs: q, T, kappa */
    if (nrhs < 3)
        mexErrMsgTxt("At least three inputs expected.");
    if ( !mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q should be a complex row vector. Try passing complex(q).");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1 )
        mexErrMsgTxt("Third input phase_shift should be a scalar.");
    if ( !mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1 )
        mexErrMsgTxt("Fourth input kappa should be a scalar.");

    D = mxGetNumberOfElements(prhs[0]);
    K = D;
    M = D;
    T = mxGetPr(prhs[1]);
    phase_shift = mxGetScalar(prhs[2]);
    kappa = (int)mxGetScalar(prhs[3]);

    /* Check values of first three inputs */
    if ( D<2 || D%2 != 0 )
        mexErrMsgTxt("Length of the first input q should be even.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( kappa != +1 && kappa != -1 )
        mexErrMsgTxt("Fourth input kappa should be +1.0 or -1.0.");

    // Default options for fnft_nsep
    opts = fnft_nsep_default_opts();

    /* Redirect FNFT error messages and warnings to Matlabs command window */
    fnft_errwarn_setprintf(mexPrintf);

    /* Check remaining inputs, if any */
    for (k=4; k<nrhs; k++) {

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
        if ( strcmp(str, "filt_none") == 0 ) {

            opts.filtering = fnft_nsep_filt_NONE;

        } else if ( strcmp(str, "filt_manual") == 0 ) {

            opts.filtering = fnft_nsep_filt_MANUAL;

            /* Extract bounding box for manual filltering */
            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
            || mxGetM(prhs[k+1]) != 1 || mxGetN(prhs[k+1]) != 4) {
                snprintf(msg, sizeof msg, "'filt_manual' should be followed by a real row vector of length four. See the help.");
                goto on_error;
            }
            double const * const tmp = mxGetPr(prhs[k+1]);
            opts.bounding_box[0] = (FNFT_REAL)tmp[0];
            opts.bounding_box[1] = (FNFT_REAL)tmp[1];
            opts.bounding_box[2] = (FNFT_REAL)tmp[2];
            opts.bounding_box[3] = (FNFT_REAL)tmp[3];

            /* Increase k to account for bounding box vector */
            k++;

        } else if ( strcmp(str, "points_per_spine") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'points_per_spine' should be followed by non-negative real number. See the help.");
                goto on_error;
            }
            opts.points_per_spine = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        } else if ( strcmp(str, "loc_mixed") == 0 ) {

            opts.localization = fnft_nsep_loc_MIXED;

        } else if ( strcmp(str, "loc_subsample_and_refine") == 0 ) {

            opts.localization = fnft_nsep_loc_SUBSAMPLE_AND_REFINE;

        } else if ( strcmp(str, "loc_gridsearch") == 0 ) {

            opts.localization = fnft_nsep_loc_GRIDSEARCH;

        } else if ( strcmp(str, "loc_max_evals") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'loc_max_evals' should be followed by non-negative real number. See the help.");
                goto on_error;
            }
            opts.max_evals = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        } else if ( strcmp(str, "loc_Dsub") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'loc_Dsub' should be followed by non-negative real number. See the help.");
                goto on_error;
            }
            opts.Dsub = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        } else if ( strcmp(str, "quiet") == 0 ) {

            fnft_errwarn_setprintf(NULL);
            
        } else if ( strcmp(str, "discr_modal") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT2_MODAL;
            
        } else if ( strcmp(str, "discr_2split2A") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT2A;
            
        } else if ( strcmp(str, "discr_2split4A") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT4A; upsampling_factor = 4;
            
        } else if ( strcmp(str, "discr_2split4B") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT4B; upsampling_factor = 2;
            
        } else if ( strcmp(str, "discr_4split4B") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_4SPLIT4B; upsampling_factor = 4;            

        } else {
            snprintf(msg, sizeof msg, "%uth input has invalid value.",
                (unsigned int)(k+1));
            goto on_error;
        }
    }

    /* Allocate memory */
    q = mxMalloc(D * sizeof(double complex));
    K = (opts.points_per_spine)*upsampling_factor*D + 1;
    M = upsampling_factor*D;
    main_spec = mxMalloc(K * sizeof(double complex));
    aux_spec = mxMalloc(M * sizeof(double complex));
    sheet_indices = mxMalloc(M * sizeof(int));
    if ( q == NULL || main_spec == NULL || aux_spec == NULL
    || sheet_indices == NULL ) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }

    /* Convert input */
    re = mxGetPr(prhs[0]);
    im = mxGetPi(prhs[0]);
    for (i=0; i<D; i++)
        q[i] = re[i] + I*im[i];

    /* Call the C routine */
    ret_code = fnft_nsep(D, q, T, phase_shift, &K, main_spec, &M, aux_spec, NULL, kappa,
        &opts);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_nsep failed (error code %i).",
            ret_code);
        goto on_error;
    }

    /* Allocate memory for the outputs */
    plhs[0] = mxCreateDoubleMatrix(1, K, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(1, M, mxCOMPLEX);

    /* Convert outputs */
    msr = mxGetPr(plhs[0]);
    msi = mxGetPi(plhs[0]);
    asr = mxGetPr(plhs[1]);
    asi = mxGetPi(plhs[1]);
    for (i=0; i<K; i++) {
        msr[i] = creal(main_spec[i]);
        msi[i] = cimag(main_spec[i]);
    }
    for (i=0; i<M; i++) {
        asr[i] = creal(aux_spec[i]);
        asi[i] = cimag(aux_spec[i]);
    }

    /* Free memory that is no longer needed */
    mxFree(q);
    mxFree(main_spec);
    mxFree(aux_spec);
    mxFree(sheet_indices);
    return;

on_error:
    mexErrMsgTxt(msg);
}
