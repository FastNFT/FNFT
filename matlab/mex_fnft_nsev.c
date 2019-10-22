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
#include "fnft_nsev.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FNFT_UINT D;
    FNFT_COMPLEX * q;
    FNFT_REAL * T;
    FNFT_UINT M;
    FNFT_COMPLEX * contspec = NULL;
    FNFT_REAL * XI;
    FNFT_UINT K;
    FNFT_COMPLEX * bound_states = NULL;
    FNFT_COMPLEX * normconsts_or_residuals = NULL;
    FNFT_INT kappa;
    FNFT_INT skip_contspec_flag = 0;
    FNFT_INT skip_bound_states_flag = 0;
    FNFT_INT skip_normconsts_flag = 0;
    FNFT_UINT i, j;
    FNFT_INT k;
    double *re, *im;
    double *csr, *csi, *bsr, *bsi, *ncr, *nci;
    char msg[128]; // buffer for error messages
    fnft_nsev_opts_t opts;
    int ret_code;

    /* Check number of inputs to avoid computing results that have not been
    requested. We treat nlhs==0 like nlhs==1 because it that case, the result
    will be stored in Matlabs ans variable. */

    if (nlhs < 2)
        skip_bound_states_flag = 1;
    if (nlhs < 3)
        skip_normconsts_flag = 1;
    
    /* Check types and dimensions of the first four inputs: q, T, XI, kappa */

    if (nrhs < 4)
        mexErrMsgTxt("At least four inputs expected.");
    if ( !mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q should be a complex row vector. Try passing complex(q).");
    if ( !mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2 )
        mexErrMsgTxt("Second input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 2 )
        mexErrMsgTxt("Third input XI should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1 )
        mexErrMsgTxt("Fourth input kappa should be a scalar.");
    
    D = mxGetNumberOfElements(prhs[0]);
    K = D;
    M = D;
    T = mxGetPr(prhs[1]);
    XI = mxGetPr(prhs[2]);
    kappa = (int)mxGetScalar(prhs[3]);
    
    /* Check values of first four inputs */

    if ( D<2 )
        mexErrMsgTxt("Length of the first input q should be at least two.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( XI[0] >= XI[1] )
        mexErrMsgTxt("XI(1) >= XI(2).");
    if ( kappa != +1 && kappa != -1 )
        mexErrMsgTxt("Fourth input kappa should be +1.0 or -1.0.");
    
    /* Default options for fnft_nsev */

    opts = fnft_nsev_default_opts();
    
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
        if ( strcmp(str, "M") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'M' should be followed by a non-negative real scalar.");
                goto on_error;
            }
            M = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        } else if ( strcmp(str, "bsloc_fasteigen") == 0 ) {
            
            opts.bound_state_localization = fnft_nsev_bsloc_FAST_EIGENVALUE;
            
        } else if ( strcmp(str, "bsloc_newton") == 0 ) {
            
            opts.bound_state_localization = fnft_nsev_bsloc_NEWTON;
            
            /* Extract initial guesses for Newtons method */
            if ( k+1 == nrhs || !mxIsComplex(prhs[k+1])
            || mxGetM(prhs[k+1]) != 1 || mxGetN(prhs[k+1]) < 1) {
                snprintf(msg, sizeof msg, "'bsloc_newton' should be followed by a complex row vector of initial guesses for Newton's method. Try passing complex(...).");
                goto on_error;
            }
            K = mxGetN(prhs[k+1]);
            bound_states = mxMalloc(K * sizeof(FNFT_COMPLEX));
            if (bound_states == NULL) {
                snprintf(msg, sizeof msg, "Out of memory.");
                goto on_error;
            }
            re = mxGetPr(prhs[k+1]);
            im = mxGetPi(prhs[k+1]);
            for (j=0; j<K; j++)
                bound_states[j] = re[j] + I*im[j];
            
            /* Increase k to account for vector of initial guesses */
            k++;
            
        } else if ( strcmp(str, "bsloc_niter") == 0 ) {
            
            /* Extract desired number of iterations */
            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
            || mxGetNumberOfElements(prhs[k+1]) != 1
                    || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'bsloc_niter' should be followed by a non-negative real scalar.");
                goto on_error;
            }
            opts.niter = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            
            /* Increase k to account for vector of initial guesses */
    	    k++;
            
        } else if ( strcmp(str, "bsloc_subsamp_refine") == 0 ) {
            
            opts.bound_state_localization = fnft_nsev_bsloc_SUBSAMPLE_AND_REFINE;

        } else if ( strcmp(str, "bsloc_Dsub") == 0 ) {

            /* Extract desired number of iterations */
            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
            || mxGetNumberOfElements(prhs[k+1]) != 1
                    || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'bsloc_Dsub' should be followed by a non-negative real scalar.");
                goto on_error;
            }
            opts.Dsub = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            
            /* Increase k to account for vector of initial guesses */
    	    k++;

            
        } else if ( strcmp(str, "bsfilt_none") == 0 ) {
            
            opts.bound_state_filtering = fnft_nsev_bsfilt_NONE;
            
        } else if ( strcmp(str, "bsfilt_basic") == 0 ) {
            
            opts.bound_state_filtering = fnft_nsev_bsfilt_BASIC;
            
        } else if ( strcmp(str, "bsfilt_full") == 0 ) {
            
            opts.bound_state_filtering = fnft_nsev_bsfilt_FULL;
            
        } else if ( strcmp(str, "discr_modal") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT2_MODAL;
            
        } else if ( strcmp(str, "discr_2split2A") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT2A;
            
        } else if ( strcmp(str, "discr_2split4A") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT4A;
            
        } else if ( strcmp(str, "discr_2split4B") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_2SPLIT4B;
            
        } else if ( strcmp(str, "discr_4split4B") == 0 ) {
            
            opts.discretization = fnft_nse_discretization_4SPLIT4B; 
            
        } else if ( strcmp(str, "RE") == 0 ) {
            
            opts.richardson_extrapolation_flag  = 1;
            
        } else if ( strcmp(str, "dstype_residues") == 0 ) {
            
            opts.discspec_type = fnft_nsev_dstype_RESIDUES;
 
        } else if ( strcmp(str, "cstype_ab") == 0 ) {
            
            opts.contspec_type = fnft_nsev_cstype_AB;
 
        } else if ( strcmp(str, "skip_cs") == 0 ) {
            
            skip_contspec_flag = 1;

        } else if ( strcmp(str, "skip_bs") == 0 ) {
            
            skip_bound_states_flag = 1;
            skip_normconsts_flag = 1; // since bound states are needed to
                                      // compute norming constants

        } else if ( strcmp(str, "skip_nc") == 0 ) {
            
            skip_normconsts_flag = 1;
          
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
   
    if (skip_contspec_flag == 0) {
        if (opts.contspec_type == fnft_nsev_cstype_AB)
            contspec = mxMalloc(2*M * sizeof(FNFT_COMPLEX));
        else
            contspec = mxMalloc(M * sizeof(FNFT_COMPLEX));      
        if (contspec == NULL) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }
    }

    if (skip_bound_states_flag == 0) {
        if (bound_states == NULL) {
            K = fnft_nsev_max_K(D, &opts);
            if (K == 0) {
                snprintf(msg, sizeof msg, "fnft_nsev_max_K returned zero.");
                goto on_error;
            }
            bound_states = mxMalloc(K * sizeof(FNFT_COMPLEX));
        }
        if (bound_states == NULL) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }
    }

    if (skip_normconsts_flag == 0) {
        normconsts_or_residuals = mxMalloc(K * sizeof(FNFT_COMPLEX));
        if (normconsts_or_residuals == NULL) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }
    }
    
    /* Convert input */

    re = mxGetPr(prhs[0]);
    im = mxGetPi(prhs[0]);
    for (i=0; i<D; i++)
        q[i] = re[i] + I*im[i];
    
    /* Call the C routine */
    
    ret_code = fnft_nsev(D, q, T, M, contspec, XI, &K, bound_states,
            normconsts_or_residuals, kappa, & opts);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_nsev failed (error code %i).",
                ret_code);
        goto on_error;
    }
    
    /* Allocate memory for the outputs */

    if (opts.contspec_type == fnft_nsev_cstype_AB)
        plhs[0] = mxCreateDoubleMatrix(1, 2*M, mxCOMPLEX);
    else
        plhs[0] = mxCreateDoubleMatrix(1, M, mxCOMPLEX);
    
    /* Allocate memory for outputs and convert results */

    if (skip_contspec_flag == 0) {
        csr = mxGetPr(plhs[0]);
        csi = mxGetPi(plhs[0]);
        if (opts.contspec_type == fnft_nsev_cstype_AB) {
            for (i=0; i<2*M; i++) {
                csr[i] = FNFT_CREAL(contspec[i]);
                csi[i] = FNFT_CIMAG(contspec[i]);
            }
        } else {
            for (i=0; i<M; i++) {
                csr[i] = FNFT_CREAL(contspec[i]);
                csi[i] = FNFT_CIMAG(contspec[i]);
            }
        }
    } else {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
    }

    if (skip_bound_states_flag == 0) {
        plhs[1] = mxCreateDoubleMatrix(1, K, mxCOMPLEX);
        bsr = mxGetPr(plhs[1]);
        bsi = mxGetPi(plhs[1]);
        for (i=0; i<K; i++) {
            bsr[i] = FNFT_CREAL(bound_states[i]);
            bsi[i] = FNFT_CIMAG(bound_states[i]);
        }
    } else if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
    }

    if (skip_normconsts_flag == 0) {
        plhs[2] = mxCreateDoubleMatrix(1, K, mxCOMPLEX);
        ncr = mxGetPr(plhs[2]);
        nci = mxGetPi(plhs[2]);
        for (i=0; i<K; i++) {
            ncr[i] = FNFT_CREAL(normconsts_or_residuals[i]);
            nci[i] = FNFT_CIMAG(normconsts_or_residuals[i]);
        }
    } else if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);

    }
    
    /* Free memory that is no longer needed */

    mxFree(q);
    mxFree(contspec);
    mxFree(bound_states);
    mxFree(normconsts_or_residuals);
    return;
    
on_error:
    mexErrMsgTxt(msg);
}
