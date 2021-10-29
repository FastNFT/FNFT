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
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
* Peter J. Prins (2021).
* Lianne de Vries (TU Delft student) 2021.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <string.h>
#include "mex.h"
#include "matrix.h"
//#include "fnft_manakov_discretization_t.h"
#include "fnft_manakovv.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    FNFT_UINT D;
    FNFT_COMPLEX * q1;
    FNFT_COMPLEX * q2;
    FNFT_REAL * T;
    FNFT_UINT M;
    FNFT_COMPLEX * contspec = NULL;
    FNFT_REAL * XI;
    FNFT_UINT K;
    FNFT_COMPLEX * bound_states = NULL;
    FNFT_COMPLEX * normconsts_or_residuals = NULL;
    FNFT_INT kappa;
    FNFT_INT skip_bound_states_flag = 0;
	FNFT_INT skip_contspec_flag = 0;
    FNFT_UINT i;
    FNFT_INT k;
    double *re_q1, *im_q1, *re_q2, *im_q2;
    double *csr, *csi, *bsr, *bsi;
    char msg[128]; // buffer for error messages
    fnft_manakovv_opts_t opts;
    int ret_code;

    /* Check number of inputs to avoid computing results that have not been
    requested. We treat nlhs==0 like nlhs==1 because it that case, the result
    will be stored in Matlabs ans variable. 
    We do not need this as we are only calculating the contspec anyway
    if (nlhs < 2)
        skip_bound_states_flag = 1;
    if (nlhs < 3)
        skip_normconsts_flag = 1;
    */
    /* Check types and dimensions of the first five inputs: q1, q2, T, XI, kappa */

    if (nrhs < 5)
        mexErrMsgTxt("At least five inputs expected.");
    if ( !mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1)
        mexErrMsgTxt("First input q1 should be a complex row vector. Try passing complex(q1).");
    if ( !mxIsComplex(prhs[1]) || mxGetM(prhs[1]) != 1)
        mexErrMsgTxt("Second input q2 should be a complex row vector. Try passing complex(q2).");
    if ( !mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 2 )
        mexErrMsgTxt("Third input T should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[3]) || mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 2 )
        mexErrMsgTxt("Fourth input XI should be a double 1x2 vector.");
    if ( !mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1 )
        mexErrMsgTxt("Fifth input kappa should be a scalar.");
    
    // Here we copy the values passed via prhs (the inputs) to variables we can use in this C file
    D = mxGetNumberOfElements(prhs[0]);     // This is the number of samples from q
    K = D;
    M = D;
    T = mxGetPr(prhs[2]);
    XI = mxGetPr(prhs[3]);
    kappa = (int)mxGetScalar(prhs[4]);
    
    /* Check values of first five inputs */

    if ( D<2 )
        mexErrMsgTxt("Length of the first input q2 should be at least two.");
    if ( !(mxGetNumberOfElements(prhs[1])==D) )       // Check if q2 has same number of elements as q1
        mexErrMsgTxt("Length of the second input q2 should be equal to length of q1.");
    if ( T[0] >= T[1] )
        mexErrMsgTxt("T(1) >= T(2).");
    if ( XI[0] >= XI[1] )
        mexErrMsgTxt("XI(1) >= XI(2).");
    if ( kappa != +1 && kappa != -1 )
        mexErrMsgTxt("Fifth input kappa should be +1.0 or -1.0.");
    
    /* Default options for fnft_manakovv */

    opts = fnft_manakovv_default_opts();
    
    /* Redirect FNFT error messages and warnings to Matlabs command window */

    fnft_errwarn_setprintf(mexPrintf);
    
    /* Check remaining inputs, if any */
    // TODO: write code for error checking here
    // Also take care of setting the opts fields according to additional inputs (cs type: rho, ab, both)
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
        if ( strcmp(str, "M") == 0 ) {

            if ( k+1 == nrhs || !mxIsDouble(prhs[k+1])
                 || mxGetNumberOfElements(prhs[k+1]) != 1
                 || mxGetScalar(prhs[k+1]) < 0.0 ) {
                snprintf(msg, sizeof msg, "'M' should be followed by a non-negative real scalar.");
                goto on_error;
            }
            M = (FNFT_UINT)mxGetScalar(prhs[k+1]);
            k++;

        }   
        else if ( strcmp(str, "RE") == 0 ) {

            opts.richardson_extrapolation_flag  = 1;

        } else if ( strcmp(str, "cstype_ab") == 0 ) {

            opts.contspec_type = fnft_manakovv_cstype_AB;

        } else if ( strcmp(str, "cstype_both") == 0 ) {

            opts.contspec_type = fnft_manakovv_cstype_BOTH;

        } else if ( strcmp(str, "skip_bs") == 0 ) {

            skip_bound_states_flag = 1;

        } else if (strcmp(str, "skip_contspec") == 0) {

			skip_contspec_flag = 1;

		} else if ( strcmp(str, "quiet") == 0 ) {

            fnft_errwarn_setprintf(NULL);
        
        // Fast discretizations
        } else if ( strcmp(str, "discr_2split3A") == 0 ) {

            opts.discretization = fnft_manakov_discretization_2SPLIT3A;

        } else if ( strcmp(str, "discr_2split3B") == 0 ) {

            opts.discretization = fnft_manakov_discretization_2SPLIT3B;

        } else if ( strcmp(str, "discr_2split4A") == 0 ) {

            opts.discretization = fnft_manakov_discretization_2SPLIT4A;

        } else if ( strcmp(str, "discr_2split4B") == 0 ) {

            opts.discretization = fnft_manakov_discretization_2SPLIT4B;

        } else if ( strcmp(str, "discr_2split6B") == 0 ) {

            opts.discretization = fnft_manakov_discretization_2SPLIT6B;

        } else if ( strcmp(str, "discr_4split4A") == 0 ) {

            opts.discretization = fnft_manakov_discretization_4SPLIT4A;

        } else if ( strcmp(str, "discr_4split4B") == 0 ) {

            opts.discretization = fnft_manakov_discretization_4SPLIT4B;

        } else if ( strcmp(str, "discr_4split6B") == 0 ) {

            opts.discretization = fnft_manakov_discretization_4SPLIT6B;

		}
		else if (strcmp(str, "discr_FTES4_4A") == 0) {

			opts.discretization = fnft_manakov_discretization_FTES4_4A;

		}
		else if (strcmp(str, "discr_FTES4_4B") == 0) {

			opts.discretization = fnft_manakov_discretization_FTES4_4B;

        } else if ( strcmp(str, "discr_FTES4_suzuki") == 0 ) {

            opts.discretization = fnft_manakov_discretization_FTES4_suzuki;


        // Slow discretizations
        } else if ( strcmp(str, "discr_BO") == 0 ) {
            
            opts.discretization = fnft_manakov_discretization_BO;
            
        } else if ( strcmp(str, "discr_CF4_2") == 0 ) {
            
            opts.discretization = fnft_manakov_discretization_CF4_2;
            
        } else {
            snprintf(msg, sizeof msg, "%uth input has invalid value.", 
                (unsigned int)(k+1));
            goto on_error;
        }
    }
    
    /* Allocate memory */

    q1 = mxMalloc(D * sizeof(FNFT_COMPLEX));
    q2 = mxMalloc(D * sizeof(FNFT_COMPLEX));
    if (q1 == NULL || q2 == NULL) {
        snprintf(msg, sizeof msg, "Out of memory.");
        goto on_error;
    }
   
    if (skip_bound_states_flag == 0) {
        if (bound_states == NULL) {
            K = fnft_manakovv_max_K(D, &opts);
            if (K == 0) {
                snprintf(msg, sizeof msg, "Discretization does not support the chosen bound state localization method.");
                goto on_error;
            }
            bound_states = mxMalloc(K * sizeof(FNFT_COMPLEX));
        }
        if (bound_states == NULL) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }
    }


    if (skip_contspec_flag == 0) {
		if (opts.contspec_type == fnft_manakovv_cstype_AB)
			contspec = mxMalloc(3 * M * sizeof(COMPLEX));
		else if (opts.contspec_type == fnft_manakovv_cstype_REFLECTION_COEFFICIENT)
			contspec = mxMalloc(2 * M * sizeof(COMPLEX));
		else    // opts.contspec == fnft_manakovv_cstype_BOTH
			contspec = mxMalloc(5 * M * sizeof(COMPLEX));
/*        if (contspec == NULL) {
            snprintf(msg, sizeof msg, "Out of memory.");
            goto on_error;
        }*/ // TODO: check this
    }   

    
    /* Convert input */

    re_q1 = mxGetPr(prhs[0]);
    im_q1 = mxGetPi(prhs[0]);
    re_q2 = mxGetPr(prhs[1]);
    im_q2 = mxGetPi(prhs[1]);
    for (i=0; i<D; i++){
        q1[i] = re_q1[i] + I*im_q1[i];
        q2[i] = re_q2[i] + I*im_q2[i];
    }
    
    /* Call the C routine */
    
    ret_code = fnft_manakovv(D, q1, q2, T, M, contspec, XI, &K, bound_states,
        normconsts_or_residuals, kappa, &opts);
    if (ret_code != FNFT_SUCCESS) {
        snprintf(msg, sizeof msg, "fnft_manakovv failed (error code %i).",
                ret_code);
        goto on_error;
    }
    
    /* Allocate memory for the outputs */
    if (opts.contspec_type == fnft_manakovv_cstype_AB)
        plhs[0] = mxCreateDoubleMatrix(1, 3*M, mxCOMPLEX);
    else if (opts.contspec_type == fnft_manakovv_cstype_REFLECTION_COEFFICIENT)
        plhs[0] = mxCreateDoubleMatrix(1, 2*M, mxCOMPLEX);
    else    // opts.contspec == fnft_manakovv_cstype_BOTH
        plhs[0] = mxCreateDoubleMatrix(1, 5*M, mxCOMPLEX);
    
    /* Allocate memory for outputs and convert results */
    if (skip_contspec_flag == 0){
        csr = mxGetPr(plhs[0]);        // original: mxGetPr
        csi = mxGetPi(plhs[0]); // mxGetPi
        if (opts.contspec_type == fnft_manakovv_cstype_AB) {
            for (i=0; i<3*M; i++) {
                csr[i] = FNFT_CREAL(contspec[i]);
                csi[i] = FNFT_CIMAG(contspec[i]);
            }
        }
        else if (opts.contspec_type == fnft_manakovv_cstype_REFLECTION_COEFFICIENT) {
            for (i=0; i<2*M; i++) {
                csr[i] = FNFT_CREAL(contspec[i]);
                csi[i] = FNFT_CIMAG(contspec[i]);
            }
        }
        else{
            for (i=0; i<5*M; i++) {
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
/*
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

    }*/
    
    /* Free memory that is no longer needed */

    mxFree(q1);
    mxFree(q2);
    mxFree(contspec);
    mxFree(bound_states);
    mxFree(normconsts_or_residuals);
    return;
    
on_error:
    mexErrMsgTxt(msg);
}
