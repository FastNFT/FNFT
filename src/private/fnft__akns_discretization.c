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
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J. Prins (TU Delft) 2018.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__akns_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 */
UINT fnft__akns_discretization_degree(akns_discretization_t
        discretization)
{
    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A: // With change of base trick
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT2_MODAL:
            return 1;
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_4SPLIT4B:
            return 2;
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
            return 3;
        case akns_discretization_2SPLIT4A:
        case akns_discretization_4SPLIT4A:
            return 4;
        case akns_discretization_2SPLIT6B:
            return 6;
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT8B:
            return 12;
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
            return 15;
        case akns_discretization_2SPLIT8A:
            return 24;
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
            return 105;
            
        default: // Unknown discretization
            return 0;
    }
}

/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__akns_discretization_boundary_coeff(akns_discretization_t discretization)
{
    
    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
        case akns_discretization_BO:
        case akns_discretization_CF4_2:
        case akns_discretization_CF4_3:
        case akns_discretization_CF5_3:
        case akns_discretization_CF6_4:
        case akns_discretization_ES4:
        case akns_discretization_TES4:
            return 0.5;
            
        default: // Unknown discretization
            return NAN;
    }
}

/**
 * This routine returns the scaling for effective number of samples based on the discretization.
 */
UINT fnft__akns_discretization_upsampling_factor(akns_discretization_t discretization)
{
    
    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
        case akns_discretization_BO:
            return 1;
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_CF4_2:
            return 2;
        case akns_discretization_CF4_3:
        case akns_discretization_CF5_3:
        case akns_discretization_ES4:
        case akns_discretization_TES4:
            return 3;
        case akns_discretization_CF6_4:
            return 4;
            
        default: // Unknown discretization
            return 0;
    }
}
/**
 * This routine returns the order of the method based on the discretization.
 */
UINT fnft__akns_discretization_method_order(akns_discretization_t discretization)
{
    
    switch (discretization) {
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
        case akns_discretization_BO:
            return 2;
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_CF4_2:
        case akns_discretization_CF4_3:
        case akns_discretization_ES4:
        case akns_discretization_TES4:
            return 4;
        case akns_discretization_CF5_3:
            return 5;
        case akns_discretization_CF6_4:
            return 6;
            
        default: // Unknown discretization
            return 0;
    }
}

/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization.
 */
INT fnft__akns_discretization_lambda_to_z(const UINT n, const REAL eps_t,
        COMPLEX * const vals, akns_discretization_t discretization)
{
    REAL degree1step;
    UINT i, upsampling_factor;
    upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = akns_discretization_degree(discretization);
    if (degree1step == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = degree1step * upsampling_factor;
    for (i = 0; i < n; i++)
        vals[i] = CEXP(2*I*vals[i]*eps_t/degree1step);
    return SUCCESS;
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization.
 */
INT fnft__akns_discretization_z_to_lambda(const UINT n, const REAL eps_t,
        COMPLEX * const vals, akns_discretization_t discretization)
{
    REAL degree1step;
    UINT i, upsampling_factor;
    upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = akns_discretization_degree(discretization);
    if (degree1step == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = degree1step * upsampling_factor;
    for (i = 0; i < n; i++)
        vals[i] = CLOG(vals[i])/(2*I*eps_t/degree1step);
    return SUCCESS;
}

/**
 * This routine computes various weights required by some methods
 * based on the discretization.
 */
INT fnft__akns_discretization_method_weights(COMPLEX **weights_ptr,
        akns_discretization_t akns_discretization)
{
    COMPLEX *weights = NULL;
    INT ret_code = SUCCESS;
    
    
    switch (akns_discretization) {
        
        case akns_discretization_BO: // Bofetta-Osborne scheme
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
            weights = malloc(1 * sizeof(COMPLEX));
            if (weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            weights[0] = 1.0;
            
            
            break;
        case akns_discretization_CF4_2:
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
            weights = malloc(4 * sizeof(COMPLEX));
            if (weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            REAL scl_factor = SQRT(3.0)/6.0;
            weights[0] = 0.25+scl_factor;
            weights[1] = 0.25-scl_factor;
            weights[2] = 0.25-scl_factor;
            weights[3] = 0.25+scl_factor;
            
            break;
        case akns_discretization_CF4_3:
            weights = malloc(9 * sizeof(COMPLEX));
            if (weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            REAL f[3][3] = {{0}};
            f[0][0] = 11.0/40.0;
            f[0][1] = 20.0/87.0;
            f[0][2] = 7.0/50.0;
            f[1][0] = 9.0/20.0;
            f[1][1] = 0.0;
            f[1][2] = -7.0/25.0;
            f[2][0] = 11.0/40.0;
            f[2][1] = -20.0/87.0;
            f[2][2] = 7.0/50.0;
            
            REAL wm[3] = {5.0/18.0,4.0/9.0,5.0/18.0};
            REAL xm[3] = {2.0*SQRT(3.0/20.0),0.0,-2.0*SQRT(3.0/20.0)};//misc_legendre_poly is defined on [-1,1]
            UINT m, n, i;
            for (m = 0; m < 3; m++){
                for (i = 0; i < 3; i++){
                    weights[i*3+m] = 0.0;
                    for (n = 0; n < 3; n++){
                        weights[i*3+m] = weights[i*3+m] + (2*n+1)*misc_legendre_poly(n,xm[m])*f[i][n];
                    }
                    weights[i*3+m] = weights[i*3+m]*wm[m];
                }
            }
            
            break;
        case akns_discretization_CF5_3:
            weights = malloc(9 * sizeof(COMPLEX));
            if (weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            
            weights[0] = ((145.0+37.0*SQRT(15.0))/900.0)+I*((5.0+3.0*SQRT(15.0))/300.0);
            weights[1] = (-1.0/45.0)+I*(1.0/15.0);
            weights[2] = ((145.0-37.0*SQRT(15.0))/900.0)+I*((5.0-3.0*SQRT(15.0))/300.0);
            weights[3] = (-2.0/45.0)+I*(-SQRT(15.0)/50.0);
            weights[4] = 22.0/45.0;
            weights[5] = CONJ(weights[3]);
            weights[6] = CONJ(weights[2]);
            weights[7] = CONJ(weights[1]);
            weights[8] = CONJ(weights[0]);
            
            break;
        case akns_discretization_CF6_4:
            weights = malloc(12 * sizeof(COMPLEX));
            if (weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            
            weights[0] = (0.245985577298764 + 0.038734389227165*I);
            weights[1] = (-0.046806149832549 + 0.012442141491185*I);
            weights[2] = (0.010894359342569 - 0.004575808769067*I);
            weights[3] = (0.062868370946917 - 0.048761268117765*I);
            weights[4] = (0.269028372054771 - 0.012442141491185*I);
            weights[5] = (-0.041970529810473 + 0.014602687659668*I);
            weights[6] = (-0.041970529810473 + 0.014602687659668*I);
            weights[7] = (0.269028372054771 - 0.012442141491185*I);
            weights[8] = (0.062868370946917 - 0.048761268117765*I);
            weights[9] = (0.010894359342569 - 0.004575808769067*I);
            weights[10] = (-0.046806149832549 + 0.012442141491185*I);
            weights[11] = (0.245985577298764 + 0.038734389227165*I);
            
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto  leave_fun;
    }
    
    
    *weights_ptr = weights;
    
    leave_fun:
        return ret_code;
}

/**
 * This routine preprocess the signal by resampling and subsampling based on the discretization.
 * The preprocessing is necessary for higher-order methods.
 */
INT fnft__akns_discretization_preprocess_signal(UINT const D,
                                                COMPLEX const * const q,
                                                COMPLEX (*r_from_q[3])(COMPLEX const),
                                                REAL const eps_t,
                                                UINT * const Dsub_ptr,
                                                COMPLEX **q_preprocessed_ptr,
                                                COMPLEX **r_preprocessed_ptr,
                                                UINT * const first_last_index,
                                                akns_discretization_t discretization)
{

    UINT i, D_effective, isub;
    INT ret_code = SUCCESS;
    COMPLEX *q_1 = NULL;
    COMPLEX *q_2 = NULL;
    COMPLEX *q_3 = NULL;

    COMPLEX *weights = NULL;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (Dsub_ptr == NULL)
        return E_INVALID_ARGUMENT(Dsub_ptr);
    if (q_preprocessed_ptr == NULL)
        return E_INVALID_ARGUMENT(q_preprocessed_ptr);
    if (r_preprocessed_ptr == NULL)
        return E_INVALID_ARGUMENT(r_preprocessed_ptr);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(esp_t);
    if (first_last_index == NULL)
        return E_INVALID_ARGUMENT(first_last_index);



    // Determine number of samples after downsampling, Dsub
    UINT Dsub = *Dsub_ptr; // desired Dsub
    if (Dsub < 2)
        Dsub = 2;
    if (Dsub > D)
        Dsub = D;
    const UINT nskip_per_step = (UINT) ROUND((REAL)D / Dsub);
    Dsub = (UINT) ROUND((REAL)D / nskip_per_step); // actual Dsub

    UINT upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0){
        ret_code =  E_INVALID_ARGUMENT(discretization);
        goto release_mem;
    }
    D_effective = Dsub * upsampling_factor;
    COMPLEX * const q_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    COMPLEX * const r_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    if (q_preprocessed == NULL || r_preprocessed == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    switch (discretization) {

        case akns_discretization_BO: // Bofetta-Osborne scheme
        case akns_discretization_2SPLIT1A:
        case akns_discretization_2SPLIT1B:
        case akns_discretization_2SPLIT2A:
        case akns_discretization_2SPLIT2B:
        case akns_discretization_2SPLIT2S:
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
        case akns_discretization_2SPLIT3S:
        case akns_discretization_2SPLIT4A:
        case akns_discretization_2SPLIT4B:
        case akns_discretization_2SPLIT5A:
        case akns_discretization_2SPLIT5B:
        case akns_discretization_2SPLIT6A:
        case akns_discretization_2SPLIT6B:
        case akns_discretization_2SPLIT7A:
        case akns_discretization_2SPLIT7B:
        case akns_discretization_2SPLIT8A:
        case akns_discretization_2SPLIT8B:
        case akns_discretization_2SPLIT2_MODAL:
            for (isub=0, i=0; isub<D_effective; isub++, i += nskip_per_step) {
                q_preprocessed[isub] = q[i];
                r_preprocessed[isub] = r_from_q[0](q[i]);
            }
            break;
        case akns_discretization_CF4_2:
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_2 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_2 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }
            REAL scl_factor = SQRT(3.0)/6.0;
            ret_code = misc_resample(D, eps_t, q, -eps_t*scl_factor*nskip_per_step, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*scl_factor*nskip_per_step, q_2);
            CHECK_RETCODE(ret_code, release_mem);


            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=2, i+= nskip_per_step) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q_2[i];
                q_preprocessed[isub+1] = weights[2]*q_1[i] + weights[3]*q_2[i];
            }
            for (isub = 0; isub<D_effective; isub++) {
                r_preprocessed[isub] = r_from_q[0](q_preprocessed[isub]);
            }

            break;
        case akns_discretization_CF4_3:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_3 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_3 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }

            ret_code = misc_resample(D, eps_t, q, -eps_t*SQRT(3.0/20.0)*nskip_per_step, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*SQRT(3.0/20.0)*nskip_per_step, q_3);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=3, i+=nskip_per_step) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q[i] + weights[2]*q_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i] + weights[4]*q[i] + weights[5]*q_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
            }
            for (i=0; i<D_effective; i++) {
                r_preprocessed[i] = r_from_q[0](q_preprocessed[i]);
            }
            break;
        case akns_discretization_CF5_3:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_3 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_3 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }

            ret_code = misc_resample(D, eps_t, q, -eps_t*SQRT(15.0)/10.0*nskip_per_step, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*SQRT(15.0)/10.0*nskip_per_step, q_3);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);


            for (isub=0, i=0; isub<D_effective; isub+=3, i+=nskip_per_step) {
                q_preprocessed[isub]   = weights[0]*q_1[i] + weights[1]*q[i] + weights[2]*q_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i] + weights[4]*q[i] + weights[5]*q_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
                r_preprocessed[isub]   = weights[0]*r_from_q[0](q_1[i])
                                         + weights[1]*r_from_q[0](q[i])
                                         + weights[2]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+1] = weights[3]*r_from_q[0](q_1[i])
                                         + weights[4]*r_from_q[0](q[i])
                                         + weights[5]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+2] = weights[6]*r_from_q[0](q_1[i])
                                         + weights[7]*r_from_q[0](q[i])
                                         + weights[8]*r_from_q[0](q_3[i]);
            }


            break;
        case akns_discretization_CF6_4:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_3 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_3 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }

            ret_code = misc_resample(D, eps_t, q, -eps_t*nskip_per_step*SQRT(15.0)/10.0, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*nskip_per_step*SQRT(15.0)/10.0, q_3);
            CHECK_RETCODE(ret_code, release_mem);

            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=4, i+=nskip_per_step) {
                q_preprocessed[isub]   = weights[0]*q_1[i] + weights[1]*q[i]  + weights[2]*q_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i] + weights[4]*q[i]  + weights[5]*q_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i]  + weights[8]*q_3[i];
                q_preprocessed[isub+3] = weights[9]*q_1[i] + weights[10]*q[i] + weights[11]*q_3[i];
                r_preprocessed[isub]   = weights[0]*r_from_q[0](q_1[i])
                                         + weights[1]*r_from_q[0](q[i])
                                         + weights[2]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+1] = weights[3]*r_from_q[0](q_1[i])
                                         + weights[4]*r_from_q[0](q[i])
                                         + weights[5]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+2] = weights[6]*r_from_q[0](q_1[i])
                                         + weights[7]*r_from_q[0](q[i])
                                         + weights[8]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+3] = weights[9]*r_from_q[0](q_1[i])
                                         + weights[10]*r_from_q[0](q[i])
                                         + weights[11]*r_from_q[0](q_3[i]);
            }
            break;
        case akns_discretization_ES4:
        case akns_discretization_TES4:
            for (isub=0, i=0; isub<D_effective; isub+=3, i+=nskip_per_step) {
                q_preprocessed[isub] = q[i];
            }

            REAL eps_t_sub = eps_t*nskip_per_step;
            REAL eps_t_sub_2 = POW(eps_t_sub,2);
            q_preprocessed[1] = (q_preprocessed[3]-0)/(2*eps_t_sub);
            q_preprocessed[2] = (q_preprocessed[3]-2*q_preprocessed[0]+0)/eps_t_sub_2;
            q_preprocessed[D_effective-2] = (0-q_preprocessed[D_effective-6])/(2*eps_t_sub);
            q_preprocessed[D_effective-1] = (0-2*q_preprocessed[D_effective-3]+q_preprocessed[D_effective-6])/eps_t_sub_2;

            for (isub=3; isub<D_effective-3; isub+=3) {
                q_preprocessed[isub+1] = (q_preprocessed[isub+3]-q_preprocessed[isub-3])/(2*eps_t_sub);
                q_preprocessed[isub+2] = (q_preprocessed[isub+3]-2*q_preprocessed[isub]+q_preprocessed[isub-3])/eps_t_sub_2;
            }
            for (isub=0; isub<D_effective; isub+=3) {
                for (i=0; i<3; i++) {
                    r_preprocessed[isub+i] = r_from_q[i](q_preprocessed[isub+i]);
                }
            }
            break;
        default: // Unknown discretization

            ret_code = E_INVALID_ARGUMENT(discretization);
            goto  release_mem;
    }

    // Original index of the first and last sample in qsub
    first_last_index[0] = 0;
    first_last_index[1] = (Dsub-1)*nskip_per_step;
    *q_preprocessed_ptr = q_preprocessed;
    *r_preprocessed_ptr = r_preprocessed;
    *Dsub_ptr = Dsub;

    release_mem:
        free(q_1);
        free(q_2);
        free(q_3);
        free(weights);
        return ret_code;
}
