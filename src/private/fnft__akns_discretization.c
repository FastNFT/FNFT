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
 * Peter J. Prins (TU Delft) 2018, 2020.
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
INT fnft__akns_discretization_method_weights(COMPLEX ** qr_weights_ptr,
                                             COMPLEX ** eps_t_weights_ptr,
                                             akns_discretization_t const akns_discretization)
{
    COMPLEX *qr_weights = NULL, *eps_t_weights = NULL;
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
            qr_weights = malloc(1 * sizeof(COMPLEX));
            eps_t_weights = malloc(1 * sizeof(COMPLEX));
            if (qr_weights == NULL || eps_t_weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            qr_weights[0] = 1.0;
            eps_t_weights[0] = 1.0;
            
            break;
        case akns_discretization_CF4_2:
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
            qr_weights = malloc(4 * sizeof(COMPLEX));
            eps_t_weights = malloc(2 * sizeof(COMPLEX));
            if (qr_weights == NULL || eps_t_weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            REAL scl_factor = 1.0/SQRT(3.0);
            qr_weights[0] = 0.5+scl_factor;
            qr_weights[1] = 0.5-scl_factor;
            qr_weights[2] = 0.5-scl_factor;
            qr_weights[3] = 0.5+scl_factor;
            eps_t_weights[0] = 0.5;
            eps_t_weights[1] = 0.5;
            
            break;
        case akns_discretization_CF4_3:
            qr_weights = malloc(9 * sizeof(COMPLEX));
            eps_t_weights = malloc(3 * sizeof(COMPLEX));
            if (qr_weights == NULL || eps_t_weights == NULL) {
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
                    qr_weights[i*3+m] = 0.0;
                    for (n = 0; n < 3; n++){
                        qr_weights[i*3+m] = qr_weights[i*3+m] + (2*n+1)*misc_legendre_poly(n,xm[m])*f[i][n];
                    }
                    qr_weights[i*3+m] = qr_weights[i*3+m]*wm[m];
                }
            }
            for (m = 0; m < 3; m++){
                eps_t_weights[m] = 0.0;
                for (i = 0; i < 3; i++)
                    eps_t_weights[m] += qr_weights[3*m + i];
                for (i = 0; i < 3; i++)
                    qr_weights[3*m + i] /= eps_t_weights[m];
            }
            
            break;
        case akns_discretization_CF5_3:
            qr_weights = malloc(9 * sizeof(COMPLEX));
            eps_t_weights = malloc(3 * sizeof(COMPLEX));
            if (qr_weights == NULL || eps_t_weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            
            qr_weights[0] = ((145.0+37.0*SQRT(15.0))/900.0)+I*((5.0+3.0*SQRT(15.0))/300.0);
            qr_weights[1] = (-1.0/45.0)+I*(1.0/15.0);
            qr_weights[2] = ((145.0-37.0*SQRT(15.0))/900.0)+I*((5.0-3.0*SQRT(15.0))/300.0);
            qr_weights[3] = (-2.0/45.0)+I*(-SQRT(15.0)/50.0);
            qr_weights[4] = 22.0/45.0;
            qr_weights[5] = CONJ(qr_weights[3]);
            qr_weights[6] = CONJ(qr_weights[2]);
            qr_weights[7] = CONJ(qr_weights[1]);
            qr_weights[8] = CONJ(qr_weights[0]);

            for (m = 0; m < 3; m++){
                eps_t_weights[m] = 0.0;
                for (i = 0; i < 3; i++)
                    eps_t_weights[m] += qr_weights[3*m + i];
                for (i = 0; i < 3; i++)
                    qr_weights[3*m + i] /= eps_t_weights[m];
            }
            
            break;
        case akns_discretization_CF6_4:
            qr_weights = malloc(12 * sizeof(COMPLEX));
            eps_t_weights = malloc(4 * sizeof(COMPLEX));
            if (qr_weights == NULL || eps_t_weights == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            
            qr_weights[0] = (0.245985577298764 + 0.038734389227165*I);
            qr_weights[1] = (-0.046806149832549 + 0.012442141491185*I);
            qr_weights[2] = (0.010894359342569 - 0.004575808769067*I);
            qr_weights[3] = (0.062868370946917 - 0.048761268117765*I);
            qr_weights[4] = (0.269028372054771 - 0.012442141491185*I);
            qr_weights[5] = (-0.041970529810473 + 0.014602687659668*I);
            qr_weights[6] = (-0.041970529810473 + 0.014602687659668*I);
            qr_weights[7] = (0.269028372054771 - 0.012442141491185*I);
            qr_weights[8] = (0.062868370946917 - 0.048761268117765*I);
            qr_weights[9] = (0.010894359342569 - 0.004575808769067*I);
            qr_weights[10] = (-0.046806149832549 + 0.012442141491185*I);
            qr_weights[11] = (0.245985577298764 + 0.038734389227165*I);

            for (m = 0; m < 4; m++){
                eps_t_weights[m] = 0.0;
                for (i = 0; i < 3; i++)
                    eps_t_weights[m] += qr_weights[3*m + i];
                for (i = 0; i < 3; i++)
                    qr_weights[3*m + i] /= eps_t_weights[m];
            }
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(akns_discretization);
            goto  leave_fun;
    }
    
leave_fun:
    if (ret_code!=SUCCESS) {
        free(qr_weights);
        free(eps_t_weights);
    } else {
        *qr_weights_ptr = qr_weights;
        *eps_t_weights_ptr = eps_t_weights;
    }
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

    COMPLEX *qr_weights = NULL;
    COMPLEX *eps_t_weights = NULL;

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

            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=2, i+= nskip_per_step) {
                q_preprocessed[isub] = qr_weights[0]*q_1[i] + qr_weights[1]*q_2[i];
                q_preprocessed[isub+1] = qr_weights[2]*q_1[i] + qr_weights[3]*q_2[i];
                r_preprocessed[isub] = qr_weights[0]*r_from_q[0](q_1[i])
                                       + qr_weights[1]*r_from_q[0](q_2[i]);
                r_preprocessed[isub+1] = qr_weights[2]*r_from_q[0](q_1[i])
                                         + qr_weights[3]*r_from_q[0](q_2[i]);
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

            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=3, i+=nskip_per_step) {
                q_preprocessed[isub] = qr_weights[0]*q_1[i] + qr_weights[1]*q[i] + qr_weights[2]*q_3[i];
                q_preprocessed[isub+1] = qr_weights[3]*q_1[i] + qr_weights[4]*q[i] + qr_weights[5]*q_3[i];
                q_preprocessed[isub+2] = qr_weights[6]*q_1[i] + qr_weights[7]*q[i] + qr_weights[8]*q_3[i];
                r_preprocessed[isub]   = qr_weights[0]*r_from_q[0](q_1[i])
                                         + qr_weights[1]*r_from_q[0](q[i])
                                         + qr_weights[2]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+1] = qr_weights[3]*r_from_q[0](q_1[i])
                                         + qr_weights[4]*r_from_q[0](q[i])
                                         + qr_weights[5]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+2] = qr_weights[6]*r_from_q[0](q_1[i])
                                         + qr_weights[7]*r_from_q[0](q[i])
                                         + qr_weights[8]*r_from_q[0](q_3[i]);
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

            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=3, i+=nskip_per_step) {
                q_preprocessed[isub]   = qr_weights[0]*q_1[i] + qr_weights[1]*q[i] + qr_weights[2]*q_3[i];
                q_preprocessed[isub+1] = qr_weights[3]*q_1[i] + qr_weights[4]*q[i] + qr_weights[5]*q_3[i];
                q_preprocessed[isub+2] = qr_weights[6]*q_1[i] + qr_weights[7]*q[i] + qr_weights[8]*q_3[i];
                r_preprocessed[isub]   = qr_weights[0]*r_from_q[0](q_1[i])
                                         + qr_weights[1]*r_from_q[0](q[i])
                                         + qr_weights[2]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+1] = qr_weights[3]*r_from_q[0](q_1[i])
                                         + qr_weights[4]*r_from_q[0](q[i])
                                         + qr_weights[5]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+2] = qr_weights[6]*r_from_q[0](q_1[i])
                                         + qr_weights[7]*r_from_q[0](q[i])
                                         + qr_weights[8]*r_from_q[0](q_3[i]);
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

            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            for (isub=0, i=0; isub<D_effective; isub+=4, i+=nskip_per_step) {
                q_preprocessed[isub]   = qr_weights[0]*q_1[i] + qr_weights[1]*q[i]  + qr_weights[2]*q_3[i];
                q_preprocessed[isub+1] = qr_weights[3]*q_1[i] + qr_weights[4]*q[i]  + qr_weights[5]*q_3[i];
                q_preprocessed[isub+2] = qr_weights[6]*q_1[i] + qr_weights[7]*q[i]  + qr_weights[8]*q_3[i];
                q_preprocessed[isub+3] = qr_weights[9]*q_1[i] + qr_weights[10]*q[i] + qr_weights[11]*q_3[i];
                r_preprocessed[isub]   = qr_weights[0]*r_from_q[0](q_1[i])
                                         + qr_weights[1]*r_from_q[0](q[i])
                                         + qr_weights[2]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+1] = qr_weights[3]*r_from_q[0](q_1[i])
                                         + qr_weights[4]*r_from_q[0](q[i])
                                         + qr_weights[5]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+2] = qr_weights[6]*r_from_q[0](q_1[i])
                                         + qr_weights[7]*r_from_q[0](q[i])
                                         + qr_weights[8]*r_from_q[0](q_3[i]);
                r_preprocessed[isub+3] = qr_weights[9]*r_from_q[0](q_1[i])
                                         + qr_weights[10]*r_from_q[0](q[i])
                                         + qr_weights[11]*r_from_q[0](q_3[i]);
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
        free(qr_weights);
        free(eps_t_weights);
        return ret_code;
}

/**
 * This routine returns the change of basis matrix from the S basis to the basis of the discretization.
 */
INT fnft__akns_discretization_change_of_basis_matrix_to_S(COMPLEX * const T,
                                                          COMPLEX const xi,
                                                          UINT const  derivative_flag, // 0- > 2x2, 1->4x4
                                                          REAL const eps_t,
                                                          akns_discretization_t const akns_discretization,
                                                          UINT vanilla_flag,
                                                          akns_pde_t const PDE)
{
    INT ret_code = SUCCESS;

    switch (PDE) {
        case akns_pde_KdV:
            if(xi==0){
                ret_code = E_DIV_BY_ZERO;
                CHECK_RETCODE(ret_code, release_mem);
            }
            switch (akns_discretization) {
                // Fill the first two columns of T for the discretization basis
                case akns_discretization_2SPLIT1A:
                case akns_discretization_2SPLIT1B:
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
                case akns_discretization_4SPLIT4A:
                case akns_discretization_4SPLIT4B:
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                case akns_discretization_ES4:
                case akns_discretization_TES4:
                    if (vanilla_flag) {
                        // T from AKNS basis to S basis for KdV:

                        T[0] = -I * 0.5 / xi;   //T_11;
                        T[1] = 0.0;             //T_12;
                        if(!derivative_flag){
                            T[2]  = -T[0];      //T_21;
                            T[3]  = 1.0;        //T_22;
                        } else {
                            T[4]  = -T[0];      //T_21;
                            T[5]  = 1.0;        //T_22;

                            T[8]  = -T[0] / xi; //T_31 = d/dxi T_11
                            T[9]  = 0.0;        //T_32 = d/dxi T_12
                            T[12] = -T[8];      //T_41 = d/dxi T_21
                            T[13] = 0.0;        //T_42 = d/dxi T_22
                        }
                    } else {
                        // T from flipped AKNS basis to S basis for KdV:

                        T[0] = 1.0;             //T_11;
                        T[1] = -I * 0.5 / xi;   //T_12;
                        if(!derivative_flag){
                            T[2]  = 0.0;        //T_21;
                            T[3]  = -T[1];      //T_22;
                        } else {
                            T[4]  = 0.0;        //T_21;
                            T[5]  = -T[1];      //T_22;

                            T[8]  = 0.0;        //T_31 = d/dxi T_11
                            T[9]  = -T[1] / xi; //T_32 = d/dxi T_12
                            T[12] = 0.0;        //T_41 = d/dxi T_21
                            T[13] = -T[9];      //T_42 = d/dxi T_22
                        }
                    }
                    break;

                case akns_discretization_2SPLIT2_MODAL:
                case akns_discretization_2SPLIT2A:
                    if (vanilla_flag) {
                        //T from modified AKNS basis to S basis for KdV. This modification reduces the degree of the polynomial scattering matrix by 1 per step.

                        T[0] = CEXP(-0.5*I*xi*eps_t) * -I * 0.5 / xi; //T_11;
                        T[1] = 0.0;                                   //T_12;
                        if(!derivative_flag){
                            T[2]  = -T[0];                            //T_21;
                            T[3]  = CEXP(+0.5*I*xi*eps_t);            //T_22;
                        } else {
                            T[4]  = -T[0];                            //T_21;
                            T[5]  = CEXP(+0.5*I*xi*eps_t);            //T_22;

                            T[8]  = -T[0] * (1.0/xi + 0.5*I*eps_t);   //T_31 = d/dxi T_11
                            T[9]  = 0.0;                              //T_32 = d/dxi T_12
                            T[12] = -T[8];                            //T_41 = d/dxi T_21
                            T[13] = 0.5 * I * eps_t * T[5];           //T_42 = d/dxi T_22
                        }
                    } else {
                        //T from flipped modified AKNS basis to S basis for KdV. This modification reduces the degree of the polynomial scattering matrix by 1 per step.

                        T[0] = CEXP(-0.5*I*xi*eps_t);                 //T_11;
                        T[1] = CEXP( 0.5*I*xi*eps_t) * -I * 0.5 / xi; //T_12;
                        if(!derivative_flag){
                            T[2]  = 0.0;                              //T_21;
                            T[3]  = -T[1];                            //T_22;
                        } else {
                            T[4]  = 0.0;                              //T_21;
                            T[5]  = -T[1];                            //T_22;

                            T[8]  = -0.5 * I * eps_t*T[0];            //T_31 = d/dxi T_11
                            T[9]  = T[1] * (-1.0/xi + 0.5*I*eps_t);   //T_32 = d/dxi T_12
                            T[12] = 0.0;                              //T_41 = d/dxi T_21
                            T[13] = -T[9];                            //T_42 = d/dxi T_22

                        }
                    }
                    break;
                default:
                    ret_code =  E_INVALID_ARGUMENT(akns_discretization);
                    CHECK_RETCODE(ret_code, release_mem);
            }
            break;
            
        case akns_pde_NSE:
            switch (akns_discretization) {
                case akns_discretization_2SPLIT2_MODAL:
                case akns_discretization_2SPLIT2A:
                    // These two discretizations make use of a modification of
                    // the AKNS basis (cf. akns_pde_KdV), which reduces the order
                    // of the required polynomials. However, they must be kept
                    // in this basis to use polynomial rootfinding for the bound
                    // states. In the end this is corrected using
                    // fnft__nse_discretization_phase_factor_rho() and /or
                    // fnft__nse_discretization_phase_factor_b(), which in essence
                    // calculalate the change of basis to E basis.
                case akns_discretization_2SPLIT1A:
                case akns_discretization_2SPLIT1B:
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
                case akns_discretization_4SPLIT4A:
                case akns_discretization_4SPLIT4B:
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                case akns_discretization_ES4:
                case akns_discretization_TES4:
                    // The AKNS basis is already the S-basis for these discretizations, return an identity matrix
                    T[0] = 1.0;            //T_11;
                    T[1] = 0.0;            //T_12;
                    if(!derivative_flag){
                        T[2]  = 0.0;       //T_21;
                        T[3]  = 1.0;       //T_22;
                    } else {
                        T[4]  = 0.0;       //T_21;
                        T[5]  = 1.0;       //T_22;

                        T[8]  = 0.0;       //T_31 = d/dxi T_11
                        T[9]  = 0.0;       //T_32 = d/dxi T_12
                        T[12] = 0.0;       //T_41 = d/dxi T_21
                        T[13] = 0.0;       //T_42 = d/dxi T_22
                    }
                    break;

                default:
                    ret_code =  E_INVALID_ARGUMENT(akns_discretization);
                    CHECK_RETCODE(ret_code, release_mem);
            }
            break;

        default:
            ret_code =  E_INVALID_ARGUMENT(PDE);
            CHECK_RETCODE(ret_code, release_mem);
    }

    if (derivative_flag) {
        // Fill the last two columns of T from the first two
        T[2]  = 0.0;  // T_13 = 0
        T[3]  = 0.0;  // T_14 = 0
        T[6]  = 0.0;  // T_23 = 0
        T[7]  = 0.0;  // T_24 = 0
        T[10] = T[0]; // T_33 = T_11
        T[11] = T[1]; // T_34 = T_12
        T[14] = T[4]; // T_43 = T_21
        T[15] = T[5]; // T_44 = T_22
    }

    release_mem:
        return ret_code;
}

/**
 * This routine returns the change of basis matrix from the S basis to the basis of the discretization.
 */
INT fnft__akns_discretization_change_of_basis_matrix_from_S(COMPLEX * const T,
                                                            COMPLEX const xi,
                                                            UINT const  derivative_flag, // 0- > 2x2, 1->4x4
                                                            REAL const eps_t,
                                                            akns_discretization_t const akns_discretization,
                                                            UINT vanilla_flag,
                                                            akns_pde_t const PDE)
{
    INT ret_code = SUCCESS;

    switch (PDE) {
        case akns_pde_KdV:
            if(xi==0){
                ret_code = E_DIV_BY_ZERO;
                CHECK_RETCODE(ret_code, release_mem);
            }
            switch (akns_discretization) {
                // Fill the first two columns of T for the discretization basis
                case akns_discretization_2SPLIT1A:
                case akns_discretization_2SPLIT1B:
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
                case akns_discretization_4SPLIT4A:
                case akns_discretization_4SPLIT4B:
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                case akns_discretization_ES4:
                case akns_discretization_TES4:
                    if (vanilla_flag) {
                        // T from S basis for KdV to AKNS basis:

                        T[0] = 2.0 * I * xi;   //T_11;
                        T[1] = 0.0;            //T_12;
                        if(!derivative_flag){
                            T[2]  = 1.0;       //T_21;
                            T[3]  = 1.0;       //T_22;
                        } else {
                            T[4]  = 1.0;       //T_21;
                            T[5]  = 1.0;       //T_22;

                            T[8]  = 2.0 * I;   //T_31 = d/dxi T_11
                            T[9]  = 0.0;       //T_32 = d/dxi T_12
                            T[12] = 0.0;       //T_41 = d/dxi T_21
                            T[13] = 0.0;       //T_42 = d/dxi T_22
                        }
                    } else {
                        // T from S basis for KdV to flipped AKNS basis:

                        T[0] = 1.0;                //T_11;
                        T[1] = 1.0;                //T_12;
                        if(!derivative_flag){
                            T[2]  = 0.0;           //T_21;
                            T[3]  = -2.0 * I * xi; //T_22;
                        } else {
                            T[4]  = 0.0;           //T_21;
                            T[5]  = -2.0 * I * xi; //T_22;

                            T[8]  = 0.0;           //T_31 = d/dxi T_11
                            T[9]  = 0.0;           //T_32 = d/dxi T_12
                            T[12] = 0.0;           //T_41 = d/dxi T_21
                            T[13] = -2.0 * I;      //T_42 = d/dxi T_22
                        }
                    }
                    break;

                case akns_discretization_2SPLIT2_MODAL:
                case akns_discretization_2SPLIT2A:
                    if (vanilla_flag) {
                        // T from S basis for KdV to modified AKNS basis. This modification reduces the degree of the polynomial scattering matrix by 1 per step.

                        T[0] = CEXP(+0.5*I*xi*eps_t) * 2.0 * I * xi;  //T_11;
                        T[1] = 0.0;                                   //T_12;
                        if(!derivative_flag){
                            T[2]  = CEXP(-0.5*I*xi*eps_t);            //T_21;
                            T[3]  = T[2];                             //T_22;
                        } else {
                            T[4]  = CEXP(-0.5*I*xi*eps_t);            //T_21;
                            T[5]  = T[4];                             //T_22;

                            T[8]  = T[0] * (1.0/xi + 0.5*I*eps_t);    //T_31 = d/dxi T_11
                            T[9]  = 0.0;                              //T_32 = d/dxi T_12
                            T[12] = T[4] * -0.5 * I * eps_t;          //T_41 = d/dxi T_21
                            T[13] = T[12];                            //T_42 = d/dxi T_22
                        }
                    } else {
                        // T from S basis for KdV to flipped modified AKNS basis. This modification reduces the degree of the polynomial scattering matrix by 1 per step.
                        T[0] = CEXP(0.5*I*xi*eps_t);                       //T_11;
                        T[1] = T[0];                                       //T_12;

                        if(!derivative_flag){
                            T[2]  = 0.0;                                   //T_21;
                            T[3]  = CEXP(-0.5*I*xi*eps_t) * -2.0 * I * xi; //T_22;
                        } else {
                            T[4]  = 0.0;                                   //T_21;
                            T[5]  = CEXP(-0.5*I*xi*eps_t) * -2.0 * I * xi; //T_22;

                            T[8]  = T[0] * 0.5 * I * eps_t;                //T_31 = d/dxi T_11
                            T[9]  = T[8];                                  //T_32 = d/dxi T_12
                            T[12] = 0.0;                                   //T_41 = d/dxi T_21
                            T[13] = T[5] * (1.0/xi - 0.5*I*eps_t);         //T_42 = d/dxi T_22
                        }
                    }
                    break;
                default:
                    ret_code =  E_INVALID_ARGUMENT(akns_discretization);
                    CHECK_RETCODE(ret_code, release_mem);
            }
            break;

        case akns_pde_NSE:
            switch (akns_discretization) {
                case akns_discretization_2SPLIT1A:
                case akns_discretization_2SPLIT1B:
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
                case akns_discretization_4SPLIT4A:
                case akns_discretization_4SPLIT4B:
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                case akns_discretization_ES4:
                case akns_discretization_TES4:
                    // The AKNS basis is already the S-basis for these discretizations, return an identity matrix
                    T[0] = 1.0;            //T_11;
                    T[1] = 0.0;            //T_12;
                    if(!derivative_flag){
                        T[2]  = 0.0;       //T_21;
                        T[3]  = 1.0;       //T_22;
                    } else {
                        T[4]  = 0.0;       //T_21;
                        T[5]  = 1.0;       //T_22;

                        T[8]  = 0.0;       //T_31 = d/dxi T_11
                        T[9]  = 0.0;       //T_32 = d/dxi T_12
                        T[12] = 0.0;       //T_41 = d/dxi T_21
                        T[13] = 0.0;       //T_42 = d/dxi T_22
                    }
                    break;

                case akns_discretization_2SPLIT2_MODAL:
                case akns_discretization_2SPLIT2A:
                    // T from S basis for NSE to modified AKNS basis. This modification reduces the degree of the polynomial scattering matrix by 1 per step.
                    // TODO: These discretizations require a modification that is currently contracted in the 'boundary coefficionts'
                    T[0] = 1.0;            //T_11;
                    T[1] = 0.0;            //T_12;
                    if(!derivative_flag){
                        T[2]  = 0.0;       //T_21;
                        T[3]  = 1.0;       //T_22;
                    } else {
                        T[4]  = 0.0;       //T_21;
                        T[5]  = 1.0;       //T_22;

                        T[8]  = 0.0;       //T_31 = d/dxi T_11
                        T[9]  = 0.0;       //T_32 = d/dxi T_12
                        T[12] = 0.0;       //T_41 = d/dxi T_21
                        T[13] = 0.0;       //T_42 = d/dxi T_22
                    }
                    break;
                default:
                    ret_code =  E_INVALID_ARGUMENT(akns_discretization);
                    CHECK_RETCODE(ret_code, release_mem);
            }
            break;

        default:
            ret_code =  E_INVALID_ARGUMENT(PDE);
            CHECK_RETCODE(ret_code, release_mem);
    }

    if (derivative_flag) {
        // Fill the last two columns of T from the first two
        T[2]  = 0.0;  // T_13 = 0
        T[3]  = 0.0;  // T_14 = 0
        T[6]  = 0.0;  // T_23 = 0
        T[7]  = 0.0;  // T_24 = 0
        T[10] = T[0]; // T_33 = T_11
        T[11] = T[1]; // T_34 = T_12
        T[14] = T[4]; // T_43 = T_21
        T[15] = T[5]; // T_44 = T_22
    }

    release_mem:
        return ret_code;
}
