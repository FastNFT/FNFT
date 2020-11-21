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
* Sander Wahls (TU Delft) 2017.
* Shrinivas Chimmalgi (TU Delft) 2017-2020.
* Peter J Prins (TU Delft) 2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__kdv_discretization.h"
#include "fnft__akns_discretization.h"

/**
 * Returns the max degree of the polynomials in a single scattering
 * matrix or zero if the discretization is unknown.
 */
UINT fnft__kdv_discretization_degree(kdv_discretization_t
        kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code;
    UINT degree1step = 0;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    degree1step = akns_discretization_degree(akns_discretization);
    leave_fun:
        return degree1step;
}


/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__kdv_discretization_boundary_coeff(kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    REAL bnd_coeff = NAN;
    INT ret_code;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    bnd_coeff = akns_discretization_boundary_coeff(akns_discretization);
    leave_fun:
        return bnd_coeff;
}

/**
 * This routine returns the scaling for effective number of samples based on the discretization.
 */
UINT fnft__kdv_discretization_upsampling_factor(kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    UINT upsampling_factor = 0;
    INT ret_code;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    upsampling_factor = akns_discretization_upsampling_factor(akns_discretization);
    leave_fun:
        return upsampling_factor;
}

/**
 * This routine returns the order of the method based on the discretization.
 */
UINT fnft__kdv_discretization_method_order(kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    UINT order = 0;
    INT ret_code;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    order = akns_discretization_method_order(akns_discretization);
    leave_fun:
        return order;
}

/**
 * This routine returns weights required by some methods based on the discretization.
 */
INT fnft__kdv_discretization_method_weights(COMPLEX **weights_ptr,
        kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    COMPLEX *weights = NULL;
    INT ret_code;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = akns_discretization_method_weights(&weights, akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    *weights_ptr = weights;
    leave_fun:
        return ret_code;
}

/**
 * This akns_discretization related to the given kdv_discretization
 */
INT fnft__kdv_discretization_to_akns_discretization(kdv_discretization_t kdv_discretization,
        akns_discretization_t * const akns_discretization)
{

    switch (kdv_discretization) {
        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT2_MODAL;
            break;
        case kdv_discretization_2SPLIT1A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT1A;
            break;
        case kdv_discretization_2SPLIT1B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT1B;
            break;
        case kdv_discretization_2SPLIT2A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT2A;
            break;
        case kdv_discretization_2SPLIT2B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT2B;
            break;
        case kdv_discretization_2SPLIT2S_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT2S;
            break;
        case kdv_discretization_2SPLIT3A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT3A;
            break;
        case kdv_discretization_2SPLIT3B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT3B;
            break;
        case kdv_discretization_2SPLIT3S_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT3S;
            break;
        case kdv_discretization_2SPLIT4A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT4A;
            break;
        case kdv_discretization_2SPLIT4B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT4B;
            break;
        case kdv_discretization_2SPLIT5A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT5A;
            break;
        case kdv_discretization_2SPLIT5B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT5B;
            break;
        case kdv_discretization_2SPLIT6A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT6A;
            break;
        case kdv_discretization_2SPLIT6B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT6B;
            break;
        case kdv_discretization_2SPLIT7A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT7A;
            break;
        case kdv_discretization_2SPLIT7B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT7B;
            break;
        case kdv_discretization_2SPLIT8A_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT8A;
            break;
        case kdv_discretization_2SPLIT8B_VANILLA:
            *akns_discretization = akns_discretization_2SPLIT8B;
            break;
        case kdv_discretization_BO_VANILLA:
            *akns_discretization = akns_discretization_BO;
            break;
        case kdv_discretization_4SPLIT4A_VANILLA:
            *akns_discretization = akns_discretization_4SPLIT4A;
            break;
        case kdv_discretization_4SPLIT4B_VANILLA:
            *akns_discretization = akns_discretization_4SPLIT4B;
            break;
        case kdv_discretization_CF4_2_VANILLA:
            *akns_discretization = akns_discretization_CF4_2;
            break;
        case kdv_discretization_CF4_3_VANILLA:
            *akns_discretization = akns_discretization_CF4_3;
            break;
        case kdv_discretization_CF5_3_VANILLA:
            *akns_discretization = akns_discretization_CF5_3;
            break;
        case kdv_discretization_CF6_4_VANILLA:
            *akns_discretization = akns_discretization_CF6_4;
            break;
        case kdv_discretization_ES4_VANILLA:
            *akns_discretization = akns_discretization_ES4;
            break;
        case kdv_discretization_TES4_VANILLA:
            *akns_discretization = akns_discretization_TES4;
            break;
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(kdv_discretization);
    }
    return SUCCESS;
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization.
 */
INT fnft__kdv_discretization_lambda_to_z(const UINT n, const REAL eps_t,
        COMPLEX * const vals, kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__akns_discretization_lambda_to_z(n, eps_t, vals, akns_discretization);
    leave_fun:
        return ret_code;
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization.
 */
INT fnft__kdv_discretization_z_to_lambda(const UINT n, const REAL eps_t,
        COMPLEX * const vals, kdv_discretization_t kdv_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__akns_discretization_z_to_lambda(n, eps_t, vals, akns_discretization);
    leave_fun:
        return ret_code;
}

/**
 * This routine returns the phase factor for reflection coefficient (rho).
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__kdv_discretization_phase_factor_rho(const REAL eps_t, const REAL T1,
        REAL * const phase_factor_rho, kdv_discretization_t kdv_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = kdv_discretization_boundary_coeff(kdv_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(kdv_discretization);
    *phase_factor_rho = -2.0*(T1 + eps_t*boundary_coeff);
    return SUCCESS;
}

/**
 * This routine returns the phase factor for a coefficient.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__kdv_discretization_phase_factor_a(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_a, kdv_discretization_t kdv_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = kdv_discretization_boundary_coeff(kdv_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(kdv_discretization);

    switch (kdv_discretization) {

        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT2A_VANILLA:
        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
            *phase_factor_a = -eps_t*D + (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            *phase_factor_a = (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        default: // Unknown discretization

            return E_INVALID_ARGUMENT(kdv_discretization);
    }
}

/**
 * This routine returns the phase factor for b coefficient.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__kdv_discretization_phase_factor_b(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_b, kdv_discretization_t kdv_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = kdv_discretization_boundary_coeff(kdv_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(kdv_discretization);

    UINT degree1step = 0;
    switch (kdv_discretization) {

        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT2A_VANILLA:
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:

            *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
            degree1step = kdv_discretization_degree(kdv_discretization);
            if (degree1step == 0)
                return E_INVALID_ARGUMENT(kdv_discretization);
            *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff) + eps_t/degree1step;
            return SUCCESS;
            break;

        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            *phase_factor_b =  - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        default: // Unknown discretization

            return E_INVALID_ARGUMENT(kdv_discretization);
    }

}


/**
 * This routine preprocess the signal by resampling and subsampling based on the discretization.
 * The preprocessing is necessary for higher-order methods.
 */
INT fnft__kdv_discretization_preprocess_signal(const UINT D, COMPLEX const * const q,
        REAL const eps_t, const INT kappa,
        UINT * const Dsub_ptr, COMPLEX **q_preprocessed_ptr, COMPLEX **r_preprocessed_ptr,
        UINT * const first_last_index,  kdv_discretization_t discretization)
{

    UINT i, D_effective, isub;
    INT ret_code = SUCCESS;
    COMPLEX *q_1 = NULL;
    COMPLEX *q_2 = NULL;
    COMPLEX *q_3 = NULL;

    COMPLEX *r_1 = NULL;
    COMPLEX *r_2 = NULL;
    COMPLEX *r_3 = NULL;
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
    if (abs(kappa) != 1)
        return E_INVALID_ARGUMENT(kappa);
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

    UINT upsampling_factor = kdv_discretization_upsampling_factor(discretization);
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

        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT2A_VANILLA:
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
            i = 0;
            for (isub=0; isub<D_effective; isub++) {
                q_preprocessed[isub] = q[i];
                r_preprocessed[isub] = -1.0;
                i += nskip_per_step;
            }
            break;
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:
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


            ret_code = kdv_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

             i = 0;
            for (isub=0; isub<D_effective; isub=isub+2) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q_2[i];
                q_preprocessed[isub+1] = weights[2]*q_1[i] + weights[3]*q_2[i];
                r_preprocessed[isub] = -1.0 * (weights[0] + weights[1]);
                r_preprocessed[isub+1] = -1.0 * (weights[2] + weights[3]);
                i += nskip_per_step;
            }


            break;
        case kdv_discretization_CF4_3_VANILLA:
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

            ret_code = kdv_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q[i] + weights[2]*q_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i] + weights[4]*q[i] + weights[5]*q_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
                r_preprocessed[isub] = -1.0 * (weights[0] + weights[1] + weights[2]);
                r_preprocessed[isub+1] = -1.0 * (weights[3] + weights[4] + weights[5]);
                r_preprocessed[isub+2] = -1.0 * (weights[6] + weights[7] + weights[8]);
                i += nskip_per_step;
            }
            break;
        case kdv_discretization_CF5_3_VANILLA:
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

            ret_code = kdv_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = weights[0]*q_1[i]  + weights[1]*q[i] + weights[2]*q_3[i];
                r_preprocessed[isub] = -1.0 * (weights[0] + weights[1] + weights[2]);
                q_preprocessed[isub+1] = weights[3]*q_1[i]+ weights[4]*q[i] + weights[5]*q_3[i];
                r_preprocessed[isub+1] = -1.0 * (weights[3] + weights[4] + weights[5]);
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] +  weights[8]*q_3[i];
                r_preprocessed[isub+2] = -1.0 * (weights[6] + weights[7] + weights[8]);
                i += nskip_per_step;
            }


            break;
        case kdv_discretization_CF6_4_VANILLA:
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

            ret_code = kdv_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            i = 0;
            for (isub=0; isub<D_effective; isub=isub+4) {
                q_preprocessed[isub] = weights[0]*q_1[i]  + weights[1]*q[i] + weights[2]*q_3[i];
                r_preprocessed[isub] = -1.0 * (weights[0] + weights[1] + weights[2]);
                q_preprocessed[isub+1] = weights[3]*q_1[i]+ weights[4]*q[i] + weights[5]*q_3[i];
                r_preprocessed[isub+1] = -1.0 * (weights[3] + weights[4] + weights[5]);
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
                r_preprocessed[isub+2] = -1.0 * (weights[6] + weights[7] + weights[8]);
                q_preprocessed[isub+3] = weights[9]*q_1[i] + weights[10]*q[i] +  weights[11]*q_3[i];
                r_preprocessed[isub+3] = -1.0 * (weights[9] + weights[10] + weights[11]);
                i += nskip_per_step;
            }
            break;
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = q[i];
                r_preprocessed[isub] = -1.0;
                i += nskip_per_step;
            }

            REAL eps_t_sub = eps_t*nskip_per_step;
            REAL eps_t_sub_2 = POW(eps_t_sub,2);
            q_preprocessed[1] = (q_preprocessed[3]-0)/(2*eps_t_sub);
            q_preprocessed[2] = (q_preprocessed[3]-2*q_preprocessed[0]+0)/eps_t_sub_2;
            q_preprocessed[D_effective-2] = (0-q_preprocessed[D_effective-6])/(2*eps_t_sub);
            q_preprocessed[D_effective-1] = (0-2*q_preprocessed[D_effective-3]+q_preprocessed[D_effective-6])/eps_t_sub_2;
            r_preprocessed[1] = (r_preprocessed[3]-0)/(2*eps_t_sub);
            r_preprocessed[2] = (r_preprocessed[3]-2*r_preprocessed[0]+0)/eps_t_sub_2;
            r_preprocessed[D_effective-2] = (0-r_preprocessed[D_effective-6])/(2*eps_t_sub);
            r_preprocessed[D_effective-1] = (0-2*r_preprocessed[D_effective-3]+r_preprocessed[D_effective-6])/eps_t_sub_2;

            for (isub=3; isub<D_effective-3; isub=isub+3) {
                q_preprocessed[isub+1] = (q_preprocessed[isub+3]-q_preprocessed[isub-3])/(2*eps_t_sub);
                q_preprocessed[isub+2] = (q_preprocessed[isub+3]-2*q_preprocessed[isub]+q_preprocessed[isub-3])/eps_t_sub_2;
                r_preprocessed[isub+1] = (r_preprocessed[isub+3]-r_preprocessed[isub-3])/(2*eps_t_sub);
                r_preprocessed[isub+2] = (r_preprocessed[isub+3]-2*r_preprocessed[isub]+r_preprocessed[isub-3])/eps_t_sub_2;
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
        free(r_1);
        free(r_2);
        free(r_3);
        free(weights);
        return ret_code;
}

/**
 * This routine returns the change of basis matrix from the basis of the discretization to the S basis.
 */
INT fnft__kdv_change_of_basis_matrix_to_S(COMPLEX * const T,
                                          COMPLEX const xi,
                                          const UINT derivative_flag, // 0- > 2x2, 1->4x4
                                          const REAL eps_t,
                                          kdv_discretization_t kdv_discretization)
{
    INT ret_code = SUCCESS;

    if(xi==0){
        ret_code = E_DIV_BY_ZERO;
        CHECK_RETCODE(ret_code, release_mem);
    }

    switch (kdv_discretization) {
        // Fill the first two columns of T for the discretization basis

        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:
        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            // T from AKNS basis to S basis:
            
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

            break;

        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
        case kdv_discretization_2SPLIT2A_VANILLA:
            //T from modified AKNS basis to S basis. This modification reduces the degree of the polynomial scattering matrix by 1 per step.

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

            break;

//        case
//            // T from companion basis to S basis
//
//            T[0] = 0.5;             //T_11;
//            T[1] = 0.5 * I / xi;    //T_12;
//            if(!derivative_flag){
//                T[2]  = 0.5;        //T_21;
//                T[3]  = -T[1];      //T_22;
//            } else {
//                T[4]  = 0.5;        //T_21;
//                T[5]  = -T[1];      //T_22;
//
//                T[8]  = 0.0;        //T_31 = d/dxi T_11
//                T[9]  = -T[1]/xi;   //T_32 = d/dxi T_12
//                T[12] = 0.0;        //T_41 = d/dxi T_21
//                T[13] = -T[9];      //T_42 = d/dxi T_22
//            }
//
//            break;
//
//        case
//            //
//            T[0] = ; //T_11;
//            T[1] = ; //T_12;
//            if(!derivative_flag){
//                T[2]  = ; //T_21;
//                T[3]  = ; //T_22;
//            } else {
//                T[4]  = ; //T_21;
//                T[5]  = ; //T_22;
//
//                T[8]  = ; //T_31 = d/dxi T_11
//                T[9]  = ; //T_32 = d/dxi T_12
//                T[12] = ; //T_41 = d/dxi T_21
//                T[13] = ; //T_42 = d/dxi T_22
//            }
//
//            break;


        default: // Unknown discretization

            return E_INVALID_ARGUMENT(kdv_discretization);
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
INT fnft__kdv_change_of_basis_matrix_from_S(COMPLEX * const T,
                                            COMPLEX const xi,
                                            UINT const derivative_flag, // 0- > 2x2, 1->4x4
                                            REAL const eps_t,
                                            kdv_discretization_t const kdv_discretization)
{
    INT ret_code = SUCCESS;

    if(xi==0){
        ret_code = E_DIV_BY_ZERO;
        CHECK_RETCODE(ret_code, release_mem);
    }

    switch (kdv_discretization) {
        // Fill the first two columns of T for the discretization basis

        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:
        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            // T from S basis to AKNS basis:

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

            break;

        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
        case kdv_discretization_2SPLIT2A_VANILLA:
            // T from S basis to modified AKNS basis. This modification reduces the degree of the polynomial scattering matrix by 1 per step.

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

            break;

//        case
//            // T from S basis to companion basis
//
//            T[0] = 1.0;            //T_11;
//            T[1] = 1.0;            //T_12;
//            if(!derivative_flag){
//                T[2]  = -I * xi;   //T_21;
//                T[3]  = I * xi;    //T_22;
//            } else {
//                T[4]  = -I * xi;   //T_21;
//                T[5]  = I * xi;    //T_22;
//
//                T[8]  = 0.0;       //T_31 = d/dxi T_11
//                T[9]  = 0.0;       //T_32 = d/dxi T_12
//                T[12] = -I;        //T_41 = d/dxi T_21
//                T[13] = I;         //T_42 = d/dxi T_22
//            }
//
//            break;
//
//        case
//            //
//            T[0] = ; //T_11;
//            T[1] = ; //T_12;
//            if(!derivative_flag){
//                T[2]  = ; //T_21;
//                T[3]  = ; //T_22;
//            } else {
//                T[4]  = ; //T_21;
//                T[5]  = ; //T_22;
//
//                T[8]  = ; //T_31 = d/dxi T_11
//                T[9]  = ; //T_32 = d/dxi T_12
//                T[12] = ; //T_41 = d/dxi T_21
//                T[13] = ; //T_42 = d/dxi T_22
//            }
//
//            break;


        default: // Unknown discretization

            return E_INVALID_ARGUMENT(kdv_discretization);
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
