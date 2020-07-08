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
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_discretization.h"
#include "fnft__akns_discretization.h"

/**
 * Returns the max degree of the polynomials in a single scattering
 * matrix or zero if the discretization is unknown.
 */
UINT fnft__nse_discretization_degree(nse_discretization_t
        nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code;
    UINT degree1step = 0;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    degree1step = akns_discretization_degree(akns_discretization); 
    leave_fun:    
        return degree1step; 
}


/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__nse_discretization_boundary_coeff(nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    REAL bnd_coeff = NAN;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);    
    bnd_coeff = akns_discretization_boundary_coeff(akns_discretization);
    leave_fun:    
        return bnd_coeff;
}

/**
 * This routine returns the scaling for effective number of samples based on the discretization.
 */
UINT fnft__nse_discretization_upsampling_factor(nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    UINT upsampling_factor = 0;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);    
    upsampling_factor = akns_discretization_upsampling_factor(akns_discretization);
    leave_fun:    
        return upsampling_factor;
}

/**
 * This routine returns the order of the method based on the discretization.
 */
UINT fnft__nse_discretization_method_order(nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    UINT order = 0;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);    
    order = akns_discretization_method_order(akns_discretization);
    leave_fun:    
        return order;
}

/**
 * This routine returns weights required by some methods based on the discretization.
 */
INT fnft__nse_discretization_method_weights(COMPLEX **weights_ptr,
        nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    COMPLEX *weights = NULL;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = akns_discretization_method_weights(&weights, akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    *weights_ptr = weights;
    leave_fun:
        return ret_code;
}

/**
 * This akns_discretization related to the given nse_discretization
 */
INT fnft__nse_discretization_to_akns_discretization(nse_discretization_t nse_discretization, 
        akns_discretization_t * const akns_discretization)
{

    switch (nse_discretization) {
        case nse_discretization_2SPLIT2_MODAL:
            *akns_discretization = akns_discretization_2SPLIT2_MODAL;
            break;
        case nse_discretization_2SPLIT1A:
            *akns_discretization = akns_discretization_2SPLIT1A;
            break;
        case nse_discretization_2SPLIT1B:
            *akns_discretization = akns_discretization_2SPLIT1B;
            break;
        case nse_discretization_2SPLIT2A:
            *akns_discretization = akns_discretization_2SPLIT2A;
            break;
        case nse_discretization_2SPLIT2B:
            *akns_discretization = akns_discretization_2SPLIT2B;
            break;
        case nse_discretization_2SPLIT2S:
            *akns_discretization = akns_discretization_2SPLIT2S;
            break;
        case nse_discretization_2SPLIT3S:
            *akns_discretization = akns_discretization_2SPLIT3S;
            break;
        case nse_discretization_2SPLIT4B:
            *akns_discretization = akns_discretization_2SPLIT4B;
            break;
        case nse_discretization_2SPLIT3A:
            *akns_discretization = akns_discretization_2SPLIT3A;
            break;
        case nse_discretization_2SPLIT3B:
            *akns_discretization = akns_discretization_2SPLIT3B;
            break;
        case nse_discretization_2SPLIT4A:
            *akns_discretization = akns_discretization_2SPLIT4A;
            break;
        case nse_discretization_2SPLIT6B:
            *akns_discretization = akns_discretization_2SPLIT6B;
            break;
        case nse_discretization_2SPLIT6A:
            *akns_discretization = akns_discretization_2SPLIT6A;
            break;
        case nse_discretization_2SPLIT8B:
            *akns_discretization = akns_discretization_2SPLIT8B;
            break;
        case nse_discretization_2SPLIT5A:
            *akns_discretization = akns_discretization_2SPLIT5A;
            break;
        case nse_discretization_2SPLIT5B:
            *akns_discretization = akns_discretization_2SPLIT5B;
            break;
        case nse_discretization_2SPLIT8A:
            *akns_discretization = akns_discretization_2SPLIT8A;
            break;
        case nse_discretization_2SPLIT7A:
            *akns_discretization = akns_discretization_2SPLIT7A;
            break;
        case nse_discretization_2SPLIT7B:
            *akns_discretization = akns_discretization_2SPLIT7B;
            break;
        case nse_discretization_BO:
            *akns_discretization = akns_discretization_BO;
            break;   
        case nse_discretization_4SPLIT4A:
            *akns_discretization = akns_discretization_4SPLIT4A;
            break;
        case nse_discretization_4SPLIT4B:
            *akns_discretization = akns_discretization_4SPLIT4B;
            break;        
        case nse_discretization_CF4_2:
            *akns_discretization = akns_discretization_CF4_2;
            break; 
        case nse_discretization_CF4_3:
            *akns_discretization = akns_discretization_CF4_3;
            break; 
        case nse_discretization_CF5_3:
            *akns_discretization = akns_discretization_CF5_3;
            break; 
        case nse_discretization_CF6_4:
            *akns_discretization = akns_discretization_CF6_4;
            break; 
        case nse_discretization_ES4:
            *akns_discretization = akns_discretization_ES4;
            break;
        case nse_discretization_TES4:
            *akns_discretization = akns_discretization_TES4;
            break;   
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(nse_discretization);
    }
    return SUCCESS;    
}
/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization. 
 */
INT fnft__nse_discretization_lambda_to_z(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fnft__akns_discretization_lambda_to_z(n, eps_t, vals, akns_discretization);
    leave_fun:    
        return ret_code;
}

/**
 * This routine maps z from the discrete-time domain to
 * lambda in the continuous-time domain based on the discretization. 
 */
INT fnft__nse_discretization_z_to_lambda(const UINT n, const REAL eps_t, 
        COMPLEX * const vals, nse_discretization_t nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
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
INT fnft__nse_discretization_phase_factor_rho(const REAL eps_t, const REAL T1,
        REAL * const phase_factor_rho, nse_discretization_t nse_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = nse_discretization_boundary_coeff(nse_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(nse_discretization);
    if (nse_discretization == nse_discretization_2SPLIT2A || nse_discretization == nse_discretization_2SPLIT2_MODAL){
        UINT degree1step = nse_discretization_degree(nse_discretization);
        if (degree1step == 0)
            return E_INVALID_ARGUMENT(nse_discretization);
        *phase_factor_rho = -2.0*(T1 + eps_t*boundary_coeff) + eps_t/degree1step;
    }
    else
        *phase_factor_rho = -2.0*(T1 + eps_t*boundary_coeff);
    return SUCCESS;
}

/**
 * This routine returns the phase factor for a coefficient.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__nse_discretization_phase_factor_a(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_a, nse_discretization_t nse_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = nse_discretization_boundary_coeff(nse_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(nse_discretization);
    
    switch (nse_discretization) {
        
        case nse_discretization_2SPLIT1A:
        case nse_discretization_2SPLIT1B:
        case nse_discretization_2SPLIT2B:
        case nse_discretization_2SPLIT2S:
        case nse_discretization_2SPLIT3A:
        case nse_discretization_2SPLIT3B:
        case nse_discretization_2SPLIT3S:
        case nse_discretization_2SPLIT4A:
        case nse_discretization_2SPLIT4B:
        case nse_discretization_2SPLIT5A:
        case nse_discretization_2SPLIT5B:
        case nse_discretization_2SPLIT6A:
        case nse_discretization_2SPLIT6B:
        case nse_discretization_2SPLIT7A:
        case nse_discretization_2SPLIT7B:
        case nse_discretization_2SPLIT8A:
        case nse_discretization_2SPLIT8B:
        case nse_discretization_4SPLIT4A:
        case nse_discretization_4SPLIT4B:
        case nse_discretization_2SPLIT2A:
        case nse_discretization_2SPLIT2_MODAL:
            *phase_factor_a = -eps_t*D + (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;
            
        case nse_discretization_BO: // Bofetta-Osborne scheme
        case nse_discretization_CF4_2:
        case nse_discretization_CF4_3:
        case nse_discretization_CF5_3:
        case nse_discretization_CF6_4:
        case nse_discretization_ES4:
        case nse_discretization_TES4:
            *phase_factor_a = (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;
            
        default: // Unknown discretization
            
            return E_INVALID_ARGUMENT(nse_discretization);
    }
}

/**
 * This routine returns the phase factor for b coefficient.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__nse_discretization_phase_factor_b(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_b, nse_discretization_t nse_discretization)
{
    REAL boundary_coeff;
    boundary_coeff = nse_discretization_boundary_coeff(nse_discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(nse_discretization);

    UINT degree1step = 0;
    switch (nse_discretization) {
        
        case nse_discretization_2SPLIT1A:
        case nse_discretization_2SPLIT1B:
        case nse_discretization_2SPLIT2B:
        case nse_discretization_2SPLIT2S:
        case nse_discretization_2SPLIT3A:
        case nse_discretization_2SPLIT3B:
        case nse_discretization_2SPLIT3S:
        case nse_discretization_2SPLIT4A:
        case nse_discretization_2SPLIT4B:
        case nse_discretization_2SPLIT5A:
        case nse_discretization_2SPLIT5B:
        case nse_discretization_2SPLIT6A:
        case nse_discretization_2SPLIT6B:
        case nse_discretization_2SPLIT7A:
        case nse_discretization_2SPLIT7B:
        case nse_discretization_2SPLIT8A:
        case nse_discretization_2SPLIT8B:
        case nse_discretization_4SPLIT4A:
        case nse_discretization_4SPLIT4B:
            *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;
            
        case nse_discretization_2SPLIT2A:
        case nse_discretization_2SPLIT2_MODAL:
            degree1step = nse_discretization_degree(nse_discretization);
            if (degree1step == 0)
                return E_INVALID_ARGUMENT(nse_discretization);
            *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff) + eps_t/degree1step;
            return SUCCESS;
            break;
            
        case nse_discretization_BO: // Bofetta-Osborne scheme
        case nse_discretization_CF4_2:
        case nse_discretization_CF4_3:
        case nse_discretization_CF5_3:
        case nse_discretization_CF6_4:
        case nse_discretization_ES4:
        case nse_discretization_TES4:
            *phase_factor_b =  - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;
            
        default: // Unknown discretization
            
            return E_INVALID_ARGUMENT(nse_discretization);
    }
    
}


/**
 * This routine preprocess the signal by resampling and subsampling based on the discretization.
 * The preprocessing is necessary for higher-order methods.
 */
INT fnft__nse_discretization_preprocess_signal(const UINT D, COMPLEX const * const q,
        REAL const eps_t, const INT kappa,
        UINT * const Dsub_ptr, COMPLEX **q_preprocessed_ptr, COMPLEX **r_preprocessed_ptr,
        UINT * const first_last_index,  nse_discretization_t discretization)
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
    const UINT nskip_per_step = ROUND((REAL)D / Dsub);
    Dsub = ROUND((REAL)D / nskip_per_step); // actual Dsub
    
    UINT upsampling_factor = nse_discretization_upsampling_factor(discretization);
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
        
        case nse_discretization_BO: // Bofetta-Osborne scheme
        case nse_discretization_2SPLIT1A:
        case nse_discretization_2SPLIT1B:
        case nse_discretization_2SPLIT2A:
        case nse_discretization_2SPLIT2B:
        case nse_discretization_2SPLIT2S:
        case nse_discretization_2SPLIT3A:
        case nse_discretization_2SPLIT3B:
        case nse_discretization_2SPLIT3S:
        case nse_discretization_2SPLIT4A:
        case nse_discretization_2SPLIT4B:
        case nse_discretization_2SPLIT5A:
        case nse_discretization_2SPLIT5B:
        case nse_discretization_2SPLIT6A:
        case nse_discretization_2SPLIT6B:
        case nse_discretization_2SPLIT7A:
        case nse_discretization_2SPLIT7B:
        case nse_discretization_2SPLIT8A:
        case nse_discretization_2SPLIT8B:
        case nse_discretization_2SPLIT2_MODAL:    
            i = 0;
            for (isub=0; isub<D_effective; isub++) {
                q_preprocessed[isub] = q[i];
                r_preprocessed[isub] = -kappa*CONJ(q[i]);
                i += nskip_per_step;
            }
            break;
        case nse_discretization_CF4_2:
        case nse_discretization_4SPLIT4A:
        case nse_discretization_4SPLIT4B:    
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
            
            
            ret_code = nse_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);
            
             i = 0;
            for (isub=0; isub<D_effective; isub=isub+2) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q_2[i];
                q_preprocessed[isub+1] = weights[2]*q_1[i] + weights[3]*q_2[i];
                r_preprocessed[isub] = -kappa*CONJ(q_preprocessed[isub]);
                r_preprocessed[isub+1] = -kappa*CONJ(q_preprocessed[isub+1]);
                i += nskip_per_step;
            }
             
             
            break;
        case nse_discretization_CF4_3:
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

            ret_code = nse_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);

            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = weights[0]*q_1[i] + weights[1]*q[i] + weights[2]*q_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i] + weights[4]*q[i] + weights[5]*q_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
                r_preprocessed[isub] = -kappa*CONJ(q_preprocessed[isub]);
                r_preprocessed[isub+1] = -kappa*CONJ(q_preprocessed[isub+1]);
                r_preprocessed[isub+2] = -kappa*CONJ(q_preprocessed[isub+2]);
                i += nskip_per_step;
            }
            break;
        case nse_discretization_CF5_3:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_3 = malloc(D * sizeof(COMPLEX));
            r_1 = malloc(D * sizeof(COMPLEX));
            r_2 = malloc(D * sizeof(COMPLEX));
            r_3 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_3 == NULL || r_1 == NULL || r_2 == NULL || r_3 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }
            
            ret_code = misc_resample(D, eps_t, q, -eps_t*SQRT(15.0)/10.0*nskip_per_step, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*SQRT(15.0)/10.0*nskip_per_step, q_3);
            CHECK_RETCODE(ret_code, release_mem);
            
            for (i=0; i < D; i++){
                r_1[i] = -kappa*CONJ(q_1[i]);
                r_2[i] = -kappa*CONJ(q[i]);
                r_3[i] = -kappa*CONJ(q_3[i]);
            }
            
            
            ret_code = nse_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);
          
            
            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = weights[0]*q_1[i]  + weights[1]*q[i] + weights[2]*q_3[i];
                r_preprocessed[isub] = weights[0]*r_1[i]  + weights[1]*r_2[i] + weights[2]*r_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i]+ weights[4]*q[i] + weights[5]*q_3[i];
                r_preprocessed[isub+1] = weights[3]*r_1[i]+ weights[4]*r_2[i] + weights[5]*r_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] +  weights[8]*q_3[i];
                r_preprocessed[isub+2] = weights[6]*r_1[i] + weights[7]*r_2[i] +  weights[8]*r_3[i];
                i += nskip_per_step;
            }
            
            
            break;
        case nse_discretization_CF6_4:
            q_1 = malloc(D * sizeof(COMPLEX));
            q_3 = malloc(D * sizeof(COMPLEX));
            r_1 = malloc(D * sizeof(COMPLEX));
            r_2 = malloc(D * sizeof(COMPLEX));
            r_3 = malloc(D * sizeof(COMPLEX));
            if (q_1 == NULL || q_3 == NULL || r_1 == NULL || r_2 == NULL || r_3 == NULL) {
                ret_code = E_NOMEM;
                goto release_mem;
            }
            
            ret_code = misc_resample(D, eps_t, q, -eps_t*nskip_per_step*SQRT(15.0)/10.0, q_1);
            CHECK_RETCODE(ret_code, release_mem);
            ret_code = misc_resample(D, eps_t, q, eps_t*nskip_per_step*SQRT(15.0)/10.0, q_3);
            CHECK_RETCODE(ret_code, release_mem);
            
            for (i=0; i < D; i++){
                r_1[i] = -kappa*CONJ(q_1[i]);
                r_2[i] = -kappa*CONJ(q[i]);
                r_3[i] = -kappa*CONJ(q_3[i]);
            }
            
            ret_code = nse_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, release_mem);
            
            i = 0;
            for (isub=0; isub<D_effective; isub=isub+4) {
                q_preprocessed[isub] = weights[0]*q_1[i]  + weights[1]*q[i] + weights[2]*q_3[i];
                r_preprocessed[isub] = weights[0]*r_1[i]  + weights[1]*r_2[i] + weights[2]*r_3[i];
                q_preprocessed[isub+1] = weights[3]*q_1[i]+ weights[4]*q[i] + weights[5]*q_3[i];
                r_preprocessed[isub+1] = weights[3]*r_1[i]+ weights[4]*r_2[i] + weights[5]*r_3[i];
                q_preprocessed[isub+2] = weights[6]*q_1[i] + weights[7]*q[i] + weights[8]*q_3[i];
                r_preprocessed[isub+2] = weights[6]*r_1[i] + weights[7]*r_2[i] +  weights[8]*r_3[i];
                q_preprocessed[isub+3] = weights[9]*q_1[i] + weights[10]*q[i] +  weights[11]*q_3[i];
                r_preprocessed[isub+3] = weights[9]*r_1[i] + weights[10]*r_2[i] +  weights[11]*r_3[i];
                i += nskip_per_step;
            }
            break;
        case nse_discretization_ES4:
        case nse_discretization_TES4:
            i = 0;
            for (isub=0; isub<D_effective; isub=isub+3) {
                q_preprocessed[isub] = q[i];
                i += nskip_per_step;
            }
            
            REAL eps_t_sub = eps_t*nskip_per_step;
            REAL eps_t_sub_2 = POW(eps_t_sub,2);        
            q_preprocessed[1] = (q_preprocessed[3]-0)/(2*eps_t_sub); 
            q_preprocessed[2] = (q_preprocessed[3]-2*q_preprocessed[0]+0)/eps_t_sub_2;
            q_preprocessed[D_effective-2] = (0-q_preprocessed[D_effective-6])/(2*eps_t_sub); 
            q_preprocessed[D_effective-1] = (0-2*q_preprocessed[D_effective-3]+q_preprocessed[D_effective-6])/eps_t_sub_2;
            

            for (isub=3; isub<D_effective-3; isub=isub+3) {
                q_preprocessed[isub+1] = (q_preprocessed[isub+3]-q_preprocessed[isub-3])/(2*eps_t_sub); 
                q_preprocessed[isub+2] = (q_preprocessed[isub+3]-2*q_preprocessed[isub]+q_preprocessed[isub-3])/eps_t_sub_2;
            }

            for (i=0; i<D_effective; i++) {
                r_preprocessed[i] = -kappa*CONJ(q_preprocessed[i]);
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