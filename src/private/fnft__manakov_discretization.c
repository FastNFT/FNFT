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
* Contributor:
* Lianne de Vries (TU Delft) 2021.
*/


#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakov_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 * 
 */
UINT fnft__manakov_discretization_degree(manakov_discretization_t
        discretization)
{
    switch (discretization) {
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
            return 6;
        case manakov_discretization_2SPLIT4B:
		case manakov_discretization_4SPLIT4B:
			return 4;
        case manakov_discretization_2SPLIT4A:            
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
            return 8;
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_4SPLIT6B:
            return 12;

                    
        default: // Unknown discretization
            return 0;
    }
}

/**
 * This routine returns the scaling for effective number of samples based on the discretization.
 */
UINT fnft__manakov_discretization_upsampling_factor(manakov_discretization_t discretization)
{
    switch (discretization) {
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
            return 2;
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT4A:            
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
        case manakov_discretization_2SPLIT6B:
            return 1;

        default: // discretization not implemented for manakov
            return 0;
    }
}

UINT fnft__manakov_discretization_method_order(manakov_discretization_t discretization){
    switch (discretization) {
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT4A:
        case manakov_discretization_2SPLIT6B:
            return 2;
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
            return 4;

        default: // discretization not implemented for manakov
            return 0;
    }
}

// This function does what nse_discretiation_preprocess_signal and akns_discretiation_preprocess_signal
// do together in the NSE case: downsample if desired (needed for Richardson extrapolation) and 
// get the required samples after interpolation if the chosen method need those.
UINT manakov_discretization_preprocess_signal(const UINT D, COMPLEX const * const q1,
        COMPLEX const * const q2, REAL const eps_t, const INT kappa,
        UINT * const Dsub_ptr, COMPLEX **q1_preprocessed_ptr, COMPLEX **q2_preprocessed_ptr, 
        UINT * const first_last_index,  manakov_discretization_t discretization){

    UINT isub, i, D_effective;

            // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q1 == NULL)
        return E_INVALID_ARGUMENT(q1);
    if (q2 == NULL)
        return E_INVALID_ARGUMENT(q2);
    if (Dsub_ptr == NULL)
        return E_INVALID_ARGUMENT(Dsub_ptr);
    if (q1_preprocessed_ptr == NULL)
        return E_INVALID_ARGUMENT(q1_preprocessed_ptr);
    if (q2_preprocessed_ptr == NULL)
        return E_INVALID_ARGUMENT(q2_preprocessed_ptr);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(esp_t);
    if (first_last_index == NULL)
        return E_INVALID_ARGUMENT(first_last_index);


    // Determine number of samples after downsampling, Dsub
    // If Dsub = D, we o not downsample
    UINT Dsub = *Dsub_ptr; // desired Dsub
    if (Dsub < 2)
        Dsub = 2;
    if (Dsub > D)
        Dsub = D;
    const UINT nskip_per_step = (UINT) ROUND((REAL)D / Dsub);
    Dsub = (UINT) ROUND((REAL)D / nskip_per_step); // actual Dsub

    UINT upsampling_factor = manakov_discretization_upsampling_factor(discretization);
    INT ret_code = SUCCESS;
    if (upsampling_factor == 0){
        ret_code =  E_INVALID_ARGUMENT(discretization);
        goto release_mem;
    }
    D_effective = Dsub * upsampling_factor;
    COMPLEX * const q1_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    COMPLEX * const q2_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    COMPLEX * const r1_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    COMPLEX * const r2_preprocessed = malloc(D_effective * sizeof(COMPLEX));
    COMPLEX *q1_c1, *q1_c2, *q2_c1, *q2_c2;
            q1_c1 = malloc(D*sizeof(COMPLEX));
            q2_c1 = malloc(D*sizeof(COMPLEX));
            q1_c2 = malloc(D*sizeof(COMPLEX));
            q2_c2 = malloc(D*sizeof(COMPLEX));
    if (q1_preprocessed == NULL || q2_preprocessed == NULL || 
            r1_preprocessed == NULL || r2_preprocessed == NULL ) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // actually getting the samples
    switch (discretization){
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT4A:            
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
        case manakov_discretization_2SPLIT6B:
            for (isub=0, i=0; isub<D_effective; isub++, i += nskip_per_step) {  // downsampling
                q1_preprocessed[isub] = q1[i];
                q2_preprocessed[isub] = q2[i];
                r1_preprocessed[isub] = conj(q1[i])*-1*kappa;
                r2_preprocessed[isub] = conj(q2[i])*-1*kappa;
            }
            break;
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
    {
        // Getting bandlimited interpolated samples. First get the samples and 
        // then subsampling
            const REAL a1 = 0.25 + sqrt(3)/6;
            const REAL a2 = 0.25 - sqrt(3)/6;
            
            // Getting non-equidistant samples. qci has sample points at T0+(n+ci)*eps_t,
            // c1 = 0.5-sqrt(3)/6, c2 = 0.5+sqrt(3)/6
            // Getting new samples:
            misc_resample(D, eps_t, q1, (0.5-SQRT(3)/6)*eps_t*nskip_per_step, q1_c1);
            misc_resample(D, eps_t, q1, (0.5+SQRT(3)/6)*eps_t*nskip_per_step, q1_c2);
            misc_resample(D, eps_t, q2, (0.5-SQRT(3)/6)*eps_t*nskip_per_step, q2_c1);
            misc_resample(D, eps_t, q2, (0.5+SQRT(3)/6)*eps_t*nskip_per_step, q2_c2);

            // now select the right values:
            for (isub=0, i=0; isub<D_effective; isub+=2, i+= nskip_per_step) {
                q1_preprocessed[isub] = a1*q1_c1[i]+a2*q1_c2[i];
                q2_preprocessed[isub] = a1*q2_c1[i]+a2*q2_c2[i];
                q1_preprocessed[isub+1] = a2*q1_c1[i]+a1*q1_c2[i];
                q2_preprocessed[isub+1] = a2*q2_c1[i]+a1*q2_c2[i];
            }
            for (i=0; i<D_effective;){
                r1_preprocessed[i] =conj(q1_preprocessed[i])*-1*kappa;
                r2_preprocessed[i] =conj(q2_preprocessed[i])*-1*kappa;
            }
    }
            break;

            default: // Unknown discretization for manakov

            ret_code = E_INVALID_ARGUMENT(discretization);
            goto  release_mem;        
    }
    first_last_index[0] = 0;
    first_last_index[1] = (Dsub-1)*nskip_per_step;
    *q1_preprocessed_ptr = q1_preprocessed;
    *q2_preprocessed_ptr = q2_preprocessed;
    *Dsub_ptr = Dsub;

    release_mem:        // TODO: free memory, do I need to take care of more variables?
        free(q1_c1);
        free(q1_c2);
        free(q2_c1);
        free(q2_c2);
        return ret_code;
}


/**
 * This routine maps lambda from continuous-time domain to
 * z in the discrete-time domain based on the discretization.
 */
INT fnft__manakov_discretization_lambda_to_z(const UINT n, const REAL eps_t,
        COMPLEX * const vals, manakov_discretization_t discretization)
{
    REAL degree1step;
    UINT i, upsampling_factor;
    upsampling_factor = manakov_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = manakov_discretization_degree(discretization);
    if (degree1step == 0)
        return E_INVALID_ARGUMENT(discretization);
    degree1step = degree1step * upsampling_factor;
    for (i = 0; i < n; i++)
        vals[i] = CEXP(2*I*vals[i]*eps_t/degree1step);
    return SUCCESS;
}


/**
 * This routine returns the phase factor for reflection coefficient (rho).
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__manakov_discretization_phase_factor_rho(const REAL eps_t, const REAL T1,
        REAL * const phase_factor_rho, manakov_discretization_t discretization)
{
    REAL boundary_coeff;
    boundary_coeff = manakov_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(discretization);

    switch (discretization) {

        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4A:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
            *phase_factor_rho = -2.0*(T1 + eps_t*boundary_coeff);
            return SUCCESS;
        break;

        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
            *phase_factor_rho = -2.0*(T1 + eps_t*boundary_coeff)-eps_t;
            return SUCCESS;
        break;

        default: // Unknown discretization

            return E_INVALID_ARGUMENT(nse_discretization);
    }

}

/**
 * This routine returns the phase factor for a coefficient.
 * It is required for applying boundary conditions to the transfer matrix values
 * based on the discretization.
 */
INT fnft__manakov_discretization_phase_factor_a(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_a, manakov_discretization_t discretization)
{
    REAL boundary_coeff;
    boundary_coeff = manakov_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(discretization);

    switch (discretization) {

        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4A:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
//            *phase_factor_a =2*(T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);   // same as in matlab, does not give better results
            *phase_factor_a = -eps_t*D + (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
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
INT fnft__manakov_discretization_phase_factor_b(const REAL eps_t, const UINT D, REAL const * const T,
        REAL * const phase_factor_b, manakov_discretization_t discretization)
{
    REAL boundary_coeff;
    boundary_coeff = manakov_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(discretization);
    UINT degree1step = 0;
    switch (discretization) {

        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4A:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
                    degree1step = manakov_discretization_degree(discretization);
            if (degree1step == 0)
                return E_INVALID_ARGUMENT(discretization);
//            *phase_factor_b =  - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
           *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
//           *phase_factor_b = -2*(T[0]-eps_t*boundary_coeff);           // matlab expression
            return SUCCESS;
            break;
        
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
            *phase_factor_b = -eps_t*(D+1) - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        default: // Unknown discretization

            return E_INVALID_ARGUMENT(nse_discretization);
    }
}

/**
 * This routine returns the boundary coefficient based on the discretization.
 */
REAL fnft__manakov_discretization_boundary_coeff(manakov_discretization_t discretization)
{
    
    switch (discretization) {
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4A:
        case manakov_discretization_2SPLIT4B:
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
        case manakov_discretization_4SPLIT6B:
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_FTES4_4A:
        case manakov_discretization_FTES4_4B:
            return 0.5;
            
        default: // Unknown discretization
            return NAN;
    }
}





