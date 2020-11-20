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
