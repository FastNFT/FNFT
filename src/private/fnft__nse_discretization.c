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
INT fnft__nse_discretization_method_weights(COMPLEX **qr_weights_ptr,
                                            COMPLEX **eps_t_weights_ptr,
                                            nse_discretization_t const nse_discretization)
{
    akns_discretization_t akns_discretization = 0;
    COMPLEX *qr_weights = NULL, *eps_t_weights = NULL;
    INT ret_code;
    ret_code = nse_discretization_to_akns_discretization(nse_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = akns_discretization_method_weights(&qr_weights, &eps_t_weights, akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    *qr_weights_ptr = qr_weights;
    *eps_t_weights_ptr = eps_t_weights;
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
    if (nse_discretization == nse_discretization_2SPLIT2A ||
        nse_discretization == nse_discretization_2SPLIT2_MODAL) {
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
COMPLEX fnft__nse_discretization_focusing_r_from_q(COMPLEX const q){return -CONJ(q);}
COMPLEX fnft__nse_discretization_defocusing_r_from_q(COMPLEX const q){return CONJ(q);}
INT fnft__nse_discretization_preprocess_signal(const UINT D, COMPLEX const * const q,
        REAL const eps_t, const INT kappa,
        UINT * const Dsub_ptr, COMPLEX **q_preprocessed_ptr, COMPLEX **r_preprocessed_ptr,
        UINT * const first_last_index,  nse_discretization_t discretization)
{
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = nse_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    COMPLEX (*r_from_q[3])(COMPLEX);
    for (UINT i = 0; i < 3; i++) {
        if (kappa == 1)
            r_from_q[i] = fnft__nse_discretization_focusing_r_from_q;
        else
            r_from_q[i] = fnft__nse_discretization_defocusing_r_from_q;
    }

    ret_code = akns_discretization_preprocess_signal(D, q, r_from_q, eps_t, Dsub_ptr, q_preprocessed_ptr, r_preprocessed_ptr, first_last_index, akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}
