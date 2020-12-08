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
        case kdv_discretization_2SPLIT2_MODAL:
            *akns_discretization = akns_discretization_2SPLIT2_MODAL;
            break;
        case kdv_discretization_2SPLIT1A_VANILLA:
        case kdv_discretization_2SPLIT1A:
            *akns_discretization = akns_discretization_2SPLIT1A;
            break;
        case kdv_discretization_2SPLIT1B_VANILLA:
        case kdv_discretization_2SPLIT1B:
            *akns_discretization = akns_discretization_2SPLIT1B;
            break;
        case kdv_discretization_2SPLIT2A_VANILLA:
        case kdv_discretization_2SPLIT2A:
            *akns_discretization = akns_discretization_2SPLIT2A;
            break;
        case kdv_discretization_2SPLIT2B_VANILLA:
        case kdv_discretization_2SPLIT2B:
            *akns_discretization = akns_discretization_2SPLIT2B;
            break;
        case kdv_discretization_2SPLIT2S_VANILLA:
        case kdv_discretization_2SPLIT2S:
            *akns_discretization = akns_discretization_2SPLIT2S;
            break;
        case kdv_discretization_2SPLIT3A_VANILLA:
        case kdv_discretization_2SPLIT3A:
            *akns_discretization = akns_discretization_2SPLIT3A;
            break;
        case kdv_discretization_2SPLIT3B_VANILLA:
        case kdv_discretization_2SPLIT3B:
            *akns_discretization = akns_discretization_2SPLIT3B;
            break;
        case kdv_discretization_2SPLIT3S_VANILLA:
        case kdv_discretization_2SPLIT3S:
            *akns_discretization = akns_discretization_2SPLIT3S;
            break;
        case kdv_discretization_2SPLIT4A_VANILLA:
        case kdv_discretization_2SPLIT4A:
            *akns_discretization = akns_discretization_2SPLIT4A;
            break;
        case kdv_discretization_2SPLIT4B_VANILLA:
        case kdv_discretization_2SPLIT4B:
            *akns_discretization = akns_discretization_2SPLIT4B;
            break;
        case kdv_discretization_2SPLIT5A_VANILLA:
        case kdv_discretization_2SPLIT5A:
            *akns_discretization = akns_discretization_2SPLIT5A;
            break;
        case kdv_discretization_2SPLIT5B_VANILLA:
        case kdv_discretization_2SPLIT5B:
            *akns_discretization = akns_discretization_2SPLIT5B;
            break;
        case kdv_discretization_2SPLIT6A_VANILLA:
        case kdv_discretization_2SPLIT6A:
            *akns_discretization = akns_discretization_2SPLIT6A;
            break;
        case kdv_discretization_2SPLIT6B_VANILLA:
        case kdv_discretization_2SPLIT6B:
            *akns_discretization = akns_discretization_2SPLIT6B;
            break;
        case kdv_discretization_2SPLIT7A_VANILLA:
        case kdv_discretization_2SPLIT7A:
            *akns_discretization = akns_discretization_2SPLIT7A;
            break;
        case kdv_discretization_2SPLIT7B_VANILLA:
        case kdv_discretization_2SPLIT7B:
            *akns_discretization = akns_discretization_2SPLIT7B;
            break;
        case kdv_discretization_2SPLIT8A_VANILLA:
        case kdv_discretization_2SPLIT8A:
            *akns_discretization = akns_discretization_2SPLIT8A;
            break;
        case kdv_discretization_2SPLIT8B_VANILLA:
        case kdv_discretization_2SPLIT8B:
            *akns_discretization = akns_discretization_2SPLIT8B;
            break;
        case kdv_discretization_BO_VANILLA:
        case kdv_discretization_BO:
            *akns_discretization = akns_discretization_BO;
            break;
        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4A:
            *akns_discretization = akns_discretization_4SPLIT4A;
            break;
        case kdv_discretization_4SPLIT4B_VANILLA:
        case kdv_discretization_4SPLIT4B:
            *akns_discretization = akns_discretization_4SPLIT4B;
            break;
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_2:
            *akns_discretization = akns_discretization_CF4_2;
            break;
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF4_3:
            *akns_discretization = akns_discretization_CF4_3;
            break;
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF5_3:
            *akns_discretization = akns_discretization_CF5_3;
            break;
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_CF6_4:
            *akns_discretization = akns_discretization_CF6_4;
            break;
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_ES4:
            *akns_discretization = akns_discretization_ES4;
            break;
        case kdv_discretization_TES4_VANILLA:
        case kdv_discretization_TES4:
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
        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
        case kdv_discretization_2SPLIT3S:
        case kdv_discretization_2SPLIT4A:
        case kdv_discretization_2SPLIT4B:
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT6B:
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
        case kdv_discretization_2SPLIT8A:
        case kdv_discretization_2SPLIT8B:
        case kdv_discretization_4SPLIT4A:
        case kdv_discretization_4SPLIT4B:
        case kdv_discretization_2SPLIT2A:
        case kdv_discretization_2SPLIT2_MODAL:
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
        case kdv_discretization_BO: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2:
        case kdv_discretization_CF4_3:
        case kdv_discretization_CF5_3:
        case kdv_discretization_CF6_4:
        case kdv_discretization_ES4:
        case kdv_discretization_TES4:
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

    switch (kdv_discretization) {

        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
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
        case kdv_discretization_2SPLIT2_MODAL:
        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2A:
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
        case kdv_discretization_2SPLIT3S:
        case kdv_discretization_2SPLIT4A:
        case kdv_discretization_2SPLIT4B:
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT6B:
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
        case kdv_discretization_2SPLIT8A:
        case kdv_discretization_2SPLIT8B:
        case kdv_discretization_4SPLIT4A:
        case kdv_discretization_4SPLIT4B:
            *phase_factor_b = -eps_t*D - (T[1]+eps_t*boundary_coeff) - (T[0]-eps_t*boundary_coeff);
            return SUCCESS;
            break;

        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
        case kdv_discretization_BO: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2:
        case kdv_discretization_CF4_3:
        case kdv_discretization_CF5_3:
        case kdv_discretization_CF6_4:
        case kdv_discretization_ES4:
        case kdv_discretization_TES4:
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
COMPLEX fnft__kdv_discretization_r_from_q(COMPLEX const q){(void)q; return -1.0;}
COMPLEX fnft__kdv_discretization_derivatives_r_from_q(COMPLEX const q){(void)q; return 0.0;}
INT fnft__kdv_discretization_preprocess_signal(const UINT D, COMPLEX const * const q,
        REAL const eps_t, const INT kappa,
        UINT * const Dsub_ptr, COMPLEX **q_preprocessed_ptr, COMPLEX **r_preprocessed_ptr,
        UINT * const first_last_index,  kdv_discretization_t discretization)
{
    (void)kappa;
    akns_discretization_t akns_discretization = 0;
    INT ret_code = SUCCESS;
    ret_code = kdv_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    COMPLEX (*r_from_q[3])(COMPLEX) = {
        fnft__kdv_discretization_r_from_q,
        fnft__kdv_discretization_derivatives_r_from_q,
        fnft__kdv_discretization_derivatives_r_from_q
    };

        ret_code = akns_discretization_preprocess_signal(D, q, r_from_q, eps_t, Dsub_ptr, q_preprocessed_ptr, r_preprocessed_ptr, first_last_index, akns_discretization);
        CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;

}

/**
 * This routine returns the change of basis matrix from the basis of the discretization to the S basis.
 */
INT fnft__kdv_discretization_change_of_basis_matrix_to_S(COMPLEX * const T,
                                                         COMPLEX const xi,
                                                         UINT const derivative_flag, // 0- > 2x2, 1->4x4
                                                         REAL const eps_t,
                                                         kdv_discretization_t const kdv_discretization)
{
    INT ret_code = SUCCESS;
    // Note: Currently all discretizations are AKNS-based, so this routine is only a wrapper

    // Fetch the vanilla flag
    UINT vanilla_flag;
    ret_code = kdv_discretization_vanilla_flag(&vanilla_flag, kdv_discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Fetch the AKNS discretizations corresponding to the KdV discretization
    akns_discretization_t akns_discretization;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Call akns_discretization_change_of_basis_matrix_to_S
    ret_code = akns_discretization_change_of_basis_matrix_to_S(T,xi,derivative_flag,eps_t, akns_discretization, vanilla_flag, akns_pde_KdV);
    CHECK_RETCODE(ret_code, release_mem);

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

release_mem:
    return ret_code;
}
/**
 * This routine returns the change of basis matrix from the S basis to the basis of the discretization.
 */
INT fnft__kdv_discretization_change_of_basis_matrix_from_S(COMPLEX * const T,
                                                           COMPLEX const xi,
                                                           UINT const derivative_flag, // 0- > 2x2, 1->4x4
                                                           REAL const eps_t,
                                                           kdv_discretization_t const kdv_discretization)
{
    INT ret_code = SUCCESS;
    // Note: Currently all discretizations are AKNS-based, so this routine is only a wrapper

    // Fetch the vanilla flag
    UINT vanilla_flag;
    ret_code = kdv_discretization_vanilla_flag(&vanilla_flag, kdv_discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Fetch the AKNS discretizations corresponding to the KdV discretization
    akns_discretization_t akns_discretization;
    ret_code = kdv_discretization_to_akns_discretization(kdv_discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Call akns_discretization_change_of_basis_matrix_from_S
    ret_code = akns_discretization_change_of_basis_matrix_from_S(T,xi,derivative_flag,eps_t, akns_discretization, vanilla_flag, akns_pde_KdV);
    CHECK_RETCODE(ret_code, release_mem);

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

release_mem:
    return ret_code;
}

/**
 * This routine tells for discretizations in AKNS basis whether they use r=-1 (vanilla) or q=-1 (not vanilla).
 */
INT fnft__kdv_discretization_vanilla_flag(UINT * const vanilla_flag,
                                          kdv_discretization_t const kdv_discretization)
{
    INT ret_code = SUCCESS;
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
        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
            *vanilla_flag = 1;
            ret_code = SUCCESS;
            CHECK_RETCODE(ret_code, release_mem);
            break;

        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
        case kdv_discretization_2SPLIT3S:
        case kdv_discretization_2SPLIT4A:
        case kdv_discretization_2SPLIT4B:
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT6B:
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
        case kdv_discretization_2SPLIT8A:
        case kdv_discretization_2SPLIT8B:
        case kdv_discretization_4SPLIT4A:
        case kdv_discretization_4SPLIT4B:
        case kdv_discretization_2SPLIT2A:
        case kdv_discretization_2SPLIT2_MODAL:
        case kdv_discretization_BO: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2:
        case kdv_discretization_CF4_3:
        case kdv_discretization_CF5_3:
        case kdv_discretization_CF6_4:
        case kdv_discretization_ES4:
        case kdv_discretization_TES4:
            *vanilla_flag = 0;
            ret_code = SUCCESS;
            CHECK_RETCODE(ret_code, release_mem);
            break;

        default: // Unknown discretization, or non-AKNS discretization
            ret_code = E_INVALID_ARGUMENT(kdv_discretization);
            CHECK_RETCODE(ret_code, release_mem);
    }
    release_mem:
        return ret_code;
}

/**
 * This routine returns the slow discretization to use for the calculation of the discrete spectrum. The returned discritization has the same error order and shares the same preprocessing of the potential.
 */
INT fnft__kdv_slow_discretization(kdv_discretization_t * const kdv_discretization_ptr) {
    INT ret_code = SUCCESS;

    switch (*kdv_discretization_ptr) {
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
        case kdv_discretization_2SPLIT2A_VANILLA:
        case kdv_discretization_2SPLIT2_MODAL_VANILLA:
            // Second order fast vanilla method
            *kdv_discretization_ptr = kdv_discretization_BO_VANILLA;
            break;

        case kdv_discretization_2SPLIT1A:
        case kdv_discretization_2SPLIT1B:
        case kdv_discretization_2SPLIT2B:
        case kdv_discretization_2SPLIT2S:
        case kdv_discretization_2SPLIT3A:
        case kdv_discretization_2SPLIT3B:
        case kdv_discretization_2SPLIT3S:
        case kdv_discretization_2SPLIT4A:
        case kdv_discretization_2SPLIT4B:
        case kdv_discretization_2SPLIT5A:
        case kdv_discretization_2SPLIT5B:
        case kdv_discretization_2SPLIT6A:
        case kdv_discretization_2SPLIT6B:
        case kdv_discretization_2SPLIT7A:
        case kdv_discretization_2SPLIT7B:
        case kdv_discretization_2SPLIT8A:
        case kdv_discretization_2SPLIT8B:
        case kdv_discretization_2SPLIT2A:
        case kdv_discretization_2SPLIT2_MODAL:
            // Second order fast non-vanilla method
            *kdv_discretization_ptr = kdv_discretization_BO;
            break;

        case kdv_discretization_4SPLIT4A_VANILLA:
        case kdv_discretization_4SPLIT4B_VANILLA:
            // Fourth order fast vanilla method
            *kdv_discretization_ptr = kdv_discretization_CF4_2_VANILLA;
            break;

        case kdv_discretization_4SPLIT4A:
        case kdv_discretization_4SPLIT4B:
            // Fourth order fast non-vanilla method
            *kdv_discretization_ptr = kdv_discretization_CF4_2;
            break;

        case kdv_discretization_BO_VANILLA: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2_VANILLA:
        case kdv_discretization_CF4_3_VANILLA:
        case kdv_discretization_CF5_3_VANILLA:
        case kdv_discretization_CF6_4_VANILLA:
        case kdv_discretization_ES4_VANILLA:
        case kdv_discretization_TES4_VANILLA:
        case kdv_discretization_BO: // Bofetta-Osborne scheme
        case kdv_discretization_CF4_2:
        case kdv_discretization_CF4_3:
        case kdv_discretization_CF5_3:
        case kdv_discretization_CF6_4:
        case kdv_discretization_ES4:
        case kdv_discretization_TES4:
            // Allready a slow discretization, leave discretization as is.
            break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(kdv_discretization_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
    }

leave_fun:
    return ret_code;
}
