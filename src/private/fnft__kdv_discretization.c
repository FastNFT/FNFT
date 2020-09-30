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
* Shrinivas Chimmalgi (TU Delft) 2017-2019.
* Peter J. Prins (TU Delft) 2018.
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
 * This akns_discretization related to the given kdv_discretization
 */
INT fnft__kdv_discretization_to_akns_discretization(kdv_discretization_t kdv_discretization,
        akns_discretization_t * const akns_discretization)
{
    switch (kdv_discretization) {
        case kdv_discretization_2SPLIT2_MODAL:
             *akns_discretization = akns_discretization_2SPLIT2_MODAL;
             break;
        case kdv_discretization_2SPLIT2_MODAL2:
             *akns_discretization = akns_discretization_2SPLIT2_MODAL2;
             break;
        case kdv_discretization_2SPLIT1A:
            *akns_discretization = akns_discretization_2SPLIT1A;
            break;
        case kdv_discretization_2SPLIT1B:
            *akns_discretization = akns_discretization_2SPLIT1B;
            break;
        case kdv_discretization_2SPLIT2A:
            *akns_discretization = akns_discretization_2SPLIT2A;
            break;
        case kdv_discretization_2SPLIT2B:
            *akns_discretization = akns_discretization_2SPLIT2B;
            break;
        case kdv_discretization_2SPLIT2S:
            *akns_discretization = akns_discretization_2SPLIT2S;
            break;
        case kdv_discretization_2SPLIT3S:
            *akns_discretization = akns_discretization_2SPLIT3S;
            break;
        case kdv_discretization_2SPLIT4B:
            *akns_discretization = akns_discretization_2SPLIT4B;
            break;
        case kdv_discretization_2SPLIT3A:
            *akns_discretization = akns_discretization_2SPLIT3A;
            break;
        case kdv_discretization_2SPLIT3B:
            *akns_discretization = akns_discretization_2SPLIT3B;
            break;
        case kdv_discretization_2SPLIT4A:
            *akns_discretization = akns_discretization_2SPLIT4A;
            break;
        case kdv_discretization_2SPLIT6B:
            *akns_discretization = akns_discretization_2SPLIT6B;
            break;
        case kdv_discretization_2SPLIT6A:
            *akns_discretization = akns_discretization_2SPLIT6A;
            break;
        case kdv_discretization_2SPLIT8B:
            *akns_discretization = akns_discretization_2SPLIT8B;
            break;
        case kdv_discretization_2SPLIT5A:
            *akns_discretization = akns_discretization_2SPLIT5A;
            break;
        case kdv_discretization_2SPLIT5B:
            *akns_discretization = akns_discretization_2SPLIT5B;
            break;
        case kdv_discretization_2SPLIT8A:
            *akns_discretization = akns_discretization_2SPLIT8A;
            break;
        case kdv_discretization_2SPLIT7A:
            *akns_discretization = akns_discretization_2SPLIT7A;
            break;
        case kdv_discretization_2SPLIT7B:
            *akns_discretization = akns_discretization_2SPLIT7B;
            break;
        case kdv_discretization_BO:
            *akns_discretization = akns_discretization_BO;
            break;
        case kdv_discretization_4SPLIT4A:
            *akns_discretization = akns_discretization_4SPLIT4A;
            break;
        case kdv_discretization_4SPLIT4B:
            *akns_discretization = akns_discretization_4SPLIT4B;
            break;
        case kdv_discretization_CF4_2:
            *akns_discretization = akns_discretization_CF4_2;
            break;
        case kdv_discretization_CF4_3:
            *akns_discretization = akns_discretization_CF4_3;
            break;
        case kdv_discretization_CF5_3:
            *akns_discretization = akns_discretization_CF5_3;
            break;
        case kdv_discretization_CF6_4:
            *akns_discretization = akns_discretization_CF6_4;
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
