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
 * Shrinivas Chimmalgi (TU Delft) 2018.
 *
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__nse_fscatter.h"
#include "fnft__misc.h"

/**
 * Returns the length of array to be allocated based on the number
 * of samples and discretization.
 */

UINT nse_fscatter_numel(UINT D, nse_discretization_t discretization)
{

    const UINT deg = nse_discretization_degree(discretization);
    if (deg == 0)
        return 0; // unknown discretization
    else
        return poly_fmult2x2_numel(deg, D);
}

INT nse_fscatter(const UINT D, COMPLEX const * const q, COMPLEX const * const r,
        const REAL eps_t, COMPLEX * const result, UINT * const deg_ptr,
        INT * const W_ptr, nse_discretization_t discretization)
{
    INT ret_code = SUCCESS;
    akns_discretization_t akns_discretization;

    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(eps_t);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (deg_ptr == NULL)
        return E_INVALID_ARGUMENT(deg_ptr);

    ret_code = nse_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = akns_fscatter(D, q, r, eps_t, result, deg_ptr, W_ptr, akns_discretization);

leave_fun:
    return ret_code;
}
