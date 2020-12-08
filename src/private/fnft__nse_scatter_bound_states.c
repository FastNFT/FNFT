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
 * Sander Wahls (TU Delft) 2017-2018, 2020.
 * Shrinivas Chimmalgi (TU Delft) 2017, 2019-2020.
 * Marius Brehler (TU Dortmund) 2018.
 * Peter J Prins (TU Delft) 2020.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 */
INT nse_scatter_bound_states(UINT const D,
                             COMPLEX const * const q,
                             COMPLEX const * const r,
                             REAL const * const T,
                             UINT const K,
                             COMPLEX * const bound_states,
                             COMPLEX * const a_vals,
                             COMPLEX * const aprime_vals,
                             COMPLEX * const b,
                             nse_discretization_t const discretization,
                             UINT const skip_b_flag)
{
    INT ret_code = SUCCESS;

    // Fetch the AKNS discretizations corresponding to the KdV discretization
    akns_discretization_t akns_discretization;
    ret_code = nse_discretization_to_akns_discretization(discretization, &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    // Call akns_scatter_bound_states
    UINT vanilla_flag = 0; // Ignored value for NSE
    ret_code = akns_scatter_bound_states(D, q, r, T, K,bound_states, a_vals, aprime_vals, b, akns_discretization, akns_pde_NSE, vanilla_flag, skip_b_flag);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return ret_code;
}
