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
 * Sander Wahls (TU Delft) 2018.
 */

#include "fnft_nsev_inverse_test_sech_defocusing.inc"

int main()
{
    REAL error_bound;
    INT ret_code;

    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();
    opts.discretization = nse_discretization_2SPLIT2A;

    opts.contspec_inversion_method
        = fnft_nsev_inverse_contspec_inversion_method_REFL_COEFF;

    error_bound = 0.0015;
    ret_code = fnft_nsev_inverse_test(error_bound, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.contspec_inversion_method
        = fnft_nsev_inverse_contspec_inversion_method_B_FROM_A;

    error_bound = 0.0015;
    ret_code = fnft_nsev_inverse_test(error_bound, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.contspec_inversion_method
        = fnft_nsev_inverse_contspec_inversion_method_B_FROM_A_WO_SPECFACT;

    error_bound = 0.0015;
    ret_code = fnft_nsev_inverse_test(error_bound, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    return SUCCESS;
}
