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

#include "fnft_nsev_inverse_test_against_forward.inc"

int main()
{
    const INT kappa = -1;
    UINT M;
    REAL error_bound;
    INT ret_code = SUCCESS;

    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();
    opts.discretization = nse_discretization_2SPLIT2_MODAL;
    opts.contspec_inversion_method
        = fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_REFL_COEFF;

    M = 10;
    error_bound = 0.017;
    ret_code = fnft_nsev_inverse_test(M, error_bound, kappa, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    M = 32;
    error_bound = 7.5e-9;
    ret_code = fnft_nsev_inverse_test(M, error_bound, kappa, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    opts.contspec_inversion_method
        = fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER;

    M = 8;
    error_bound = 1.1e-14;
    ret_code = fnft_nsev_inverse_test(M, error_bound, kappa, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
