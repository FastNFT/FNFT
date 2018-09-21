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
 * Shrinivas Chimmalgi (TU Delft) 2018.
 */

#include "fnft_nsev_inverse_test_truncated_soliton.inc"


int main()
{
    UINT M;
    UINT D;
    REAL error_bound;
    INT ret_code = SUCCESS;

    fnft_nsev_inverse_opts_t opts = fnft_nsev_inverse_default_opts();
    opts.discretization = nse_discretization_2SPLIT2_MODAL;
    opts.discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANTS;


    D = 512;
    M = 4*D;
    error_bound = 0.0033;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D = D*2;
    M =4*D;
    error_bound = error_bound/2;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts);
    CHECK_RETCODE(ret_code, leave_fun);



leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
