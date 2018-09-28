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

#include "fnft_nsev_inverse_test_against_forward_w_discrete.inc"


int main()
{
    UINT M;
    UINT D;
    REAL error_bound;
    INT ret_code = SUCCESS;

    fnft_nsev_inverse_opts_t opts_inv = fnft_nsev_inverse_default_opts();
    opts_inv.discretization = nse_discretization_2SPLIT2A;
    opts_inv.discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANTS;      

    D = 512;
    M = 2*D;
    error_bound = 0.014;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D = D*2;
    M = 2*D;
    error_bound = error_bound/4;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);

    opts_inv.discspec_type = fnft_nsev_inverse_dstype_RESIDUES;      

    D = 512;
    M = 2*D;
    error_bound = 0.014;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);
    
    D = D*2;
    M = 2*D;
    error_bound = error_bound/4;
    ret_code = fnft_nsev_inverse_test(D, M, error_bound, &opts_inv);
    CHECK_RETCODE(ret_code, leave_fun);



leave_fun:
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
