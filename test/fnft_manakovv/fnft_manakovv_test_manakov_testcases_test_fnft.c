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
* Contributor:
* Lianne de Vries (TU Delft student) 2021
* Sander Wahls (TU Delft) 2021
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__manakovv_testcases.h"
#include "fnft_manakovv.h"

INT main(){

    manakovv_testcases_t tc = manakovv_testcases_SECH_FOCUSING;
    UINT D = 16;
    REAL error_bounds[7] = {200, 200, 200, 200, 200, 200, 200};
    fnft_manakovv_opts_t opts = fnft_manakovv_default_opts();
    INT ret_code = manakovv_testcases_test_fnft(tc, D, error_bounds, &opts);
    if (ret_code == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

