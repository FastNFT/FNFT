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
 * Sander Wahls (TU Delft) 2018, 2020.
 */

#include "fnft__kdv_finvscatter_test.inc"

int main(void) {
    const kdv_discretization_t discretization = kdv_discretization_2SPLIT2A;

    UINT D = 8;
    REAL error_bound = 1.8e12*FNFT_EPSILON;
    if (kdv_finvscatter_test(D, error_bound, discretization) != SUCCESS)
        return EXIT_FAILURE;

    D = 16384;
    error_bound = 4e9*FNFT_EPSILON;
    if (kdv_finvscatter_test(D, error_bound, discretization) != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
