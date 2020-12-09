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

#include "fnft__nse_finvscatter_test.inc"

int main()
{
    const INT kappa = -1;
    const nse_discretization_t discretization
        = fnft_nse_discretization_2SPLIT2_MODAL;

    UINT D = 8;
    REAL error_bound = 5.0*FNFT_EPSILON;
    if (nse_finvscatter_test(D, kappa, error_bound, discretization) != SUCCESS)
        return EXIT_FAILURE;

    D = 16384;
    error_bound
#ifdef HAVE_FFTW3
        = 483.0*FNFT_EPSILON;
#else
        = 1253.0*FNFT_EPSILON;
#endif
    if (nse_finvscatter_test(D, kappa, error_bound, discretization) != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
