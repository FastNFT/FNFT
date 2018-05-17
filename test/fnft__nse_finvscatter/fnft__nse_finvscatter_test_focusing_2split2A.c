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

#include "fnft__nse_finvscatter_test.inc"

int main()
{
    const INT kappa = +1;
    const REAL error_bound
#ifdef HAVE_FFTW3
        = 61.0*FNFT_EPSILON;
#else
        = 1262.0*FNFT_EPSILON;
#endif    
    const nse_discretization_t discretization
        = fnft_nse_discretization_2SPLIT2A;

    if (nse_finvscatter_test(kappa, error_bound, discretization) != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

