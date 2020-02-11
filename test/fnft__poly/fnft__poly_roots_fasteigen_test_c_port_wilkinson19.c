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
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__poly_roots_fasteigen.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#ifdef DEBUG
#include <stdio.h>
#endif

// Test case is a Wilkinson polynomial with roots 1,2,3,...,19

INT poly_roots_fasteigen_test_c_port_wilkinson19()
{
    const UINT deg = 19;
    COMPLEX p[20] = {
        1,
        -190,
        16815,
        -920550,
        34916946,
        -973941900,
        20692933630,
        -342252511900,
        4465226757381,
        -46280647751910,
        381922055502195,
        -2503858755467550,
        12953636989943896,
        -52260903362512720,
        161429736530118960,
        -371384787345228000,
        610116075740491776,
        -668609730341153280,
        431565146817638400,
        -121645100408832000
    };
    COMPLEX roots[19] = {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
    };
    COMPLEX roots_fortran[19];
    COMPLEX roots_c_port[19];
    INT ret_code;

    ret_code = poly_roots_fasteigen_(deg, p, roots_fortran);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);

    ret_code = poly_roots_fasteigen(deg, p, roots_c_port);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);

    // check that C port gives same result as Fortran version
    for (UINT i=0; i<deg; i++) {
        const REAL diff = CABS(roots_fortran[i] - roots_c_port[i]);
        const REAL tol = 10000*EPSILON*CABS(roots_fortran[i]);
        if (!(diff <= tol)) {
#ifdef DEBUG
            printf("diff = %g > tol = %g\n", diff, tol);
#endif
            return E_TEST_FAILED;
        }
    }

    // check that roots are close to results from Octave
    const REAL dist = misc_hausdorff_dist(deg, roots_c_port, deg, roots);
    if (!(dist <= 0.0155)) { // seems ok, Matlab roots results in 0.0101
#ifdef DEBUG
        printf("dist = %g\n", dist);
#endif
        return E_TEST_FAILED;
    }

    return SUCCESS;
}

INT main()
{
    if ( poly_roots_fasteigen_test_c_port_wilkinson19() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
