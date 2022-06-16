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
* Shrinivas Chimmalgi (TU Delft) 2020.
* Sander Wahls (TU Delft) 2022.
*/

/**
 * @file fnft__global_root_finding_algorithm.h
 * @brief Global root finding algorithm.
 * @ingroup roots
 */

#ifndef FNFT__GLOBAL_ROOT_POLE_FINDING_ALGORITHM
#define FNFT__GLOBAL_ROOT_POLE_FINDING_ALGORITHM

#include "fnft.h"

/**
 * @brief Global root pole finding algorithm.
 * @ingroup roots
 *
 * Finds the roots and poles of a function within a finite rectangular
 * region of the complex plane.
 *
 * This code is based on the Matlab library GRPF by Piotr Kowalczyk from
 * Gdansk University of Technology that is available under the MIT license
 * at https://github.com/PioKow/GRPF. For more information, please visit
 * the website or see the papers
 *
 * - P. Kowalczyk, “Complex Root Finding Algorithm Based on Delaunay
 *   Triangulation”, ACM Transactions on Mathematical Software, vol. 41,
 *   no. 3, art. 19, pp. 1-13, June 2015
 * - P. Kowalczyk, "Global Complex Roots and Poles Finding Algorithm Based
 *   on Phase Analysis for Propagation and Radiation Problems," IEEE
 *   Transactions on Antennas and Propagation, vol. 66, no. 12, pp. 7198-7205,
 *   Dec. 2018
 *
 * @param[in,out] K_ptr Pointer to a number. Upon entry it should contain the
 *      size of the array roots. Upon exit, it contains the number of found
 *      roots.
 * @param [out] roots Pre-allocated array in which the routine will store the
 *      found roots.
 * @param [in] fun Pointer to a routine that evaluates the function whose roots
 *      are to be found. The first argument is number of desired function
 *      evaluations, the second argument contains the desired arguments for
 *      the evaluations, and the third argument is where the function values
 *      will be stored. The last argument may contain parameters that will
 *      passed through to the routine when it is called.
 * @param [in,out] params_ptr Pointer to parameters that we be provided to
 *      the routine fun everytime it is called.
 * @param [in] NodesMax TODO
 * @param [in] bounding_box_local Specifies the retangular region that is
 *      searched for roots. The first two entries speficy the smallest and
 *      largest real part, the last two entries the smallest and largest
 *      real part.
 * @param [in] tol Desired bound on the root approximation error, used to
 *      stop the algorithm.
 * @param [in] niter Maximum number of iterations.
 */
FNFT_INT fnft__global_root_pole_finding_algorithm(UINT * const K_ptr,
        COMPLEX * roots, INT fun (UINT, COMPLEX *, COMPLEX *, void *),
        void * params_ptr, UINT NodesMax, const REAL bounding_box_local[4],
        const REAL Tol, const UINT niter);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define global_root_pole_finding_algorithm(...) fnft__global_root_pole_finding_algorithm(__VA_ARGS__)
#endif

#endif
