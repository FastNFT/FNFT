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
* Peter J Prins (TU Delft) 2017-2018.
* Marius Brehler (TU Dortmund) 2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include <string.h> // for memcpy
#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_chirpz.h"
#include "fnft__kdv_fscatter.h"
#include "fnft__kdv_discretization.h"
#include "fnft_kdvv.h"

/**
 * Stores additional options for the routine fnft_kdvv.
 */
kdvv_opts_t default_opts = {
    .discretization = kdv_discretization_2SPLIT8B
};

/**
 * Creates a new options variable for fnft_kdvv with default settings.
 */
kdvv_opts_t fnft_kdvv_default_opts()
{
    return default_opts;
}

/**
 * Declare auxiliary routine used by the main routine fnft_kdvv.
 * Its body follows below.
 */
static INT tf2contspec_negxi(UINT deg,
    COMPLEX *transfer_matrix, REAL const * const T,
    const UINT D, REAL const * const XI, const UINT M,
    COMPLEX * result, fnft_kdvv_opts_t * opts_ptr);

/**
 * Fast nonlinear Fourier transform for the Korteweg-de Vries equation with
 * vanishing boundary conditions.
 */
INT fnft_kdvv(const UINT D,
    COMPLEX * const u,
    REAL const * const T,
    const UINT M, 
    COMPLEX * const contspec,
    REAL const * const XI,
    UINT * const K_ptr,
    COMPLEX * const bound_states,
    COMPLEX * const normconsts_or_residues,
    fnft_kdvv_opts_t * opts_ptr) 
{
    COMPLEX *transfer_matrix = NULL;
    UINT deg;
    INT ret_code = SUCCESS;
    INT W = 0, *W_ptr = NULL;
    (void) W; //To prevent complier warnings

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (u == NULL)
        return E_INVALID_ARGUMENT(u);
    if (T == NULL || T[0] >= T[1])
        return E_INVALID_ARGUMENT(T);
    if (contspec == NULL)
        return E_INVALID_ARGUMENT(contspec);
    if (XI == NULL || XI[0] >= XI[1])
        return E_INVALID_ARGUMENT(XI);
    if (K_ptr != NULL)
        return E_NOT_YET_IMPLEMENTED(K_ptr, Please pass "NULL".);
    if (bound_states != NULL)
        return E_NOT_YET_IMPLEMENTED(bound_states, Please pass "NULL".);
    if (normconsts_or_residues != NULL)
        return E_NOT_YET_IMPLEMENTED(normconsts_or_residues, Please pass "NULL".);

    if (opts_ptr == NULL)
        opts_ptr = &default_opts;

    // Allocate memory for the transfer matrix
    transfer_matrix = malloc(kdv_fscatter_numel(D,opts_ptr->discretization)*sizeof(COMPLEX));
    if (transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Determine step size
    const REAL eps_t = (T[1] - T[0])/(D - 1);

    // Compute the transfer matrix 
    ret_code = kdv_fscatter(D, u, eps_t, transfer_matrix, &deg,
        W_ptr, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);


    // Compute the continuous spectrum
    ret_code = tf2contspec_negxi(deg, transfer_matrix, T, D, XI, M,
                                     contspec, opts_ptr);
    CHECK_RETCODE(ret_code, release_mem);

release_mem:
    free(transfer_matrix);

    return ret_code;
}

// Auxiliary funnction: Computes continuous spectrum on a frequency grid
// from a given transfer matrix, using the negative frequencies.
static INT tf2contspec_negxi(UINT deg,
                                 COMPLEX *transfer_matrix, REAL const * const T,
                                 const UINT D, REAL const * const XI, const UINT M,
                                 COMPLEX *result, fnft_kdvv_opts_t * opts_ptr)
{
    COMPLEX *H_vals = NULL;
    COMPLEX *H11_vals = NULL;
    COMPLEX *H12_vals = NULL;
    COMPLEX *H21_vals = NULL;
    COMPLEX *H22_vals = NULL;
    COMPLEX A, V, sqrt_z;
    REAL boundary_coeff, degree1step;
    REAL xi;
    UINT i;
    INT ret_code = SUCCESS;

    if (opts_ptr == NULL)
        return E_INVALID_ARGUMENT(opts_ptr);
    
    // Determine discretization-specific coefficients
    degree1step = kdv_discretization_degree(opts_ptr->discretization);
    boundary_coeff = kdv_discretization_boundary_coeff(opts_ptr->discretization);
    if (degree1step == NAN || boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);

    // Allocate memory
    H_vals = malloc(4*M * sizeof(COMPLEX));
    if (H_vals == NULL)
        return E_NOMEM;
    H11_vals = H_vals;
    H12_vals = H11_vals + M;
    H21_vals = H12_vals + M;
    H22_vals = H21_vals + M;
    
    // Set step sizes
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    const REAL eps_xi = (XI[1] - XI[0])/(M - 1);
    
    // Evaluate the entries of the transfer matrix on the frequency grid
    // xi(i) = -(XI1 + i*eps_xi), where i=0,...,M-1.
    // Since z=exp(2*j*xi*eps_t/degree1step), we find that the z at which z the
    // transfer matrix has to be evaluated are given by z(i)=1/(A*V^-i), where:
       
    V = CEXP(-2.0*I*eps_xi * eps_t / degree1step );
    A = CEXP( 2.0*I*XI[0] * eps_t / degree1step );
   
//    ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, H11_vals);
//    CHECK_RETCODE(ret_code, release_mem);

    ret_code = poly_chirpz(deg, transfer_matrix + (deg+1), A, V, M, H12_vals);
    CHECK_RETCODE(ret_code, release_mem);

//    ret_code = poly_chirpz(deg, transfer_matrix + 2*(deg+1), A, V, M,
//                           H21_vals);
//    CHECK_RETCODE(ret_code, release_mem);

    ret_code = poly_chirpz(deg, transfer_matrix + 3*(deg+1), A, V, M,
                           H22_vals);
    CHECK_RETCODE(ret_code, release_mem);

    if (opts_ptr->discretization==kdv_discretization_2SPLIT2A){
        // Correct H12_vals and H21_vals for trick that implements 2split2A with
        // first order polynomials instead of second order polynomials.
        for (i=0; i<M; i++) {
            xi = -XI[0] - i*eps_xi;
            sqrt_z = CEXP( I*xi*eps_t / degree1step);
            H12_vals[i] /= sqrt_z;
            H21_vals[i] *= sqrt_z;
        }
    }
    
    // Compute the continuous spectrum
    for (i=0; i<M; i++) {
        xi = -XI[0] - i*eps_xi;
        result[i] = CEXP( 2.0*I*xi * (T[1] + boundary_coeff*eps_t) )
            * H12_vals[i];
        result[i] /= 2.0*I*xi * H22_vals[i] - H12_vals[i];
    }
    
    // Release memory and return
release_mem:
    free(H_vals);
    return ret_code;
}
