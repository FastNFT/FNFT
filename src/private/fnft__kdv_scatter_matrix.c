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
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J Prins (TU Delft) 2020.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__kdv_scatter.h"

/**
 * If derivative_flag=0 returns [S11 S12 S21 S22] in result where
 * S = [S11, S12; S21, S22] is the scattering matrix computed using the 
 * chosen scheme.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in 
 * result where S11' is the derivative of S11 w.r.t to lambda.
 * Result should be preallocated with size 4*K or 8*K accordingly.
 */
INT kdv_scatter_matrix(const UINT D, COMPLEX const * const q,
    COMPLEX const * const r, const REAL eps_t, const INT kappa,
    const UINT K, COMPLEX const * const lambda,
    COMPLEX * const result, kdv_discretization_t discretization,
    const UINT derivative_flag)
{
    (void) &kappa; // Suppress compiler warning
    INT ret_code = SUCCESS;
    UINT i;
    akns_discretization_t akns_discretization;

    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (!(eps_t > 0))
        return E_INVALID_ARGUMENT(eps_t);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (lambda == NULL)
        return E_INVALID_ARGUMENT(lambda);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    
    ret_code = kdv_discretization_to_akns_discretization(discretization,
            &akns_discretization);
    CHECK_RETCODE(ret_code, leave_fun);


    INT vanilla_flag;
    kdv_discretization_vanilla_flag(&vanilla_flag, discretization);

    // Calculate the change of state matrix in AKNS basis across the potential
    if(vanilla_flag) {
        ret_code = akns_scatter_matrix(D, q, r, eps_t, K, lambda, result,
            akns_discretization, derivative_flag);
        CHECK_RETCODE(ret_code, leave_fun);
    } else {
        ret_code = akns_scatter_matrix(D, r, q, eps_t, K, lambda, result,
        akns_discretization, derivative_flag);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Change the basis of the change of state matrix to the S-basis
    UINT const N = derivative_flag ? 4 : 2; // size of the matrices: 2x2 or 4x4
    COMPLEX Tmx[4][4], Mmx[4][4], H[4*4];
    UINT const H_idx_to_result_idx[4*4] = { 0,1,8,8,
                                            2,3,8,8,
                                            4,5,0,1,
                                            6,7,2,3 };
    COMPLEX const W[8][16] = {
        {0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0} };
    for (i=0; i<K; i++) { // Loop over lambda samples
        // Copy the values of the change of state matrix in AKNS basis to a 2 dimensional array
        if (!derivative_flag) {
            for (UINT n=0; n<N*N; n++){
                H[n] = result[i*4+n];
            }
        } else {
            for (UINT n=0; n<N*N; n++){
                if (H_idx_to_result_idx[n]==8){
                    H[n] = 0.0;
                } else {
                    H[n] = result[i*8+H_idx_to_result_idx[n]];
                }
            }
        }

        // Fetch the change of basis matrix from S to AKNS
        ret_code = kdv_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],lambda[i],derivative_flag,eps_t,discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        // Right-multiply the change of state matrix by the change of basis matrix from S to AKNS
        misc_matrix_mult(N,N,N,&H[0],&Tmx[0][0],&Mmx[0][0]);

        // Fetch the change of basis matrix from AKNS to S
        ret_code = kdv_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],lambda [i],derivative_flag,eps_t,discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        // Left-multiply the change of state matrix by the change of basis matrix from AKNS to S
        misc_matrix_mult(N,N,N,&Tmx[0][0],&Mmx[0][0],&H[0]);

        // Copy to the result
        if (!derivative_flag) {
            for (UINT n=0; n<4; n++){
                result[i*4+n] = H[n];
            }
        } else {
            for (UINT n=0; n<8; n++){
                // The line below looks in Matlab notation like result(8*i+n) = W(n,:) * H(:);
                misc_matrix_mult(1,16,1,&W[n][0],&H[0],&result[i*8+n]);
            }
        }
    }
    
    leave_fun:
        return ret_code;
}
