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
 * Sander Wahls (KIT) 2026
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_kdvv_inverse.h"
#include <stdio.h> // for printf


// INT input_multisoliton_E_pairs(
//     UINT const K,
//     COMPLEX const * const bound_states,
//     COMPLEX const * const normconsts_or_residues,
//     const UINT D,
//     COMPLEX * const q,
//     REAL const * const T)
// {
//     INT ret_code = SUCCESS;

//     // Initialize local variables
//     COMPLEX tmp;
//     COMPLEX * bnd_states = NULL;
//     COMPLEX * norm_consts = NULL;

//     bnd_states = malloc(K * sizeof(COMPLEX));
//     CHECK_NOMEM(bnd_states, ret_code, leave_fun);
//     norm_consts = malloc(K * sizeof(COMPLEX));
//     CHECK_NOMEM(norm_consts, ret_code, leave_fun);

//     //sorting bnd_states and according norm_const in descending order based on magnitude of imaginary part
//     for (UINT i = 0; i < K; ++i){
//         for (UINT j = i + 1; j < K; ++j){
//             if (CIMAG(bnd_states[i]) < CIMAG(bnd_states[j])){
//                 tmp =  bnd_states[i];
//                 bnd_states[i] = bnd_states[j];
//                 bnd_states[j] = tmp;
//                 tmp =  norm_consts[i];
//                 norm_consts[i] = norm_consts[j];
//                 norm_consts[j] = tmp;
//             }
//         }
//     }

    

//     // Apply normconst
//     //TODO

//     // Classic Crum transformation

//     UINT N_rest = K;
    
//     while(N_rest > 0){

//         // Scale trajectories
//         //TODO  

//         ret_code = add_two_solitons_E(q, x_grid, );
//         CHECK_RETCODE(ret_code, leave_fun);
//     }

// leave_fun:
//     free(bound_states);
//     free(norm_consts);

//     return ret_code;
// }


INT add_one_soliton_E(
    COMPLEX const * const bound_state,
    REAL const * const theta_E1,
    REAL const * const theta_E2,
    REAL const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    INT ret_code = SUCCESS;

    COMPLEX * M_min1_11 = NULL;
    M_min1_11 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(M_min1_11, ret_code, leave_fun);

    // Transformation into a real number.
    COMPLEX k1 = bound_state[0]*I;

    // printf("Marker1");

    // Determine where x>0
    UINT i_pos;
    
    for (UINT i = 0; i < D; i++){
        if (x_grid[i] > 0){
            i_pos = i;
            break;
        }
    }

    UINT n_pos=D-i_pos;
    UINT n_neg=i_pos;

    // Divide x_grid into positive and negative
    REAL * x_pos = NULL;
    REAL * x_neg = NULL;

    x_pos = malloc(n_pos * sizeof(REAL));
    CHECK_NOMEM(x_pos, ret_code, leave_fun);

    x_neg = malloc(n_neg * sizeof(REAL));
    CHECK_NOMEM(x_neg, ret_code, leave_fun);

    for (UINT i = 0; i < D; i++){
        if (i < i_pos){
            x_neg[i] = x_grid[i];
        } else {
            x_pos[i - i_pos] = x_grid[i];
        }
    }

    // divide theta_E
    COMPLEX * th1_pos = NULL;
    COMPLEX * th2_pos = NULL;
    COMPLEX * th1_neg = NULL;
    COMPLEX * th2_neg = NULL;

    th1_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th1_pos, ret_code, leave_fun);
    th2_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th2_pos, ret_code, leave_fun);
    th1_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th1_neg, ret_code, leave_fun);
    th2_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th2_neg, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        th1_neg[i] = theta_E1[i];
        th2_neg[i] = theta_E2[i];
    }

    for (UINT i=0; i<n_pos; i++){
        th1_pos[i] = theta_E1[i + i_pos];
        th2_pos[i] = theta_E2[i + i_pos];
    }

    // Calculation for positive x
    
    COMPLEX * w_p_inv = NULL;
    w_p_inv = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(w_p_inv, ret_code, leave_fun);

    for (UINT i=0; i<n_pos; i++){
        w_p_inv[i] = 1/(th1_pos[i] + th2_pos[i]*CEXP(2*k1*x_pos[i]));
    }

    COMPLEX * prefactor_pos = NULL;
    prefactor_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor_pos, ret_code, leave_fun);

    for (UINT i=0; i<n_pos; i++){
        prefactor_pos[i] = 2*k1*k1*th1_pos[i]*th2_pos[i]*w_p_inv[i]*w_p_inv[i];
    }

    for (UINT i=0; i<n_pos; i++){
        M_min1_11[i+n_neg]=prefactor_pos[i]*CEXP(2*k1*x_pos[i]);
    }
    

    // Calculation for negative x

    COMPLEX * w_n_inv = NULL;
    w_n_inv = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(w_n_inv, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        w_n_inv[i] = 1/(th1_neg[i]*CEXP(-2*k1*x_neg[i]) + th2_neg[i]);
    }

    COMPLEX * prefactor_neg = NULL;
    prefactor_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor_neg, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        prefactor_neg[i] = 2*k1*k1*th1_neg[i]*th2_neg[i]*w_n_inv[i]*w_n_inv[i];
    }

    for (UINT i=0; i<n_neg; i++){
        M_min1_11[i] = prefactor_neg[i]*CEXP(-2*k1*x_neg[i]);
    }

    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = -q[i] + 4 * M_min1_11[i];      
    }

leave_fun:
    free(x_pos);
    free(x_neg);
    free(th1_pos);
    free(th2_pos);
    free(th1_neg);
    free(th2_neg);
    free(w_p_inv);
    free(prefactor_pos);
    free(w_n_inv);
    free(prefactor_neg);

    return ret_code;
}


INT fnft_kdvv_inverse(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    UINT const K,
    COMPLEX const * const bound_states,
    COMPLEX const * const normconsts_or_residues,
    const UINT D,
    COMPLEX * const q,
    REAL const * const T,
    void *opts_ptr)
{
    
    /*
    D       number of spatial points
    out     output of transformation

    */

    /* Tests
    number of normconsts = number bound_states
    */
   
    INT ret_code = SUCCESS;

    // Initialize q
    for (UINT n=0; n<D; n++){
        q[n] = 0;
    }

    // Initialize spatial grid with D points between T[0] and T[1]
    REAL * x_grid = NULL;
    x_grid = malloc(D * sizeof(REAL));
    CHECK_NOMEM(x_grid, ret_code, leave_fun);

    const REAL eps_t = (T[1] - T[0])/(D - 1);

    for (UINT n=0; n<D; n++) {
        x_grid[n]= T[0] + n*eps_t;
    }

    // Initialize theta_E
    REAL * theta_E1 = NULL;
    theta_E1 = malloc(D * sizeof(REAL));
    CHECK_NOMEM(theta_E1, ret_code, leave_fun);

    REAL * theta_E2 = NULL;
    theta_E2 = malloc(D * sizeof(REAL));
    CHECK_NOMEM(theta_E2, ret_code, leave_fun);

    //TODO: only temporary solution!
    for (UINT i=0; i<D; i++){
        theta_E1[i]=1;
        theta_E2[i]=10;
    }

    // Add solitions
    for (UINT i=0; i<K; i++){
        ret_code = add_one_soliton_E(bound_states, theta_E1, theta_E2, x_grid, D, q);
    }



    // For Debugging:
    // for (UINT n=0; n<D; n++){
    //     q[n] = x_grid[n];
    // }
        

leave_fun:
    free(x_grid);
    free(theta_E1);
    free(theta_E2);
    return ret_code;
}

