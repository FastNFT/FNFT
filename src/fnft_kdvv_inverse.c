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


INT add_one_soliton_E(
    COMPLEX const * const bound_state,
    COMPLEX const * const theta_E1,
    COMPLEX const * const theta_E2,
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

    fnft_printf_ptr_t printf_ptr = fnft_errwarn_getprintf();
    printf_ptr("Marker! %d \n", ret_code);

    // Determine where x>0
    UINT i_pos;
    
    for (UINT i = 0; i < D; i++){
        if (x_grid[i] > 0){
            i_pos = i;
            break;
        }
    }

    COMPLEX * w_inv = NULL;
    w_inv = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(w_inv, ret_code, leave_fun);

    COMPLEX * prefactor = NULL;
    prefactor = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor, ret_code, leave_fun);
    

    // Calculation for positive x
    
    for (UINT i = i_pos; i<D; i++){
        w_inv[i] = 1/(theta_E1[i] + theta_E2[i]*CEXP(2*k1*x_grid[i]));
        prefactor[i] = 2*k1*k1*theta_E1[i]*theta_E2[i]*w_inv[i]*w_inv[i];
        M_min1_11[i]=prefactor[i]*CEXP(2*k1*x_grid[i]);

    }    

    // Calculation for negative x

    for (UINT i = 0; i<i_pos; i++){
        w_inv[i] = 1/(theta_E1[i]*CEXP(-2*k1*x_grid[i]) + theta_E2[i]);
        prefactor[i] = 2*k1*k1*theta_E1[i]*theta_E2[i]*w_inv[i]*w_inv[i];
        M_min1_11[i]=prefactor[i]*CEXP(-2*k1*x_grid[i]);

    }

    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = -q[i] + 4 * M_min1_11[i];      
    }

    // -- Update Jost Solution --
    COMPLEX * jZ0 = NULL; // Eigenvalues to add here
    COMPLEX * jZ = NULL; // Eigenvalues Left
    UINT N_jZ = sizeof(jZ);

    // Create C_E Matrix (for Updating Jost Solution)
    UINT n_C_E = 2*2*D*N_jZ; 
    COMPLEX * C_E = NULL;
    C_E = malloc(n_C_E*sizeof(COMPLEX));
    CHECK_NOMEM(C_E, ret_code, leave_fun);

    for (UINT i=0; i<n_C_E; i++){
        C_E[i] = 0;
    }

    COMPLEX * C_E_1_1 = &C_E[0*0*D*N_jZ];
    COMPLEX * C_E_2_2 = &C_E[1*1*D*N_jZ];

    for (UINT i=0; i<D; i++){
        for (UINT j=0; j<N_jZ; j++){
            C_E_1_1[i*N_jZ+j] = jZ[j];
            C_E_2_2[i*N_jZ+j] = -jZ[j];
        }
    }

    // Compute a newcprefactor for C_E (for Updating Jost Solution)
    // COMPLEX * prefactor_C_E_pos = NULL;
    // prefactor_C_E_pos = malloc(n_pos * sizeof(COMPLEX));
    // CHECK_NOMEM(prefactor_C_E_pos, ret_code, leave_fun);

    // COMPLEX * prefactor_C_E_neg = NULL;
    // prefactor_C_E_neg = malloc(n_neg * sizeof(COMPLEX));
    // CHECK_NOMEM(prefactor_C_E_neg, ret_code, leave_fun);

    // for (UINT i=0; i<n_pos; i++){
    //     prefactor_C_E_pos[i] = jZ0[0] * (th1_pos[i] - th2_pos[i] * CEXP(2*jZ0[0]*x_pos[i])) * w_p_inv[i];
    // }

    // for (UINT i=0; i<n_neg; i++){
    //     prefactor_C_E_neg[i] = jZ0[0] * (th1_neg[i] * CEXP(-2*jZ0[0]*x_neg[i]) - th2_neg[i]) * w_n_inv[i];
    // }

    // for (UINT i=i_pos; i<D; i++){
    //     for (UINT j=0; j<N_jZ; j++){
    //         C_E_1_1[i*N_jZ+j] = C_E_1_1[i*N_jZ+j] + prefactor_C_E_pos[i-i_pos] + M_min1_11;
    //         C_E_2_2[i*N_jZ+j] = C_E_2_2[i*N_jZ+j] + prefactor_C_E_pos[i-i_pos];
    //     }
    // }







leave_fun:
    free(w_inv);
    free(prefactor);
    free(C_E);
    // free(prefactor_C_E_pos);
    // free(prefactor_C_E_neg);

    return ret_code;
}


INT add_two_solitons_E(
    COMPLEX const * const bound_states,
    COMPLEX const * const theta_E11,
    COMPLEX const * const theta_E12,
    COMPLEX const * const theta_E21,
    COMPLEX const * const theta_E22,
    REAL const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    INT ret_code = SUCCESS;

    COMPLEX * M_min1_11 = NULL;
    M_min1_11 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(M_min1_11, ret_code, leave_fun);

    // Transformation into a real number.
    COMPLEX k1 = bound_states[0]*I;
    COMPLEX k2 = bound_states[1]*I;

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

    // theta_E for 1.eigenvalue
    COMPLEX * th11_pos = NULL;
    COMPLEX * th12_pos = NULL;
    COMPLEX * th11_neg = NULL;
    COMPLEX * th12_neg = NULL;

    th11_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th11_pos, ret_code, leave_fun);
    th12_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th12_pos, ret_code, leave_fun);
    th11_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th11_neg, ret_code, leave_fun);
    th12_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th12_neg, ret_code, leave_fun);

    // theta_E for 2.eigenvalue
    COMPLEX * th21_pos = NULL;
    COMPLEX * th22_pos = NULL;
    COMPLEX * th21_neg = NULL;
    COMPLEX * th22_neg = NULL;

    th21_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th21_pos, ret_code, leave_fun);
    th22_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(th22_pos, ret_code, leave_fun);
    th22_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th22_neg, ret_code, leave_fun);
    th21_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(th21_neg, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        th11_neg[i] = theta_E11[i];
        th21_neg[i] = theta_E21[i];
        th12_neg[i] = theta_E12[i];
        th22_neg[i] = theta_E22[i];
    }

    for (UINT i=0; i<n_pos; i++){
        th11_pos[i] = theta_E11[i + i_pos];
        th21_pos[i] = theta_E21[i + i_pos];
        th12_pos[i] = theta_E12[i + i_pos];
        th22_pos[i] = theta_E22[i + i_pos];
    }

    // Intermediate variables
    COMPLEX * t1_pos = NULL;
    COMPLEX * t2_pos = NULL;
    COMPLEX * t1_neg = NULL;
    COMPLEX * t2_neg = NULL;

    t1_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(t1_pos, ret_code, leave_fun);
    t2_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(t2_pos, ret_code, leave_fun);
    t1_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(t1_neg, ret_code, leave_fun);
    t2_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(t2_neg, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        t1_pos[i] = th11_pos[i] + th12_pos[i] * CEXP(-2*k1*x_pos[i]);
        t2_pos[i] = th21_pos[i] + th22_pos[i] * CEXP(-2*k2*x_pos[i]);
    }

    for (UINT i=0; i<n_pos; i++){
        t1_neg[i] = th11_neg[i] + th12_neg[i] * CEXP(2*k1*x_neg[i]);
        t2_neg[i] = th21_neg[i] + th22_neg[i] * CEXP(2*k2*x_neg[i]);
    }


    // Calculation for positive x
    
    COMPLEX * w_p_inv = NULL;
    w_p_inv = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(w_p_inv, ret_code, leave_fun);

    for (UINT i=0; i<n_pos; i++){
        w_p_inv[i] = 1/(t1_pos[i] * k2 * (th21_pos[i] - th22_pos[i]*CEXP(-2*k2*x_pos[i])) - t2_pos[i]*k1*(th11_pos[i]-th12_pos[i]*CEXP(-2*k1*x_pos[i])));
    }

    COMPLEX * prefactor_pos = NULL;
    prefactor_pos = malloc(n_pos * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor_pos, ret_code, leave_fun);

    for (UINT i=0; i<n_pos; i++){
        prefactor_pos[i] = 2*k1*k1*th11_pos[i]*th21_pos[i]*w_p_inv[i]*w_p_inv[i];
    }

    for (UINT i=0; i<n_pos; i++){
        M_min1_11[i+n_neg]=prefactor_pos[i]*CEXP(2*k1*x_pos[i]);
    }
    

    // Calculation for negative x

    COMPLEX * w_n_inv = NULL;
    w_n_inv = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(w_n_inv, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        w_n_inv[i] = 1/(th11_neg[i]*CEXP(-2*k1*x_neg[i]) + th21_neg[i]);
    }

    COMPLEX * prefactor_neg = NULL;
    prefactor_neg = malloc(n_neg * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor_neg, ret_code, leave_fun);

    for (UINT i=0; i<n_neg; i++){
        prefactor_neg[i] = 2*k1*k1*th11_neg[i]*th21_neg[i]*w_n_inv[i]*w_n_inv[i];
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
    free(th11_pos);
    free(th21_pos);
    free(th11_neg);
    free(th21_neg);
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

    // Initialize new storage for bound states and normconsts for later sorting
    COMPLEX tmp;

    COMPLEX * bound_states_sorted = NULL;
    bound_states_sorted = malloc(K * sizeof(COMPLEX));
    CHECK_NOMEM(bound_states_sorted, ret_code, leave_fun);

    COMPLEX * normconsts_sorted = NULL;
    normconsts_sorted = malloc(K * sizeof(COMPLEX));
    CHECK_NOMEM(normconsts_sorted, ret_code, leave_fun);

    for (UINT i=0; i<K; i++){
        bound_states_sorted[i] = bound_states[i];
        normconsts_sorted[i] = normconsts_or_residues[i];
    }

    //sorting bound_states and according norm_const in descending order based on magnitude of imaginary part
    for (UINT i = 0; i < K; ++i){
        for (UINT j = i + 1; j < K; ++j){
            if (CIMAG(bound_states_sorted[i]) < CIMAG(bound_states_sorted[j])){
                tmp =  bound_states_sorted[i];
                bound_states_sorted[i] = bound_states_sorted[j];
                bound_states_sorted[j] = tmp;
                tmp =  normconsts_sorted[i];
                normconsts_sorted[i] = normconsts_sorted[j];
                normconsts_sorted[j] = tmp;
            }
        }
    }

    // Shift Sign of norm_consts if K even
    if (K%2 == 0){
        for (UINT i=0; i<K; i++){
            normconsts_sorted[i] = -1*normconsts_sorted[i];
        }
    }    

    // Declare theta_E
    COMPLEX * theta_E11 = NULL;
    theta_E11 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E11, ret_code, leave_fun);

    COMPLEX * theta_E12 = NULL;
    theta_E12 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E12, ret_code, leave_fun);

    COMPLEX * theta_E21 = NULL;
    theta_E21 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E21, ret_code, leave_fun);

    COMPLEX * theta_E22 = NULL;
    theta_E22 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E22, ret_code, leave_fun);

    UINT N_rest = K;
    UINT step_idx = 0;

    // Add solitions
    while (N_rest > 0){

        // If the number of eigenvalues left is odd, then only one solition 
        // should be added
        if (N_rest%2==1){

            for (UINT j=0; j<D; j++){
                theta_E11[j] = 1;
                theta_E12[j] = normconsts_sorted[step_idx];
            }

            // Crum-Transformation step
            ret_code = add_one_soliton_E(&bound_states_sorted[step_idx], theta_E11, theta_E12, x_grid, D, q);
            CHECK_RETCODE(ret_code, leave_fun);

            step_idx++;
            N_rest--;
        } 
        // If the number of eigenvalues left is even, then two solitions
        // should be added        
        else {
            
        }

    }

    // for (UINT i=0; i<K; i++){

    //     // Assign values to theta_E dependend on next eigenvalue to add
    //     for (UINT j=0; j<D; j++){
    //         theta_E1_1[j] = 1;
    //         theta_E2_1[j] = normconsts_sorted[i];
    //     }

    //     // Crum-Transformation step
    //     ret_code = add_one_soliton_E(&bound_states_sorted[i], theta_E1_1, theta_E2_1, x_grid, D, q);
    //     CHECK_RETCODE(ret_code, leave_fun);
    // }


leave_fun:
    free(x_grid);
    free(bound_states_sorted);
    free(normconsts_sorted);
    free(theta_E11);
    free(theta_E12);
    return ret_code;
}

