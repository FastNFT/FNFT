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

INT add_one_soliton(
    COMPLEX const * const bound_state_to_add,
    COMPLEX const * const bound_states_left,
    UINT const N_bound_states_left,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    COMPLEX const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    // Initialize variables and arrays
    INT ret_code = SUCCESS;

    COMPLEX * M_min1_11 = NULL;
    M_min1_11 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(M_min1_11, ret_code, leave_fun);

    COMPLEX * w_inv = NULL;
    w_inv = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(w_inv, ret_code, leave_fun);

    COMPLEX * prefactor = NULL;
    prefactor = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor, ret_code, leave_fun);

    // Defining theta vectors out of theta_E vectors, belonging to the bound state to add
    // These vectors are used in the Crum-Transformation step
    COMPLEX * const th1 = &theta_E1[0];
    COMPLEX * const th2 = &theta_E2[0];

    // Transformation into a real number
    COMPLEX k1 = bound_state_to_add[0];

    // Determine where x>0
    UINT i_pos;
    
    for (UINT i = 0; i < D; i++){
        if (CREAL(x_grid[i]) > 0){
            i_pos = i;
            break;
        }
    }

    // -- Crum Transform --
    // Calculation for positive x
    for (UINT i = i_pos; i<D; i++){
        w_inv[i] = 1/(th1[i] + th2[i] * CEXP(2*k1*x_grid[i]));
        prefactor[i] = 2 * k1*k1 * th1[i] * th2[i] * w_inv[i]*w_inv[i];
        M_min1_11[i] = prefactor[i] * CEXP(2*k1*x_grid[i]);
    }    

    // Calculation for negative x
    for (UINT i = 0; i<i_pos; i++){
        w_inv[i] = 1/(th1[i] * CEXP(-2*k1*x_grid[i]) + th2[i]);
        prefactor[i] = 2 * k1*k1 * th1[i] * th2[i] * w_inv[i]*w_inv[i];
        M_min1_11[i] = prefactor[i] *  CEXP(-2*k1*x_grid[i]);

    }
    
    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = -q[i] + 4 * M_min1_11[i];      
    }
    
    // -- Update Jost Solution --
    // Check, if Jost has to be updated. Only when there are bound states left, that are not added yet
    UINT const N_jZ = N_bound_states_left;        // Only because shorter name
    UINT const is_Jost_to_update = N_jZ != 0;

    // theta vectors according to the bound states left, which corresponds to the Jost solutions
    COMPLEX * th1_jZ = &theta_E1[D];
    COMPLEX * th2_jZ = &theta_E2[D];

    // Create C_E Matrix (for Updating Jost Solution)
    UINT n_C_E = 2*2*D*N_jZ; 
    COMPLEX * C_E = NULL;
    C_E = malloc(n_C_E*sizeof(COMPLEX));
    CHECK_NOMEM(C_E, ret_code, leave_fun);

    for (UINT i=0; i<n_C_E; i++){
        C_E[i] = 0;
    }

    COMPLEX * C_E_1_1 = &C_E[0*D*N_jZ];
    COMPLEX * C_E_1_2 = &C_E[1*D*N_jZ];
    COMPLEX * C_E_2_1 = &C_E[2*D*N_jZ];
    COMPLEX * C_E_2_2 = &C_E[3*D*N_jZ];

    for (UINT i=0; i<D; i++){
        for (UINT j=0; j<N_jZ; j++){
            C_E_1_1[j*D+i] = -bound_states_left[j];
            C_E_2_2[j*D+i] = bound_states_left[j];
        }
    }

    // Compute a new prefactor for C_E (for Updating Jost Solution)
    COMPLEX * prefactor_C_E = NULL;
    prefactor_C_E = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor_C_E, ret_code, leave_fun);

    // Calculation for positive x
    for (UINT i=i_pos; i<D && is_Jost_to_update; i++){
        prefactor_C_E[i] = k1 * (th1[i] - th2[i] * CEXP(2*k1*x_grid[i])) * w_inv[i];

        for (UINT j=0; j<N_jZ; j++){
            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor_C_E[i] + M_min1_11[i] * 1/bound_states_left[j];
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * CEXP(2*(k1+bound_states_left[j])*x_grid[i])*1/bound_states_left[j];
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] - prefactor[i] * CEXP(2*(k1-bound_states_left[j])*x_grid[i])*1/bound_states_left[j];
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] + prefactor_C_E[i] - M_min1_11[i] * 1/bound_states_left[j];
        }
    }
    
    // Calculation for negative x
    for (UINT i=0; i<i_pos && is_Jost_to_update; i++){
        prefactor_C_E[i] = k1 * (th1[i] * CEXP(-2*k1*x_grid[i]) - th2[i]) * w_inv[i];
        
        for (UINT j=0; j<N_jZ; j++){
            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor_C_E[i] + M_min1_11[i] * 1/bound_states_left[j];
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * CEXP(-2*(k1-bound_states_left[j])*x_grid[i])*1/bound_states_left[j];
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] - prefactor[i] * CEXP(-2*(k1+bound_states_left[j])*x_grid[i])*1/bound_states_left[j];
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] + prefactor_C_E[i] - M_min1_11[i] * 1/bound_states_left[j];
        }
    }
    
    // Map Jost solution
    COMPLEX tmp_th1;
    COMPLEX tmp_th2;

    for (UINT i=0; i<D && is_Jost_to_update; i++){
        for (UINT j=0; j<N_jZ; j++){
            tmp_th1 = th1_jZ[j*D+i];
            tmp_th2 = th2_jZ[j*D+i];

            th1_jZ[j*D+i] = C_E_1_1[j*D+i] * tmp_th1 + C_E_1_2[j*D+i] * tmp_th2;
            th2_jZ[j*D+i] = C_E_2_1[j*D+i] * tmp_th1 + C_E_2_2[j*D+i] * tmp_th2;
        }
    }

leave_fun:
    free(M_min1_11);
    free(w_inv);
    free(prefactor);
    free(C_E);
    free(prefactor_C_E);

    return ret_code;
}


INT add_two_solitons(
    COMPLEX const * const bound_states_to_add,
    COMPLEX const * const bound_states_left,
    UINT const N_bound_states_left,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    COMPLEX const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    INT ret_code = SUCCESS;

    COMPLEX * w_inv = NULL;
    w_inv = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(w_inv, ret_code, leave_fun);

    COMPLEX * prefactor = NULL;
    prefactor = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor, ret_code, leave_fun);

    COMPLEX * dq = NULL;
    dq = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(dq, ret_code, leave_fun);

    // Initialize pm0_jZ
    UINT N_pm0_jZ = 2 * N_bound_states_left + 1;
    COMPLEX * pm0_jZ = NULL;
    pm0_jZ = malloc((2*N_bound_states_left+1) * sizeof(COMPLEX));
    CHECK_NOMEM(pm0_jZ, ret_code, leave_fun);
    
    for (UINT i=0; i<N_bound_states_left; i++){
        pm0_jZ[i] = -bound_states_left[i];
        pm0_jZ[2*i+1] = bound_states_left[i];
    }
    
    pm0_jZ[2*N_bound_states_left] = 0;
    
    // dimension of s is dependend on number of bound states left
    COMPLEX * s = NULL;
    s = malloc(D * N_pm0_jZ * sizeof(COMPLEX));
    CHECK_NOMEM(s, ret_code, leave_fun);

    // Change Sign of the bound states to add
    COMPLEX k1 = -bound_states_to_add[0];
    COMPLEX k2 = -bound_states_to_add[1];

    // Determine where x>0
    UINT i_pos;
    
    for (UINT i = 0; i < D; i++){
        if (CREAL(x_grid[i]) > 0){
            i_pos = i;
            break;
        }
    }

    // Defining theta vectors out of theta_E vectors, belonging to the bound states to add
    // These vectors are used in the Crum-Transformation step
    // 1. eigenvalue
    COMPLEX * th11 = &theta_E1[0];
    COMPLEX * th12 = &theta_E2[0];

    //2. eigenvalue
    COMPLEX * th21 = &theta_E1[D];
    COMPLEX * th22 = &theta_E2[D];

    // Intermediate variables
    COMPLEX * t1 = NULL;
    COMPLEX * t2 = NULL;

    t1 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t1, ret_code, leave_fun);
    t2 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t2, ret_code, leave_fun);

    // -- Crum-Transformation --
    // for positive x
    for (UINT i=i_pos; i<D; i++){
        t1[i] = th11[i] + th12[i] * CEXP(-2*k1*x_grid[i]);
        t2[i] = th21[i] + th22[i] * CEXP(-2*k2*x_grid[i]);
        w_inv[i] = 1/(t1[i] * k2 * (th21[i] - th22[i] * CEXP(-2*k2*x_grid[i])) - 
                    t2[i] * k1 * (th11[i] - th12[i] * CEXP(-2*k1*x_grid[i])));
        prefactor[i] = 2 * (k2*k2 - k1*k1) * w_inv[i]*w_inv[i];

        for (UINT j=0; j<N_pm0_jZ; j++) {
            s[j*D+i] =  t2[i]*t2[i] * k1*k1 * th11[i] * th12[i] * CEXP(-2*(k1-pm0_jZ[j])*x_grid[i]) -
                        t1[i]*t1[i] * k2*k2 * th21[i] * th22[i] * CEXP(-2*(k2-pm0_jZ[j])*x_grid[i]);
        }

        dq[i] = prefactor[i] * s[(N_pm0_jZ-1)*D+i];
    }

    // for negative x
    for (UINT i=0; i<i_pos; i++){
        t1[i] = th11[i] * CEXP(2*k1*x_grid[i]) + th12[i];
        t2[i] = th21[i] * CEXP(2*k2*x_grid[i]) + th22[i];
        w_inv[i] = 1/(t1[i] * k2 * (th21[i] * CEXP(2*k2*x_grid[i]) - th22[i]) - 
                    t2[i] * k1 * (th11[i] * CEXP(2*k1*x_grid[i]) - th12[i]));
        prefactor[i] = 2 * (k2*k2 - k1*k1) * w_inv[i]*w_inv[i];

        for (UINT j=0; j<N_pm0_jZ; j++) {
            s[j*D+i] =  t2[i]*t2[i] * k1*k1 * th11[i] * th12[i] * CEXP(2*(k1+pm0_jZ[j])*x_grid[i]) -
                        t1[i]*t1[i] * k2*k2 * th21[i] * th22[i] * CEXP(2*(k2+pm0_jZ[j])*x_grid[i]);
        }

        dq[i] = prefactor[i] * s[(N_pm0_jZ-1)*D+i];
    }

    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = q[i] + 4 * dq[i];      
    }


    // -- Update Jost Solution --

    // Check, if Jost has to be updated. Only when there are bound states left, that are not added yet
    UINT const N_jZ = N_bound_states_left;        // Only because shorter name
    UINT const is_Jost_to_update = N_jZ != 0;

    // theta vectors according to the bound states left, which corresponds to the Jost solutions
    COMPLEX * th1_jZ = &theta_E1[D];
    COMPLEX * th2_jZ = &theta_E2[D];

    // Create C_E Matrix (for Updating Jost Solution)
    UINT n_C_E = 2*2*D*N_jZ; 
    COMPLEX * C_E = NULL;
    C_E = malloc(n_C_E*sizeof(COMPLEX));
    CHECK_NOMEM(C_E, ret_code, leave_fun);

    for (UINT i=0; i<n_C_E; i++){
        C_E[i] = 0;
    }

    COMPLEX * C_E_1_1 = &C_E[0*D*N_jZ];
    COMPLEX * C_E_1_2 = &C_E[1*D*N_jZ];
    COMPLEX * C_E_2_1 = &C_E[2*D*N_jZ];
    COMPLEX * C_E_2_2 = &C_E[3*D*N_jZ];

    for (UINT i=0; i<D; i++){
        for (UINT j=0; j<N_jZ; j++){
            C_E_1_1[j*D+i] = bound_states_left[j]*bound_states_left[j];
            C_E_2_2[j*D+i] = bound_states_left[j]*bound_states_left[j];
        }
    }

    COMPLEX prefactor1;
    COMPLEX M_0_11;
    COMPLEX prefactor2;
    COMPLEX mpjZinv;
    COMPLEX mpjZinv_0;
    COMPLEX M_min1_11_jZinv;

    // Calculation for positive x
    for (UINT i=i_pos; i<D && is_Jost_to_update; i++){
        prefactor1 = (k2*k2 - k1*k1) * t1[i] * t2[i] * w_inv[i];
        M_0_11 = -(k1*k1 + k2*k2)/2 + prefactor[i] * (k2*k2 - k1*k1)/4 * t1[i]*t1[i] * t2[i]*t2[i];
        prefactor2 = prefactor[i] * k1 * k2;
        

        mpjZinv_0 = (  k2 * th21[i] * th22[i] * (th11[i]*th11[i] - th12[i]*th12[i] * 
                    CEXP(-4*k1*x_grid[i])) * CEXP(-2*(k2) * x_grid[i]) + 
                    -k1 * th11[i] * th12[i] * (th21[i]*th21[i] - th22[i]*th22[i] * 
                    CEXP(-4*k2*x_grid[i])) * CEXP(-2*(k1) * x_grid[i]));


        for (UINT j=0; j<N_jZ; j++){
            mpjZinv = (  k2 * th21[i] * th22[i] * (th11[i]*th11[i] - th12[i]*th12[i] * 
                        CEXP(-4*k1*x_grid[i])) * CEXP(-2*(k2-pm0_jZ[j]) * x_grid[i]) + 
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] - th22[i]*th22[i] * 
                        CEXP(-4*k2*x_grid[i])) * CEXP(-2*(k1-pm0_jZ[j]) * x_grid[i])  
                    ) * 1/bound_states_left[j];
            

            M_min1_11_jZinv = prefactor2 * mpjZinv_0 * 1/bound_states_left[j];

            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor1 * bound_states_left[j] + M_0_11 + M_min1_11_jZinv;
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * s[j*D+i] + prefactor2 * mpjZinv;
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] + prefactor[i] * s[(j+N_jZ)*D+i] - prefactor2 * mpjZinv;
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] - prefactor1 * bound_states_left[j] + M_0_11 - M_min1_11_jZinv;
        }
    }

    // Calculation for negative x
    for (UINT i=i_pos; i<D && is_Jost_to_update; i++){
        prefactor1 = (k2*k2 - k1*k1) * t1[i] * t2[i] * w_inv[i];
        M_0_11 = -(k1*k1 + k2*k2)/2 + prefactor[i] * (k2*k2 - k1*k1)/4 * t1[i]*t1[i] * t2[i]*t2[i];
        prefactor2 = prefactor[i] * k1 * k2;
        

        mpjZinv = (  k2 * th21[i] * th22[i] * (th11[i]*th11[i] * CEXP(4*k1*x_grid[i]) - th12[i]*th12[i]) *               
                    CEXP(2*(k2) * x_grid[i]) +
                    -k1 * th11[i] * th12[i] * (th21[i]*th21[i] * CEXP(4*k2*x_grid[i]) - th22[i]*th22[i]) *
                    CEXP(2*(k1) * x_grid[i]) );



        for (UINT j=0; j<N_jZ; j++){
            mpjZinv = (  k2 * th21[i] * th22[i] * (th11[i]*th11[i] * CEXP(4*k1*x_grid[i]) - 
                        th12[i]*th12[i]) * CEXP(2*(k2+pm0_jZ[j]) * x_grid[i]) +
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] * CEXP(4*k2*x_grid[i]) - 
                        th22[i]*th22[i]) * CEXP(2*(k1+pm0_jZ[j]) * x_grid[i])
                    ) * 1/bound_states_left[j];

            M_min1_11_jZinv = prefactor2 * mpjZinv_0 * 1/bound_states_left[j];

            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor1 * bound_states_left[j] + M_0_11 + M_min1_11_jZinv;
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * s[j*D+i] + prefactor2 * mpjZinv;
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] + prefactor[i] * s[(j+N_jZ)*D+i] - prefactor2 * mpjZinv;
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] - prefactor1 * bound_states_left[j] + M_0_11 - M_min1_11_jZinv;
        }
    }


    // Map Jost solution
    for (UINT i=0; i<D && is_Jost_to_update; i++){
        for (UINT j=0; j<N_jZ; j++){
            th1_jZ[j*D+i] = C_E_1_1[j*D+i] * th1_jZ[j*D+i] + C_E_1_2[j*D+i] * th2_jZ[j*D+i];
            th2_jZ[j*D+i] = C_E_2_1[j*D+i] * th1_jZ[j*D+i] + C_E_2_2[j*D+i] * th2_jZ[j*D+i];
        }
    }

leave_fun:
    free(w_inv);
    free(prefactor);
    free(dq);
    free(pm0_jZ);
    free(s);
    free(t1);
    free(t2);
    free(C_E);
    
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
    COMPLEX * x_grid = NULL;
    x_grid = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(x_grid, ret_code, leave_fun);

    const COMPLEX eps_t = (T[1] - T[0])/(D - 1);

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
        bound_states_sorted[i] = bound_states[i]*I; // bound states should be real for further computing
        normconsts_sorted[i] = normconsts_or_residues[i];
    }

    //sorting bound_states and according norm_const in descending order based on magnitude of imaginary part
    for (UINT i = 0; i < K; ++i){
        for (UINT j = i + 1; j < K; ++j){
            if (CABS(bound_states_sorted[i]) < CABS(bound_states_sorted[j])){
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
    COMPLEX * theta_E1 = NULL;
    theta_E1 = malloc(D * K * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E1, ret_code, leave_fun);

    COMPLEX * theta_E2 = NULL;
    theta_E2 = malloc(D * K * sizeof(COMPLEX));
    CHECK_NOMEM(theta_E2, ret_code, leave_fun);

    for (UINT i=0; i<D; i++){
        for (UINT j=0; j<K; j++){
            theta_E1[j*D+i] = 1;
            theta_E2[j*D+i] = normconsts_sorted[j];
        }
    }

    UINT N_rest = K;
    UINT step_idx = 0;
    UINT N_step;

    // Add solitions
    while (N_rest > 0){

        // Select indices for this step and for the remaining eigenvalues to add
        // If the number of eigenvalues left is odd, then only one solition should be added
        if (N_rest%2==1){ N_step = 1; }
        else { N_step = 2; }
        
        // Scale trajectories (magnitudes of th1 and th2 symmetric around 0)
        COMPLEX tmp_scaling_factor;

        for (UINT i=0; i<D; i++){
            for (UINT j=step_idx; j<step_idx+N_step; j++){
                tmp_scaling_factor = CPOW(2, -ROUND( (LOG2(CABS(theta_E1[j*D+i])) + LOG2(CABS(theta_E2[j*D+i])))/2 ));
                theta_E1[j*D+i] = theta_E1[j*D+i] * tmp_scaling_factor;
                theta_E2[j*D+i] = theta_E2[j*D+i] * tmp_scaling_factor;
            }
        }

        // Crum-Transformation step
        if (N_step==1){
            ret_code = add_one_soliton(&bound_states_sorted[step_idx], 
                                        &bound_states_sorted[step_idx+1],
                                        N_rest-1,
                                        &theta_E1[step_idx*D],
                                        &theta_E2[step_idx*D],
                                        x_grid, 
                                        D, 
                                        q);
            CHECK_RETCODE(ret_code, leave_fun);
        } else {
            ret_code = add_two_solitons(&bound_states_sorted[step_idx], 
                                        &bound_states_sorted[step_idx+2],
                                        N_rest-2,
                                        &theta_E1[step_idx*D],
                                        &theta_E2[step_idx*D],
                                        x_grid, 
                                        D, 
                                        q);
            CHECK_RETCODE(ret_code, leave_fun);
        }

        N_rest = N_rest-N_step;
        step_idx = step_idx + N_step;

        // Filter out errors that result in negative values
        for (UINT i=0; i<D; i++){
            if (CREAL(q[i]) < 0){ q[i] = 0; }
        }

    }

leave_fun:
    free(x_grid);
    free(bound_states_sorted);
    free(normconsts_sorted);
    free(theta_E1);
    free(theta_E2);
    return ret_code;
}

