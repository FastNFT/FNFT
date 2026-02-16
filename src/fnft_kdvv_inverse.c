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

static INT find_first_positive_value(
    COMPLEX const * const array,
    const UINT length,
    UINT * const i_pos_ptr)
{
    INT ret_code = SUCCESS;
    
    for (UINT i = 0; i < length; i++){
        if (CREAL(array[i]) > 0){
            *i_pos_ptr = i;
            break;
        }
    }

    return ret_code;
}

static INT one_solition_crum_transformation(
    COMPLEX const k1,
    COMPLEX const * const th1,
    COMPLEX const * const th2,
    const UINT D,
    const UINT i_pos,
    COMPLEX const * const x_grid,
    COMPLEX * const w_inv,
    COMPLEX * const prefactor,
    COMPLEX * const M_min1_11)
{
    INT ret_code = SUCCESS;

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

    return ret_code;
}

static INT two_solitions_crum_transformation(
    COMPLEX const k1,
    COMPLEX const k2,
    const UINT N_pm0_jZ,
    COMPLEX const * const pm0_jZ,
    COMPLEX const * const theta_E1,
    COMPLEX const * const theta_E2,
    COMPLEX * const t1,
    COMPLEX * const t2,
    const UINT D,
    const UINT i_pos,
    COMPLEX const * const x_grid,
    COMPLEX * const w_inv,
    COMPLEX * const prefactor,
    COMPLEX * const dq,
    COMPLEX * const s)
{
    INT ret_code = SUCCESS;

    // Defining theta vectors out of theta_E vectors, belonging to the bound states to add
    // These vectors are used in the Crum-Transformation step
    // thetas according to 1. eigenvalue
    COMPLEX const * const th11 = &theta_E1[0];
    COMPLEX const * const th12 = &theta_E2[0];

    // thetas according to 2. eigenvalue
    COMPLEX const * const th21 = &theta_E1[D];
    COMPLEX const * const th22 = &theta_E2[D];

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

leave_fun:

    return ret_code;
}



static INT one_solition_update_Jost(
    const COMPLEX k1,
    const UINT N_jZ,
    COMPLEX const * const bound_states_left,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    const UINT D,
    const UINT i_pos,
    COMPLEX const * const x_grid,
    COMPLEX const * const w_inv,
    COMPLEX const * const prefactor,
    COMPLEX const * const M_min1_11)
{
    INT ret_code = SUCCESS;

    // Check, if Jost has to be updated. Only when there are bound states left, that are not added yet
    UINT const is_Jost_to_update = N_jZ != 0;

    // Defining theta vectors out of theta_E vectors, belonging to the bound state to add
    COMPLEX const * const th1 = &theta_E1[0];
    COMPLEX const * const th2 = &theta_E2[0];

    // theta vectors according to the bound states left, which corresponds to the Jost solutions
    COMPLEX * const th1_jZ = &theta_E1[D];
    COMPLEX * const th2_jZ = &theta_E2[D];

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

    // Compute a new prefactor for C_E
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
    free(C_E);
    free(prefactor_C_E);

    return ret_code;
}

static INT two_solitions_update_Jost(
    const COMPLEX k1,
    const COMPLEX k2,
    const UINT N_jZ,
    COMPLEX const * const bound_states_left,
    COMPLEX const * const pm0_jZ,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    COMPLEX const * const t1,
    COMPLEX const * const t2,
    const UINT D,
    const UINT i_pos,
    COMPLEX const * const x_grid,
    COMPLEX const * const w_inv,
    COMPLEX const * const prefactor,
    COMPLEX const * const s)
{
    INT ret_code = SUCCESS;

    // Check, if Jost has to be updated. Only when there are bound states left, that are not added yet
    UINT const is_Jost_to_update = N_jZ != 0;

    // Defining theta vectors out of theta_E vectors, belonging to the bound states to add
    // These vectors are used in the Crum-Transformation step
    // thetas according to 1. eigenvalue
    COMPLEX const * const th11 = &theta_E1[0];
    COMPLEX const * const th12 = &theta_E2[0];

    // thetas according to 2. eigenvalue
    COMPLEX const * const th21 = &theta_E1[D];
    COMPLEX const * const th22 = &theta_E2[D];

    // theta vectors according to the bound states left, which corresponds to the Jost solutions
    COMPLEX * const th1_jZ = &theta_E1[2*D];
    COMPLEX * const th2_jZ = &theta_E2[2*D];

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
    COMPLEX mpjZinv_0;
    COMPLEX mpjZinv_1;
    COMPLEX mpjZinv_2;
    COMPLEX M_min1_11_jZinv;

    // Calculation for positive x
    for (UINT i=i_pos; i<D && is_Jost_to_update; i++){
        prefactor1 = (k2*k2 - k1*k1) * t1[i] * t2[i] * w_inv[i];
        M_0_11 = -(k1*k1 + k2*k2)/2 + prefactor[i] * (k2*k2 - k1*k1)/4 * t1[i]*t1[i] * t2[i]*t2[i];
        prefactor2 = prefactor[i] * k1 * k2;
        
        mpjZinv_0 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] - th12[i]*th12[i] * 
                    CEXP(-4*k1*x_grid[i])) * CEXP(-2*(k2) * x_grid[i]) + 
                    -k1 * th11[i] * th12[i] * (th21[i]*th21[i] - th22[i]*th22[i] * 
                    CEXP(-4*k2*x_grid[i])) * CEXP(-2*(k1) * x_grid[i]));

        for (UINT j=0; j<N_jZ; j++){
            mpjZinv_1 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] - th12[i]*th12[i] * 
                        CEXP(-4*k1*x_grid[i])) * CEXP(-2*(k2-pm0_jZ[j]) * x_grid[i]) + 
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] - th22[i]*th22[i] * 
                        CEXP(-4*k2*x_grid[i])) * CEXP(-2*(k1-pm0_jZ[j]) * x_grid[i])  
                    ) * 1/bound_states_left[j];

            mpjZinv_2 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] - th12[i]*th12[i] * 
                        CEXP(-4*k1*x_grid[i])) * CEXP(-2*(k2-pm0_jZ[j+N_jZ]) * x_grid[i]) + 
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] - th22[i]*th22[i] * 
                        CEXP(-4*k2*x_grid[i])) * CEXP(-2*(k1-pm0_jZ[j+N_jZ]) * x_grid[i])  
                    ) * 1/bound_states_left[j];

            M_min1_11_jZinv = prefactor2 * mpjZinv_0 * 1/bound_states_left[j];

            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor1 * bound_states_left[j] + M_0_11 + M_min1_11_jZinv;
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * s[j*D+i] + prefactor2 * mpjZinv_1;
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] + prefactor[i] * s[(j+N_jZ)*D+i] - prefactor2 * mpjZinv_2;
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] - prefactor1 * bound_states_left[j] + M_0_11 - M_min1_11_jZinv;
        }
    }

    // Calculation for negative x
    for (UINT i=0; i<i_pos && is_Jost_to_update; i++){
        prefactor1 = (k2*k2 - k1*k1) * t1[i] * t2[i] * w_inv[i];
        M_0_11 = -(k1*k1 + k2*k2)/2 + prefactor[i] * (k2*k2 - k1*k1)/4 * t1[i]*t1[i] * t2[i]*t2[i];
        prefactor2 = prefactor[i] * k1 * k2;

        mpjZinv_0 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] * CEXP(4*k1*x_grid[i]) - th12[i]*th12[i]) *               
                    CEXP(2*(k2) * x_grid[i]) +
                    -k1 * th11[i] * th12[i] * (th21[i]*th21[i] * CEXP(4*k2*x_grid[i]) - th22[i]*th22[i]) *
                    CEXP(2*(k1) * x_grid[i]) );

        for (UINT j=0; j<N_jZ; j++){
            mpjZinv_1 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] * CEXP(4*k1*x_grid[i]) - 
                        th12[i]*th12[i]) * CEXP(2*(k2+pm0_jZ[j]) * x_grid[i]) +
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] * CEXP(4*k2*x_grid[i]) - 
                        th22[i]*th22[i]) * CEXP(2*(k1+pm0_jZ[j]) * x_grid[i])
                    ) * 1/bound_states_left[j];

            mpjZinv_2 = (k2 * th21[i] * th22[i] * (th11[i]*th11[i] * CEXP(4*k1*x_grid[i]) - 
                        th12[i]*th12[i]) * CEXP(2*(k2+pm0_jZ[j+N_jZ]) * x_grid[i]) +
                        -k1 * th11[i] * th12[i] * (th21[i]*th21[i] * CEXP(4*k2*x_grid[i]) - 
                        th22[i]*th22[i]) * CEXP(2*(k1+pm0_jZ[j+N_jZ]) * x_grid[i])
                    ) * 1/bound_states_left[j];

            M_min1_11_jZinv = prefactor2 * mpjZinv_0 * 1/bound_states_left[j];

            C_E_1_1[j*D+i] = C_E_1_1[j*D+i] + prefactor1 * bound_states_left[j] + M_0_11 + M_min1_11_jZinv;
            C_E_1_2[j*D+i] = C_E_1_2[j*D+i] + prefactor[i] * s[j*D+i] + prefactor2 * mpjZinv_1;
            C_E_2_1[j*D+i] = C_E_2_1[j*D+i] + prefactor[i] * s[(j+N_jZ)*D+i] - prefactor2 * mpjZinv_2;
            C_E_2_2[j*D+i] = C_E_2_2[j*D+i] - prefactor1 * bound_states_left[j] + M_0_11 - M_min1_11_jZinv;
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
    free(C_E);

    return ret_code;
}


static INT add_one_soliton(
    COMPLEX const * const bound_states,
    UINT const N_bound_states,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    COMPLEX const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    // Initialize variables and arrays
    INT ret_code = SUCCESS;

    COMPLEX * M_min1_11 = NULL;
    COMPLEX * w_inv = NULL;
    COMPLEX * prefactor = NULL;
    M_min1_11 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(M_min1_11, ret_code, leave_fun);
    w_inv = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(w_inv, ret_code, leave_fun);
    prefactor = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor, ret_code, leave_fun);

    // Determine where x>0
    UINT i_pos;
    ret_code = find_first_positive_value(x_grid, D, &i_pos);
    CHECK_RETCODE(ret_code, leave_fun);

    // -- Crum Transform --
    COMPLEX bound_state_added_here = bound_states[0];
    ret_code = one_solition_crum_transformation(    bound_state_added_here,
                                                    &theta_E1[0],           // thetas belonging to the bound state to add
                                                    &theta_E2[0],           // thetas belonging to the bound state to add
                                                    D,
                                                    i_pos,
                                                    x_grid,
                                                    w_inv,
                                                    prefactor,
                                                    M_min1_11);
    CHECK_RETCODE(ret_code, leave_fun);                                                
    
    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = -q[i] + 4 * M_min1_11[i];      
    }
    
    // Update Jost Solution only when there are further bound states left, which are not added yet
    UINT N_bound_states_left = N_bound_states - 1; 
    if (N_bound_states_left) {
        COMPLEX const * const bound_states_left = &bound_states[1];
        ret_code = one_solition_update_Jost(    bound_state_added_here,
                                                N_bound_states_left,        
                                                bound_states_left,
                                                theta_E1,
                                                theta_E2,
                                                D,
                                                i_pos,
                                                x_grid,
                                                w_inv,
                                                prefactor,
                                                M_min1_11);
        CHECK_RETCODE(ret_code, leave_fun);                                            
    }

leave_fun:
    free(M_min1_11);
    free(w_inv);
    free(prefactor);

    return ret_code;
}


INT add_two_solitons(
    COMPLEX const * const bound_states,
    UINT const N_bound_states,
    COMPLEX * const theta_E1,
    COMPLEX * const theta_E2,
    COMPLEX const * const x_grid,
    const UINT D,
    COMPLEX * const q)
{
    INT ret_code = SUCCESS;

    COMPLEX * w_inv = NULL;
    COMPLEX * prefactor = NULL;
    COMPLEX * dq = NULL;
    w_inv = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(w_inv, ret_code, leave_fun);
    prefactor = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(prefactor, ret_code, leave_fun);
    dq = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(dq, ret_code, leave_fun);

    // For the bound states, which are not 
    const UINT N_bound_states_left = N_bound_states - 2;
    UINT const is_Jost_to_update = N_bound_states_left != 0;

    COMPLEX const * bound_states_left = NULL;

    if (is_Jost_to_update) {
        bound_states_left = &bound_states[2];
    }

    // Initialize pm0_jZ
    UINT N_pm0_jZ = 2 * N_bound_states_left + 1;
    COMPLEX * pm0_jZ = NULL;
    pm0_jZ = malloc((2*N_bound_states_left+1) * sizeof(COMPLEX));
    CHECK_NOMEM(pm0_jZ, ret_code, leave_fun);
    
    for (UINT i=0; i<N_bound_states_left && is_Jost_to_update; i++){
        pm0_jZ[i] = bound_states_left[i];
        pm0_jZ[i+N_bound_states_left] = -bound_states_left[i];
    }
    
    pm0_jZ[2*N_bound_states_left] = 0;

    // dimension of s is dependend on number of bound states left
    COMPLEX * s = NULL;
    s = malloc(D * N_pm0_jZ * sizeof(COMPLEX));
    CHECK_NOMEM(s, ret_code, leave_fun);

    // Intermediate variables
    COMPLEX * t1 = NULL;
    COMPLEX * t2 = NULL;

    t1 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t1, ret_code, leave_fun);
    t2 = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(t2, ret_code, leave_fun);
    
    // Change Sign of the bound states to add
    COMPLEX k1 = -bound_states[0];
    COMPLEX k2 = -bound_states[1];

    // Determine where x>0
    UINT i_pos;
    ret_code = find_first_positive_value(x_grid, D, &i_pos);
    CHECK_RETCODE(ret_code, leave_fun);

    // -- Crum-Transformation --
    ret_code = two_solitions_crum_transformation(   k1,
                                                    k2,
                                                    N_pm0_jZ,
                                                    pm0_jZ,
                                                    theta_E1,
                                                    theta_E2,
                                                    t1,
                                                    t2,
                                                    D,
                                                    i_pos,
                                                    x_grid,
                                                    w_inv,
                                                    prefactor,
                                                    dq,
                                                    s   );
    CHECK_RETCODE(ret_code, leave_fun);


    // Update output
    for (UINT i=0; i<D; i++){
        q[i] = q[i] + 4 * dq[i];      
    }

    // -- Update Jost Solution --
    if (is_Jost_to_update) {
        ret_code = two_solitions_update_Jost(   k1,
                                                k2,
                                                N_bound_states_left,
                                                bound_states_left,
                                                pm0_jZ,
                                                theta_E1,
                                                theta_E2,
                                                t1,
                                                t2,
                                                D,
                                                i_pos,
                                                x_grid,
                                                w_inv,
                                                prefactor,
                                                s   );
        CHECK_RETCODE(ret_code, leave_fun);
    }

leave_fun:
    free(w_inv);
    free(prefactor);
    free(dq);
    free(pm0_jZ);
    free(s);
    free(t1);
    free(t2);
    
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
            ret_code = add_one_soliton( &bound_states_sorted[step_idx], 
                                        N_rest,
                                        &theta_E1[step_idx*D],
                                        &theta_E2[step_idx*D],
                                        x_grid, 
                                        D, 
                                        q);
            CHECK_RETCODE(ret_code, leave_fun);
        } else {
            ret_code = add_two_solitons(&bound_states_sorted[step_idx], 
                                        N_rest,
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

