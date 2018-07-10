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

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_nsev_inverse.h"
#include "fnft__errwarn.h"
#include "fnft__nse_discretization.h"
#include "fnft__nse_fscatter.h"
#include "fnft__poly_chirpz.h"
#include "fnft__nse_finvscatter.h"
#include "fnft__fft_wrapper.h"
#include "fnft__poly_specfact.h"

static fnft_nsev_inverse_opts_t default_opts = {
    .discretization = nse_discretization_2SPLIT2A,
    .contspec_type = fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT,
    .contspec_inversion_method = fnft_nsev_inverse_csmethod_DEFAULT,
    .max_iter = 100,
    .oversampling_factor = 8
};

fnft_nsev_inverse_opts_t fnft_nsev_inverse_default_opts()
{
    return default_opts;
}

INT fnft_nsev_inverse_XI(const UINT D, REAL const * const T,
                         const UINT M, REAL * const XI,
                         const nse_discretization_t discretization)
{
    if (M == 0)
        return E_INVALID_ARGUMENT(M);
    if (XI == NULL)
        return E_INVALID_ARGUMENT(XI);
    if (T == NULL || !(T[0] < T[1]))
        return E_INVALID_ARGUMENT(T);

    const REAL degree1step = nse_discretization_degree(discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(discretization);
    const REAL map_coeff = 2/degree1step;
    // Note that z = exp(2.0*PI*I*i/M) = exp(map_coeff*I*xi*eps_t).
    // The range below results the same values as a normal M-point FFT,
    // but in a different order.
    const REAL eps_t = (T[1] - T[0]) / (D - 1);
    XI[0] = CLOG( CEXP(2.0*FNFT_PI*I * (M/2 + 1)/M) ) / (map_coeff*I*eps_t);
    XI[1] = CLOG( CEXP(2.0*FNFT_PI*I * (M/2)/M) ) / (map_coeff*I*eps_t);
    return SUCCESS;
}

// Auxliary functions. Function bodies follow below.
static INT transfer_matrix_from_reflection_coefficient(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr);
static INT transfer_matrix_from_B_of_tau(
    const UINT M,
    COMPLEX * const contspec,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr);

INT fnft_nsev_inverse(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    UINT const K,
    COMPLEX const * const bound_states,
    COMPLEX const * const normconsts_or_residues,
    const UINT D,
    COMPLEX * const q,
    REAL const * const T,
    const INT kappa,
    fnft_nsev_inverse_opts_t *opts_ptr)
{
    if (M%2 != 0)
        return E_INVALID_ARGUMENT(M);
    if (M < D)
        return E_INVALID_ARGUMENT(M);
    if (contspec == NULL)
        return E_INVALID_ARGUMENT(contspec);
    if (K > 0)
        return E_NOT_YET_IMPLEMENTED(K>0, Please pass K=0.);
    if (bound_states != NULL)
        return E_NOT_YET_IMPLEMENTED(bound_states!=NULL, Please pass NULL.)
    if (normconsts_or_residues != NULL)
        return E_NOT_YET_IMPLEMENTED(normconsts_or_residues!=NULL, Please pass NULL.)
    if (D < 2 || (D&(D - 1)) != 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL || !(T[0] < T[1]))
        return E_INVALID_ARGUMENT(T);
    if (kappa != +1 && kappa != -1)
        return E_INVALID_ARGUMENT(kappa);
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->discretization != nse_discretization_2SPLIT2A
    && opts_ptr->discretization != nse_discretization_2SPLIT2_MODAL)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);

    INT ret_code = SUCCESS;

    // To avoid compiler warnings about unused variables
    (void) bound_states;
    (void) normconsts_or_residues;

    // Allocated memory and initialize values that we will need later

    const UINT deg = D*nse_discretization_degree(opts_ptr->discretization);
    if (deg == 0)
        return E_INVALID_ARGUMENT(discretization);

    COMPLEX * const transfer_matrix = malloc(4*(deg+1) * sizeof(COMPLEX));
    if (transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }



    // Step 1: Construct the transfer matrix from the provided representation
    // of the continuous spectrum.

    switch (opts_ptr->contspec_type) {

    case fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT:
        ret_code = transfer_matrix_from_reflection_coefficient(M, contspec, XI,
                                                               D, T, deg,
                                                               transfer_matrix,
                                                               kappa, opts_ptr);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    case fnft_nsev_inverse_cstype_B_OF_TAU:
        ret_code = transfer_matrix_from_B_of_tau(M, contspec, D, T, deg,
                                                 transfer_matrix, kappa,
                                                 opts_ptr);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    default:

        ret_code = E_INVALID_ARGUMENT(opts_ptr->contspec_type);
        goto leave_fun;

    }

    // Step 2: Recover the time domain signal from the transfer matrix.

    const REAL eps_t = (T[1] - T[0]) / (D - 1);
    ret_code = nse_finvscatter(deg, transfer_matrix, q, eps_t, kappa,
        opts_ptr->discretization);

leave_fun:
    free(transfer_matrix);
    return ret_code;
}

// This auxiliary function removes the phase factors that result from
// the boundary conditions from the continuous spectrum.
// It furthermore creates a reordered version of the corrected spectrum
// that is compatible with the ordering used by the FFT.
static inline INT remove_boundary_conds_and_reorder_for_fft(
    const UINT M,
    COMPLEX * const contspec,
    COMPLEX * const contspec_reordered,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const nse_discretization_t discretization)
{
    UINT i;
    INT ret_code = SUCCESS;
    const REAL eps_t = (T[1] - T[0]) / (D - 1);
    const REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    REAL phase_factor_rho = 0.0;
    const REAL bnd_coeff = nse_discretization_boundary_coeff(discretization);
    const REAL degree1step = nse_discretization_degree(discretization);
    if (degree1step == NAN || bnd_coeff == NAN)
        return E_INVALID_ARGUMENT(discretization);
    if (discretization == nse_discretization_2SPLIT2A || discretization == nse_discretization_2SPLIT2_MODAL)
        phase_factor_rho = -2.0*(T[1] + eps_t*bnd_coeff) + eps_t/degree1step;
    else
        phase_factor_rho = -2.0*(T[1] + eps_t*bnd_coeff);


    for (i=0; i<M; i++) {
        const REAL xi = XI[0] + i*eps_xi;
        contspec[i] *= CEXP(-I*xi*phase_factor_rho);
    }
    for (i=0; i<=M/2; i++)
        contspec_reordered[i] = contspec[i + (M/2 - 1)];
    for (i=M/2+1; i<M; i++)
        contspec_reordered[i] = contspec[i - (M/2 + 1)];
    
    return ret_code;
    
}

// This auxiliary function builds a transfer matrix from a given reflection
// coefficient in which B(z) is build from the FFT of the reflection
// coefficient and A(z)=1. See, e.g., Skaar et al, J Quantum Electron 37(2),
// 2001.
static inline INT
transfer_matrix_from_reflection_coefficient_TFMATRIX_CONTAINS_REFL_COEFF(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const nse_discretization_t discretization,
    const INT kappa)
{
    fft_wrapper_plan_t plan = fft_wrapper_safe_plan_init();
    INT ret_code = SUCCESS;
    UINT i;

    // Allocate some memory

    COMPLEX * const contspec_reordered = fft_wrapper_malloc(M * sizeof(COMPLEX));
    COMPLEX * const b_coeffs = fft_wrapper_malloc(M * sizeof(COMPLEX));
    if (contspec_reordered == NULL || b_coeffs == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    // Construct the polynomial B(z) that interpolates the given continuous
    // spectrum using a M-point FFT. Only the first deg+1 coefficients will
    // later be used to build the transfer matrix. Note that the order of
    // the computed coefficients is reversed.

    ret_code = remove_boundary_conds_and_reorder_for_fft(M, contspec, contspec_reordered,
                                              XI, D, T, discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_create_plan(&plan, M, contspec_reordered, b_coeffs,
                                       -1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_execute_plan(plan, contspec_reordered, b_coeffs);
    CHECK_RETCODE(ret_code, leave_fun);

    // Build the transfer matrix with B(z) computed above and A(z)=0.

    for (i=0; i<=deg; i++) {
        transfer_matrix[i] = 0.0;
        transfer_matrix[1*(deg+1) + i] = -kappa*CONJ(b_coeffs[M-1-deg+i]/M);
        transfer_matrix[2*(deg+1) + i] = b_coeffs[deg - i]/M;
        transfer_matrix[3*(deg+1) + i] = 0.0;
    }

    // Change A(z)=0 to A(z)=1.

    transfer_matrix[deg] = 1.0;
    transfer_matrix[3*(deg+1)] = 1.0;
    
leave_fun:
    fft_wrapper_destroy_plan(&plan);
    fft_wrapper_free(contspec_reordered);
    fft_wrapper_free(b_coeffs);
    return ret_code;
}

// This auxiliary function builds a transfer matrix from a given reflection
// coefficient using Algorithm 1 in http://arxiv.org/abs/1607.01305v2
// (defocusing case only).
static inline INT
transfer_matrix_from_reflection_coefficient_TFMATRIX_CONTAINS_AB_FROM_ITER(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const nse_discretization_t discretization,
    const INT kappa,
    const UINT max_iter)
{
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    fft_wrapper_plan_t plan_inv = fft_wrapper_safe_plan_init();
    INT ret_code = SUCCESS;
    UINT i, iter;

    if (D < 2 || (D&(D-1)) != 0)
        return E_INVALID_ARGUMENT(D);
    if (M != D)
        return E_INVALID_ARGUMENT(M);
    if (D != deg)
        return E_ASSERTION_FAILED;
    if (kappa != -1)
        return E_INVALID_ARGUMENT(kappa);

    // Allocate some memory and initial (inverse) FFT plans

    COMPLEX * const contspec_reordered = fft_wrapper_malloc(M * sizeof(COMPLEX));
    COMPLEX * const fft_in = fft_wrapper_malloc(M * sizeof(COMPLEX));
    COMPLEX * const fft_out = fft_wrapper_malloc(M * sizeof(COMPLEX));
    COMPLEX * const a_coeffs = fft_wrapper_malloc(D * sizeof(COMPLEX));
    COMPLEX * const b_coeffs = fft_wrapper_malloc(D * sizeof(COMPLEX));
    if (contspec_reordered == NULL || fft_in == NULL || fft_out == NULL
    || a_coeffs == NULL || b_coeffs == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    ret_code = fft_wrapper_create_plan(&plan_fwd, D, fft_in, b_coeffs, -1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_create_plan(&plan_inv, D, fft_in, fft_out, +1);
    CHECK_RETCODE(ret_code, leave_fun);

    // Remove phase factors due to initial conditions from contspec,
    // reorder if necessary

    ret_code = remove_boundary_conds_and_reorder_for_fft(M, contspec, contspec_reordered,
                                              XI, D, T, discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL prev_phase_change = FNFT_INF;
    REAL prev_phase_change_diff = FNFT_INF;
    for (iter=0; iter<max_iter; iter++) {

        // Construct the coefficients of B(z)

        for (i=0; i<D; i++) {
            const COMPLEX q = contspec_reordered[i];
            const REAL q_abs = CABS(q);
            fft_in[i] = q / SQRT( 1.0 + kappa*q_abs*q_abs ) / D;
        }

        ret_code = fft_wrapper_execute_plan(plan_fwd, fft_in, b_coeffs);
        CHECK_RETCODE(ret_code, leave_fun);
        for (i=0; i<D/2; i++) {
            const COMPLEX tmp = b_coeffs[i];
            b_coeffs[i] = b_coeffs[D - 1 - i];
            b_coeffs[D - 1 - i] = tmp;
        }

        // Construct the coefficients of A(z)

        ret_code = poly_specfact(D-1, b_coeffs, a_coeffs, 32, kappa);
        CHECK_RETCODE(ret_code, leave_fun);

        // Update the phases of the continuous spectrum

        for (i=0; i<D; i++)
            fft_in[i] = a_coeffs[D - 1 - i];
        ret_code = fft_wrapper_execute_plan(plan_inv, fft_in, fft_out);
        CHECK_RETCODE(ret_code, leave_fun);

        REAL cur_phase_change = 0.0;
        for (i=0; i<D; i++) {
            fft_out[i] = CARG( fft_out[i] );
            cur_phase_change += CABS( fft_out[i] ) / D;
        }
        for (i=0; i<=M/2; i++)
            contspec_reordered[i] = contspec[i+(M/2-1)]*CEXP(I*fft_out[i]);
        for (i=M/2+1; i<M; i++)
            contspec_reordered[i] = contspec[i-(M/2+1)]*CEXP(I*fft_out[i]);

        // Stop iterating if the phases stop changing or the correctios stop decreasing.

        const REAL cur_phase_change_diff =
            FABS( cur_phase_change - prev_phase_change );
        if (cur_phase_change_diff < 10*FNFT_EPSILON)
            break;
        prev_phase_change = cur_phase_change;
        if (cur_phase_change_diff > 0.9*prev_phase_change_diff)
            break;
        prev_phase_change_diff = cur_phase_change_diff;
    }
    if (iter == max_iter)
        WARN("Maximum number of iterations reached when constructing transfer matrix.")

    // Build the tranfer matrix

    for (i=0; i<D; i++) {
        transfer_matrix[1 + i] = a_coeffs[i];
        transfer_matrix[1*(deg+1) + i] = -kappa*CONJ(b_coeffs[D-1 - i]);
        transfer_matrix[2*(deg+1) + 1 + i] = b_coeffs[i];
        transfer_matrix[3*(deg+1) + i] = a_coeffs[D-1 - i];
    }
    transfer_matrix[0] = 0.0;
    transfer_matrix[2*(deg+1) - 1] = 0.0;
    transfer_matrix[2*(deg+1)] = 0.0;
    transfer_matrix[4*(deg+1) - 1] = 0.0;
        

leave_fun:
    fft_wrapper_free(contspec_reordered);
    fft_wrapper_free(fft_in);
    fft_wrapper_free(fft_out);
    fft_wrapper_free(a_coeffs);
    fft_wrapper_free(b_coeffs);
    fft_wrapper_destroy_plan(&plan_fwd);
    fft_wrapper_destroy_plan(&plan_inv);
    return ret_code;
}

// Auxiliary function. Calls the right method to build the transfer matrix
// for the inverse scattering stage if the continuous spectrum is provided
// in form of (samples of) a reflection coefficient.
static INT transfer_matrix_from_reflection_coefficient(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr)
{
    if (XI == NULL)
        return E_INVALID_ARGUMENT(XI);

    INT ret_code = SUCCESS;

    switch (opts_ptr->contspec_inversion_method) {

    case fnft_nsev_inverse_csmethod_DEFAULT:
        // fall through
    case fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_REFL_COEFF:

        ret_code = transfer_matrix_from_reflection_coefficient_TFMATRIX_CONTAINS_REFL_COEFF(
            M, contspec, XI, D, T, deg, transfer_matrix, opts_ptr->discretization, kappa);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    case fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER:

        ret_code = transfer_matrix_from_reflection_coefficient_TFMATRIX_CONTAINS_AB_FROM_ITER(
            M, contspec, XI, D, T, deg, transfer_matrix, opts_ptr->discretization, kappa,
            opts_ptr->max_iter);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    default:

        ret_code = E_INVALID_ARGUMENT(opts_ptr->contspec_inversion_method);
        goto leave_fun;

    }

leave_fun:
    return ret_code;
}

// Auxiliary function. Builds a transfer matrix from samples of the
// inverse Fourier transform B(tau) of b(xi). Uses the approach
// described in https://doi.org/10.1109/ECOC.2017.8346231
static INT transfer_matrix_from_B_of_tau(
    const UINT M,
    COMPLEX * const contspec,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr)
{
    if (M != D)
        return E_INVALID_ARGUMENT(M);
    if (T[0] != -T[1])
        return E_INVALID_ARGUMENT(T);
    if (opts_ptr->contspec_inversion_method
    != fnft_nsev_inverse_csmethod_DEFAULT)
        return E_INVALID_ARGUMENT(opts_ptr->contspec_inversion_method);

    INT ret_code = SUCCESS;
    UINT i;
    const REAL degree1step = nse_discretization_degree(opts_ptr->discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(discretization);
    COMPLEX * const b_coeffs = transfer_matrix + 2*(deg+1) + 1;
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    for (i=0; i<D; i++)
        b_coeffs[i] = 2*eps_t*contspec[i]/degree1step;
    b_coeffs[0] *= 0.5;
    b_coeffs[D-1] *= 0.5;

    COMPLEX * const a_coeffs = transfer_matrix + 1;
    ret_code = poly_specfact(D-1, b_coeffs, a_coeffs,
                             opts_ptr->oversampling_factor, kappa);
    CHECK_RETCODE(ret_code, leave_fun);

    for (i=0; i<D; i++) {
        transfer_matrix[1*(deg+1) + i] = -kappa*CONJ(b_coeffs[D-1 - i]);
        transfer_matrix[3*(deg+1) + i] = a_coeffs[D-1 - i];
    }
    transfer_matrix[0] = 0.0;
    transfer_matrix[2*(deg+1) - 1] = 0.0;
    transfer_matrix[2*(deg+1)] = 0.0;
    transfer_matrix[4*(deg+1) - 1] = 0.0;
    
leave_fun:
    return ret_code;
}
