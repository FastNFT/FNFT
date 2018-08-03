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
#include "fnft__misc.h"


static fnft_nsev_inverse_opts_t default_opts = {
    .discretization = nse_discretization_2SPLIT2A,
    .contspec_type = fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT,
    .contspec_inversion_method = fnft_nsev_inverse_csmethod_DEFAULT,
    .discspec_type = fnft_nsev_inverse_dstype_NORMING_CONSTANTS,
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

    INT ret_code = SUCCESS;
    COMPLEX XIC[2];
    // Note that z = exp(2.0*PI*I*i/M) = exp(map_coeff*I*xi*eps_t).
    // The range below results the same values as a normal M-point FFT,
    // but in a different order.
    const REAL eps_t = (T[1] - T[0]) / (D - 1);
    XIC[0] = CEXP(2.0*FNFT_PI*I * (M/2 + 1)/M);
    XIC[1] = CEXP(2.0*FNFT_PI*I * (M/2)/M) ;
    ret_code = nse_z_to_lambda(2, eps_t, &XIC[0], discretization);
    XI[0] = CREAL(XIC[0]);
    XI[1] = CREAL(XIC[1]);
    return ret_code;
}

// Auxliary functions. Function bodies follow below.
static INT transfer_matrix_from_reflection_coefficient(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT D,
    REAL const * const T,
    const UINT deg,
    COMPLEX * transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr);
static inline INT transfer_matrix_from_b_of_xi(
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
    COMPLEX * transfer_matrix,
    const INT kappa,
    fnft_nsev_inverse_opts_t const * const opts_ptr);
static INT add_discrete_spectrum(
    UINT const K,
    COMPLEX const * const bound_states,
    COMPLEX const * const normconsts_or_residues,
    const UINT D,
    COMPLEX * const q,
    REAL const * const T,
    const UINT contspec_flag,
    fnft_nsev_inverse_opts_t const * const opts_ptr);
static INT compute_eigenfunctions(
    UINT const K,
    COMPLEX const * const bound_states,
    const UINT D,
    COMPLEX * const q,
    REAL const * const T,
    COMPLEX * const phi,
    COMPLEX * const psi);
static INT precompensate_for_cdt_phaseshifts(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT K,
    COMPLEX const * const bound_states);

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
    if (M > 0 && contspec == NULL)
        return E_INVALID_ARGUMENT(contspec);
    if (contspec != NULL && M%2 != 0)
        return E_INVALID_ARGUMENT(M);
    if (contspec != NULL && M < D)
        return E_INVALID_ARGUMENT(M);
    if (D < 2 || (D&(D - 1)) != 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL || !(T[0] < T[1]))
        return E_INVALID_ARGUMENT(T);
    if (kappa != +1 && kappa != -1)
        return E_INVALID_ARGUMENT(kappa);
    if (K > 0 && kappa != +1)
        return E_SANITY_CHECK_FAILED(Discrete spectrum is present only in the focussing case(kappa=1).);
    if (K > 0 && bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    UINT i;
    for (i = 0; i < K; i++) {
        if (CIMAG(bound_states[i]) <= 0)
            return E_SANITY_CHECK_FAILED(bound_states should be stricly in the upper-half complex-plane.);
    }
    if (K > 0 && normconsts_or_residues == NULL)
        return E_INVALID_ARGUMENT(normconsts_or_residues);
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->discretization != nse_discretization_2SPLIT2A
            && opts_ptr->discretization != nse_discretization_2SPLIT2_MODAL)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);
    if (contspec == NULL && K == 0)
        return E_SANITY_CHECK_FAILED(Neither contspec nor discspec provided.);
    if (XI == NULL && contspec != NULL &&
        opts_ptr->contspec_type != fnft_nsev_inverse_cstype_B_OF_TAU)
        return E_INVALID_ARGUMENT(XI);

    INT ret_code = SUCCESS;
    INT contspec_flag = 0;
    COMPLEX *transfer_matrix = NULL;

    if (contspec != NULL) {

        // This flag is used by the add_discrete_spectrum to select method
        // contspec == NULL.
        contspec_flag = 1;

        // Allocate memory and initialize values that we will need later

        const UINT deg = D*nse_discretization_degree(opts_ptr->discretization);
        if (deg == 0)
            return E_INVALID_ARGUMENT(discretization);

        transfer_matrix = malloc(4*(deg+1) * sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }

        // Step 1: Construct the transfer matrix from the provided
        // representation of the continuous spectrum.

        switch (opts_ptr->contspec_type) {

        case fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT:
            ret_code = precompensate_for_cdt_phaseshifts(M, contspec, XI, K,
                                                          bound_states);
            CHECK_RETCODE(ret_code, leave_fun);
            ret_code = transfer_matrix_from_reflection_coefficient(M, contspec,
                       XI, D, T, deg, transfer_matrix, kappa, opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        case fnft_nsev_inverse_cstype_B_OF_XI:
            ret_code = transfer_matrix_from_b_of_xi(M, contspec, XI, D, T, deg,
                                                    transfer_matrix, kappa,
                                                    opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        case fnft_nsev_inverse_cstype_B_OF_TAU:
            ret_code = transfer_matrix_from_B_of_tau(M, contspec, D, T, deg,
                       transfer_matrix, kappa, opts_ptr);
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
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Availability of bound_states and normconsts_or_residues has
    // already been checked.
    if (K > 0) {
        ret_code = add_discrete_spectrum(K, bound_states,
                                         normconsts_or_residues,
                                         D, q, T, contspec_flag, opts_ptr);
        CHECK_RETCODE(ret_code, leave_fun);
    }

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
        fnft_nsev_inverse_opts_t const * const opts_ptr)
{
    UINT i;
    INT ret_code = SUCCESS;
    const REAL eps_t = (T[1] - T[0]) / (D - 1);
    const REAL eps_xi = (XI[1] - XI[0]) / (M - 1);
    REAL phase_factor = 0.0;

    switch (opts_ptr->contspec_type) {
    case fnft_nsev_inverse_cstype_REFLECTION_COEFFICIENT:
        ret_code = nse_phase_factor_rho(eps_t, T[1], &phase_factor,
                                        opts_ptr->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    case fnft_nsev_inverse_cstype_B_OF_XI:
        ret_code = nse_phase_factor_b(eps_t, D, T, &phase_factor,
                                      opts_ptr->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        break;

    default:
        ret_code = E_INVALID_ARGUMENT(opts->contspec_type);
        goto leave_fun;
    }

    for (i=0; i<M; i++) {
        const REAL xi = XI[0] + i*eps_xi;
        contspec[i] *= CEXP(-I*xi*phase_factor);
    }

    for (i=0; i<=M/2; i++)
        contspec_reordered[i] = contspec[i + (M/2 - 1)];
    for (i=M/2+1; i<M; i++)
        contspec_reordered[i] = contspec[i - (M/2 + 1)];

leave_fun:
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
        const INT kappa,
        fnft_nsev_inverse_opts_t const * const opts_ptr)
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

    ret_code = remove_boundary_conds_and_reorder_for_fft(M, contspec,
                                                         contspec_reordered,
                                                         XI, D, T, opts_ptr);
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
        const INT kappa,
        fnft_nsev_inverse_opts_t const * const opts_ptr)
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

    ret_code = remove_boundary_conds_and_reorder_for_fft(M, contspec,
                                                         contspec_reordered,
                                                         XI, D, T, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);

    REAL prev_phase_change = FNFT_INF;
    REAL prev_phase_change_diff = FNFT_INF;
    const UINT max_iter = opts_ptr->max_iter;
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

        // Build the transfer matrix

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
                M, contspec, XI, D, T, deg, transfer_matrix, kappa, opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        case fnft_nsev_inverse_csmethod_TFMATRIX_CONTAINS_AB_FROM_ITER:

            ret_code = transfer_matrix_from_reflection_coefficient_TFMATRIX_CONTAINS_AB_FROM_ITER(
                M, contspec, XI, D, T, deg, transfer_matrix, kappa, opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        default:

            ret_code = E_INVALID_ARGUMENT(opts_ptr->contspec_inversion_method);
            goto leave_fun;

    }

    leave_fun:
        return ret_code;
}

// Auxiliary function. Constructs a transfer matrix given b(xi). The transfer
// is build from A(z) and B(z), where first B(z) is formed from the given values
// of b(xi) using a FFT, and then A(z) is found using spectral factorization.
static inline INT transfer_matrix_from_b_of_xi(
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

    // Construct the polynomial B(z) that interpolates the given b(xi)
    // using a M-point FFT. Only the first deg+1 coefficients will
    // later be used to build the transfer matrix. Note that the order of
    // the computed coefficients is reversed.

    ret_code = remove_boundary_conds_and_reorder_for_fft(M, contspec,
                                                         contspec_reordered,
                                                         XI, D, T, opts_ptr);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_create_plan(&plan, M, contspec_reordered, b_coeffs,
                                       -1);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = fft_wrapper_execute_plan(plan, contspec_reordered, b_coeffs);
    CHECK_RETCODE(ret_code, leave_fun);

    // Add B(z) to the transfer matrix

    for (i=0; i<=deg; i++) {
        transfer_matrix[1*(deg+1) + i] = -kappa*CONJ(b_coeffs[M-1-deg+i]/M);
        transfer_matrix[2*(deg+1) + i] = b_coeffs[deg - i]/M;
    }

    // Compute A(z) and add it to the transfer matrix.

    ret_code = poly_specfact(deg, transfer_matrix+2*(deg+1), transfer_matrix,
                             opts_ptr->oversampling_factor, kappa);
    CHECK_RETCODE(ret_code, leave_fun);

    for (i=0; i<=deg; i++)
        transfer_matrix[3*(deg+1) + i] = transfer_matrix[deg-i];

    leave_fun:
        fft_wrapper_destroy_plan(&plan);
        fft_wrapper_free(contspec_reordered);
        fft_wrapper_free(b_coeffs);
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

static INT add_discrete_spectrum(
        UINT const K,
        COMPLEX const * const bound_states,
        COMPLEX const * const normconsts_or_residues,
        const UINT D,
        COMPLEX * const q,
        REAL const * const T,
        const UINT contspec_flag,
        fnft_nsev_inverse_opts_t const * const opts_ptr)
{
    if (K <= 0)
        return E_INVALID_ARGUMENT(K);
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (normconsts_or_residues == NULL)
        return E_INVALID_ARGUMENT(normconsts_or_residues);

    COMPLEX * bnd_states = NULL;
    COMPLEX * norm_consts = NULL;
    COMPLEX * bnd_states_conj = NULL;
    COMPLEX * bnd_states_diff = NULL;
    COMPLEX tmp;
    COMPLEX * phi = NULL;
    COMPLEX * psi = NULL;
    UINT i,j;
    UINT zc_point = 0; //index of zero-crossing of time
    INT ret_code = SUCCESS;
    REAL * t = NULL;
    REAL const eps_t = (T[1] - T[0])/(D-1);
    t = malloc(D * sizeof(REAL));
    bnd_states = malloc(K * sizeof(COMPLEX));
    norm_consts = malloc(K * sizeof(COMPLEX));
    bnd_states_conj = malloc(K * sizeof(COMPLEX));
    bnd_states_diff = malloc(K * sizeof(COMPLEX));
    if (t == NULL || bnd_states == NULL ||
            bnd_states_conj == NULL || bnd_states_diff == NULL ||
            norm_consts == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    //Building time vector
    // zero-crossing point is used to improve numerical conditioning
    j = 0; //Using j as flag to find zero-crossing
    for (i = 0; i < D; i++){
        t[i] = T[0] + eps_t*i;
        if (t[i] >= 0.0 && j == 0){
            zc_point = i;
            j = 1;//clear flag
        }
    }

    //Creating local copies of bound_states and normconsts_or_residues
    for (i = 0; i < K; i++){
        bnd_states[i] = bound_states[i];
        norm_consts[i] = normconsts_or_residues[i];
    }
    //sorting bnd_states in descending order based on magnitude of imaginary part
    for (i = 0; i < K; ++i){
        for (j = i + 1; j < K; ++j){
            if (CIMAG(bnd_states[i]) < CIMAG(bnd_states[j])){
                tmp =  bnd_states[i];
                bnd_states[i] = bnd_states[j];
                bnd_states[j] = tmp;
                tmp =  norm_consts[i];
                norm_consts[i] = norm_consts[j];
                norm_consts[j] = tmp;
            }
        }
    }
    //Checking for multiplicity of bound_states.
    for (i = 0; i < K-1; i++){
        if (bnd_states[i+1] == bnd_states[i]){
            ret_code = E_SANITY_CHECK_FAILED(Bound_states should be simple (multiplicity should be 1).);
            goto leave_fun;
        }
    }

    //Computing few values which will be used repeatdely in the code below

    for (i = 0; i < K; i++){
        bnd_states_conj[i] = CONJ(bnd_states[i]);
        bnd_states_diff[i] = 2*I*CIMAG(bnd_states[i]);
    }
    //CDT works with norming constants so residues need to be converted to
    // norming constants
    if (opts_ptr->discspec_type == fnft_nsev_inverse_dstype_RESIDUES){
        for (i = 0; i < K; i++){
            tmp = 1;
            for (j = 0; j < K; j++){
                if (j != i){
                    tmp = tmp*(bnd_states[i]-bnd_states[j])/(bnd_states[i]-CONJ(bnd_states[j]));
                }
            }
            norm_consts[i]=(norm_consts[i]/bnd_states_diff[i])*tmp;
        }

    }

    // In absence of contspec some tricks can be used to
    // speed up computation and improve numerical conditioning
    // of the multi-soliton solution
    if (contspec_flag == 0 &&
            opts_ptr->contspec_inversion_method != fnft_nsev_inverse_csmethod_USE_SEED_POTENTIAL_INSTEAD){

        COMPLEX qt, rho, rhok[K], rhoc, f;
        UINT n;
        for (n = zc_point; n < D; n++){
            qt = 0.0;
            for (i = 0; i < K; i++){
                rhok[i] = norm_consts[i]*CEXP(2*I*bnd_states[i]*t[n]);
            }
            for (i = 0; i < K; i++){
                rho = rhok[i];
                rhoc = CONJ(rho);
                f = bnd_states_diff[i]/(1+CABS(rho)*CABS(rho));
                qt = qt + 2*rhoc*f*I;
                for (j = i+1; j < K; j++){
                    rhok[j] = (((bnd_states[j]-bnd_states[i])*rhok[j]+
                            (rhok[j]-rho)*f)/(bnd_states[j]-bnd_states_conj[i]-(1+rhoc*rhok[j])*f));
                }
            }
            q[n] = qt;
        }

        for (i = 0; i < K; i++)
            norm_consts[i] = 1/norm_consts[i];

        for (n = 0; n < zc_point; n++){
            qt = 0.0;
            for (i = 0; i < K; i++){
                rhok[i] = norm_consts[i]*CEXP(-2*I*bnd_states[i]*t[n]);
            }
            for (i = 0; i < K; i++){
                rho = rhok[i];
                rhoc = CONJ(rho);
                f = bnd_states_diff[i]/(1+CABS(rho)*CABS(rho));
                qt = qt + 2*rhoc*f*I;
                for (j = i+1; j < K; j++){
                    rhok[j] = (((bnd_states[j]-bnd_states[i])*rhok[j]+
                            (rhok[j]-rho)*f)/(bnd_states[j]-bnd_states_conj[i]-(1+rhoc*rhok[j])*f));
                }
            }
            q[n] = CONJ(qt);
        }

    }
    else if ((contspec_flag == 0 &&
            opts_ptr->contspec_inversion_method == fnft_nsev_inverse_csmethod_USE_SEED_POTENTIAL_INSTEAD) ||
            (contspec_flag == 1 &&
            opts_ptr->contspec_inversion_method != fnft_nsev_inverse_csmethod_USE_SEED_POTENTIAL_INSTEAD)){

        // Allocate memory for computing eigenfuctions
        // In MATLAB notation
        //phi = zeros(2,D,K)
        //Here in C phi = [phi(1,1,1),phi(1,2,1),...,phi(1,D,1),phi(1,1,2),...
        // phi(1,D,K), phi(2,1,1),phi(2,2,1),...,phi(2,D,1),phi(2,1,2),...
        // phi(2,D,K)]
        phi = malloc((D*K*2) * sizeof(COMPLEX));
        psi = malloc((D*K*2) * sizeof(COMPLEX));
        if (phi == NULL || psi == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
//         // TODO: Find reason why this is needed
        INT sgn_fac = 1;
        if (K%2 != 0 &&
                opts_ptr->contspec_inversion_method == fnft_nsev_inverse_csmethod_USE_SEED_POTENTIAL_INSTEAD){
            sgn_fac = -1;
            for (i = 0; i < K; i++)
                norm_consts[i] = -norm_consts[i];
        }
//         //
        ret_code = compute_eigenfunctions(K, bnd_states, D, q, T, phi, psi);
        CHECK_RETCODE(ret_code, leave_fun);
        COMPLEX S1[K], S2[K], phi1, phi2, psi1, psi2, beta, qn;
        UINT n;
        for (n = 0; n < D; n++){
            qn = q[n];
            for (i = 0; i < K; i++){
                phi1 = phi[n + i*D];
                phi2 = phi[n + (K+i)*D];
                psi1 = psi[n + i*D];
                psi2 = psi[n + (K+i)*D];
                if (i > 0){
                    for (j = 0; j < i; j++){
                        tmp = (bnd_states[i]-S1[j])*phi1 -S2[j]*phi2;
                        phi2 = CONJ(S2[j])*phi1 + (bnd_states[i]-CONJ(S1[j]))*phi2;
                        phi1 = tmp;
                        tmp = (bnd_states[i]-S1[j])*psi1 -S2[j]*psi2;
                        psi2 = CONJ(S2[j])*psi1 + (bnd_states[i]-CONJ(S1[j]))*psi2;
                        psi1 = tmp;
                    }
                }
                beta = (phi1 - norm_consts[i]*psi1)/(phi2 - norm_consts[i]*psi2);
                tmp = CABS(beta)*CABS(beta);
                S1[i] = (tmp*bnd_states[i] + CONJ(bnd_states[i]))/(1 + tmp);
                S2[i] = (2*I*CIMAG(bnd_states[i])*beta)/(1 + tmp);
                qn = qn - 2*I*S2[i];

            }
            q[n] = sgn_fac*qn;
        }
    }
    else
        return E_INVALID_ARGUMENT(opts_ptr->contspec_inversion_method);

    leave_fun:
        free(t);
        free(phi);
        free(psi);
        free(bnd_states);
        free(bnd_states_conj);
        free(bnd_states_diff);
        free(norm_consts);
        return ret_code;
}

static INT compute_eigenfunctions(
        UINT const K,
        COMPLEX const * const bound_states,
        const UINT D,
        COMPLEX * const q,
        REAL const * const T,
        COMPLEX * phi,
        COMPLEX * psi)
{
    INT ret_code = SUCCESS;
    UINT i, n;
    COMPLEX * phi1_ptr, * phi2_ptr;
    COMPLEX * psi1_ptr, * psi2_ptr;
    COMPLEX l, ks, k, scl, ch, sh, u1, tmp1, tmp2;
    const REAL eps_t = ((T[1] - T[0])/(D - 1))/2; //Note this is half of the
    // time-step defined elsewhere
    // computing-phi
    for (i = 0; i < K; i++) { // iterate over bound_states
        l = bound_states[i];
        // Assign local pointers
        phi1_ptr = phi + i*D;
        phi2_ptr = phi1_ptr + K*D;
        //Apply boundary condition
        phi1_ptr[0] = CEXP(-I*l*T[0]);
        phi2_ptr[0] = 0.0;

        ks = (-(CABS(q[0])*CABS(q[0]))-(l*l));
        k = CSQRT(ks);
        ch = CCOSH(k*eps_t);
        sh = CSINH(k*eps_t)/k;
        u1 = l*sh*I;
        for (n = 1; n < D; n++){
            //First half-step
            if (ks != 0){
                phi1_ptr[n] = (ch-u1)*phi1_ptr[n-1] + (q[n-1]*sh)*phi2_ptr[n-1];
                phi2_ptr[n] = (-CONJ(q[n-1])*sh)*phi1_ptr[n-1] + (ch + u1)*phi2_ptr[n-1];
            }
            else{
                phi1_ptr[n] = phi1_ptr[n-1];
                phi2_ptr[n] = phi2_ptr[n-1];
            }
            //Second half-step
            ks = (-(CABS(q[n])*CABS(q[n]))-(l*l));
            k = CSQRT(ks);
            ch = CCOSH(k*eps_t);
            sh = CSINH(k*eps_t)/k;
            u1 = l*sh*I;
            if (ks != 0){
                tmp1 = (ch-u1)*phi1_ptr[n] + (q[n]*sh)*phi2_ptr[n];
                tmp2 = (-CONJ(q[n])*sh)*phi1_ptr[n] + (ch + u1)*phi2_ptr[n];
                phi1_ptr[n] = tmp1;
                phi2_ptr[n] = tmp2;
            }
        }
    }
    // computing-psi
    for (i = 0; i < K; i++) { // iterate over bound_states
        l = bound_states[i];
        // Assign local pointers
        psi1_ptr = psi + i*D;
        psi2_ptr = psi1_ptr + K*D;
        //Apply boundary condition
        psi1_ptr[D-1] = 0.0;
        psi2_ptr[D-1] = CEXP(I*l*T[1]);

        ks = (-(CABS(q[D-1])*CABS(q[D-1]))-(l*l));
        k = CSQRT(ks);
        ch = CCOSH(k*eps_t);
        sh = CSINH(k*eps_t)/k;
        u1 = l*sh*I;
        for (n = D-1; n > 0; n--){
            //First half-step
            if (ks != 0){
                scl = ((ch-u1)*(ch + u1))-(( -CONJ(q[n])*sh)*(q[n]*sh));
                psi1_ptr[n-1] = ((ch + u1)*psi1_ptr[n] + (-q[n]*sh)*psi2_ptr[n])/scl;
                psi2_ptr[n-1] = ((CONJ(q[n])*sh)*psi1_ptr[n] + (ch - u1)*psi2_ptr[n])/scl;
            }
            else{
                psi1_ptr[n-1] = psi1_ptr[n];
                psi2_ptr[n-1] = psi2_ptr[n];
            }
            //Second half-step
            ks = (-(CABS(q[n-1])*CABS(q[n-1]))-(l*l));
            k = CSQRT(ks);
            ch = CCOSH(k*eps_t);
            sh = CSINH(k*eps_t)/k;
            u1 = l*sh*I;
            if (ks != 0){
                scl = ((ch-u1)*(ch + u1))-(( -CONJ(q[n-1])*sh)*(q[n-1]*sh));
                tmp1 = ((ch + u1)*psi1_ptr[n-1] + (-q[n-1]*sh)*psi2_ptr[n-1])/scl;
                tmp2 = ((CONJ(q[n-1])*sh)*psi1_ptr[n-1] + (ch - u1)*psi2_ptr[n-1])/scl;
                psi1_ptr[n-1] = tmp1;
                psi2_ptr[n-1] = tmp2;
            }
        }
    }
    return ret_code;

}

// Auxiliary function. The Darboux transform has the side-effect of multiplying
// the a(xi) of the seed potentials with Blaschke factors when adding solitons
// on top of a given seed potential. This routine precompensates for this effect
// by multiplying these Blaschke factors to the reflection coefficient that
// is used to generate the seed for the Darboux transform in the full INFT.
static INT precompensate_for_cdt_phaseshifts(
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    const UINT K,
    COMPLEX const * const bound_states)
{
    UINT i, k;

    if (K == 0) // nothing to do
        return SUCCESS;

    const REAL eps_xi = (XI[1] - XI[0])/(M - 1);
    for (i=0; i<M; i++) {
        const REAL xi = XI[0] + i*eps_xi;
        for (k=0; k<K; k++)
            contspec[i] *= (xi - bound_states[k])/(xi - CONJ(bound_states[k]));
    }

    return SUCCESS;
}
