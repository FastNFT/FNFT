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
 * Shrinivas Chimmalgi (TU Delft) 2017-2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include <string.h> // for memcpy
#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_chirpz.h"
#include "fnft_nsev.h"
#include "fnft__nse_fscatter.h"
#include "fnft__nse_scatter.h"
#include "fnft__nse_discretization.h"
#include "fnft__misc.h" // for l2norm

static fnft_nsev_opts_t default_opts = {
    .bound_state_filtering = nsev_bsfilt_FULL,
    .bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE,
    .niter = 10,
    .Dsub = 0, // auto
    .discspec_type = nsev_dstype_NORMING_CONSTANTS,
    .contspec_type = nsev_cstype_REFLECTION_COEFFICIENT,
    .normalization_flag = 1,
    .discretization = nse_discretization_2SPLIT4B
};

/**
 * Creates a new options variable for fnft_nsev with default settings.
 * See the header file for a detailed description.
 */
fnft_nsev_opts_t fnft_nsev_default_opts()
{
    return default_opts;
}

/**
 * Returns the maximum number of bound states that can be detected by
 * fnft_nsev. See header file for details.
 */
UINT fnft_nsev_max_K(const UINT D, fnft_nsev_opts_t const * const opts)
{
    if (opts != NULL)
        return nse_discretization_degree(opts->discretization) * D;
    else
        return nse_discretization_degree(default_opts.discretization) * D;
}

/**
 * Declare auxiliary routines used by the main routine fnft_nsev.
 * Their bodies follow below.
 */
static inline INT tf2contspec(
    const UINT deg,
    const INT W,
    COMPLEX * const transfer_matrix,
    REAL const * const T,
    const UINT D,
    REAL const * const XI,
    const UINT M,
    COMPLEX *result,
    fnft_nsev_opts_t * const opts);

static inline INT tf2boundstates(
    UINT D,
    COMPLEX const * const q,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    REAL const * const T,
    const REAL eps_t,
    UINT * const K_ptr,
    COMPLEX * const bound_states,
    fnft_nsev_opts_t * const opts);

static inline INT tf2normconsts_or_residues(
    const UINT D,
    COMPLEX const * const q,
    REAL const * const T,
    const UINT K,
    COMPLEX * const transfer_matrix,
    const UINT deg,
    COMPLEX * const bound_states,
    COMPLEX * const normconsts_or_residues,
    fnft_nsev_opts_t * const opts);

static inline INT refine_roots_newton(
    const UINT D,
    COMPLEX const * const q,
    REAL const * const T,
    const UINT K,
    COMPLEX * bound_states,
    nse_discretization_t discretization,
    const UINT niter);

/**
 * Fast nonlinear Fourier transform for the nonlinear Schroedinger
 * equation with vanishing boundary conditions.
 * See the header file for documentation.
 */
INT fnft_nsev(
    const UINT D,
    COMPLEX * const q,
    REAL const * const T,
    const UINT M,
    COMPLEX * const contspec,
    REAL const * const XI,
    UINT * const K_ptr,
    COMPLEX * const bound_states,
    COMPLEX * const normconsts_or_residues,
    const INT kappa,
    fnft_nsev_opts_t *opts)
{
    COMPLEX *transfer_matrix = NULL;
    COMPLEX *qsub = NULL;
    REAL eps_t;
    UINT deg;
    INT W = 0, *W_ptr = NULL;
    INT ret_code = SUCCESS;
    UINT i;
    
    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL || T[0] >= T[1])
        return E_INVALID_ARGUMENT(T);
    if (contspec != NULL) {
        if (XI == NULL || XI[0] >= XI[1])
            return E_INVALID_ARGUMENT(XI);
    }
    if (abs(kappa) != 1)
        return E_INVALID_ARGUMENT(kappa);
    if (bound_states != NULL) {
        if (K_ptr == NULL)
            return E_INVALID_ARGUMENT(K_ptr);
    }
    if (opts == NULL)
        opts = &default_opts;
    
    // Allocate memory for the transfer matrix. Note that after computation
    // of the transfer matrix, the second and fourth quarter of the
    // array carry redundant information that is not used. These quarters
    // are therefore used as buffers and may be overwritten at some point.
    i = nse_fscatter_numel(D, opts->discretization);
    if (i == 0) { // size D>=2, this means unknown discretization
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto release_mem;
    }
    transfer_matrix = malloc(i*sizeof(COMPLEX));
    if (transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }
    
    // Determine step size
    eps_t = (T[1] - T[0])/(D - 1);

    // Compute the transfer matrix
    if (opts->normalization_flag)
        W_ptr = &W;
    ret_code = nse_fscatter(D, q, eps_t, kappa, transfer_matrix, &deg, W_ptr,
        opts->discretization);
    CHECK_RETCODE(ret_code, release_mem);
    
    // Compute the continuous spectrum
    if (contspec != NULL && M > 0) {
        ret_code = tf2contspec(deg, W, transfer_matrix, T, D, XI, M,
            contspec, opts);
        CHECK_RETCODE(ret_code, release_mem);
    }
    
    // Compute the discrete spectrum
    if (kappa == +1 && bound_states != NULL) {
        
        // Compute the bound states
        if (opts->bound_state_localization == nsev_bsloc_SUBSAMPLE_AND_REFINE) {
            // the mixed method gets special treatment
            
            // First step: Find initial guesses for the bound states using the
            // fast eigenvalue method. To bound the complexity, a subsampled
            // version of q, qsub, will be passed to the fast eigenroutine.
            UINT Dsub = opts->Dsub;
            if (Dsub == 0) // The user wants us to determine Dsub
                Dsub = SQRT(D * LOG2(D) * LOG2(D));
            UINT first_last_index[2];
            ret_code = misc_downsample(D, q, &Dsub, &qsub, first_last_index);
            CHECK_RETCODE(ret_code, release_mem);
            REAL const Tsub[2] = { T[0] + first_last_index[0]*eps_t,
                T[0] + first_last_index[1]*eps_t };
          
            // Fixed bound states of qsub using the fast eigenvalue method
            opts->bound_state_localization = nsev_bsloc_FAST_EIGENVALUE;
            ret_code = fnft_nsev(Dsub, qsub, Tsub, 0, NULL, XI, K_ptr,
                bound_states, NULL, kappa, opts);
            CHECK_RETCODE(ret_code, release_mem);
           
            // Second step: Refine the found bound states using Newton's method
            // on the full signal.
            opts->bound_state_localization = nsev_bsloc_NEWTON;
            ret_code = tf2boundstates(D, q, deg, transfer_matrix, T,
                    eps_t, K_ptr, bound_states, opts);
            CHECK_RETCODE(ret_code, release_mem);
           
            // Restore original state of opts
            opts->bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE;
            
        } else { // any other method is handled directly by the subroutine
            
            ret_code = tf2boundstates(D, q, deg, transfer_matrix, T,
                    eps_t, K_ptr, bound_states, opts);
            CHECK_RETCODE(ret_code, release_mem);

        }
        // Norming constants and/or residues)
        if (normconsts_or_residues != NULL && *K_ptr != 0) {

            ret_code = tf2normconsts_or_residues(D, q, T, *K_ptr,
                transfer_matrix, deg, bound_states, normconsts_or_residues,
                opts);
            CHECK_RETCODE(ret_code, release_mem);

        }
    } else if (K_ptr != NULL) {
        *K_ptr = 0;
    }
    
release_mem:
    free(transfer_matrix);
    free(qsub);
        
    return ret_code;
}

// Auxiliary function: Computes continuous spectrum on a frequency grid
// from a given transfer matrix.
static inline INT tf2contspec(
    const UINT deg,
    const INT W,
    COMPLEX * const transfer_matrix,
    REAL const * const T,
    const UINT D,
    REAL const * const XI,
    const UINT M,
    COMPLEX * const result,
    fnft_nsev_opts_t * const opts)
{
    COMPLEX *a_vals;
    COMPLEX A, V;
    REAL eps_t, eps_xi;
    REAL xi, map_coeff, bnd_coeff, scale;
    REAL phase_factor_rho, phase_factor_a, phase_factor_b;
    INT ret_code;
    UINT i, offset = 0;
    
    // We reuse an unused part of the transfer matrix as buffer, if possible
    if (M <= deg+1) {
        a_vals = transfer_matrix + (deg+1);
    } else {
        a_vals = malloc(M*sizeof(COMPLEX));
        if (a_vals == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
    }
    
    // Set step sizes
    eps_t = (T[1] - T[0])/(D - 1);
    eps_xi = (XI[1] - XI[0])/(M - 1);

    // Retrieve some discretization-specific constants
    map_coeff = nse_discretization_mapping_coeff(opts->discretization);
    bnd_coeff = nse_discretization_boundary_coeff(opts->discretization);
    if (map_coeff == NAN || bnd_coeff == NAN) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }

    // Prepare the use of the chirp transform. The entries of the transfer
    // matrix that correspond to a and b will be evaluated on the frequency
    // grid xi(i) = XI1 + i*eps_xi, where i=0,...,M-1. Since
    // z=exp(map_coeff*j*xi*eps_t), we find that the z at which z the transfer
    // matrix has to be evaluated are given by z(i) = 1/(A * V^-i), where:
    V = CEXP(map_coeff*I*eps_xi*eps_t);
    A = CEXP(-map_coeff*I*XI[0]*eps_t);

    // Compute the continuous spectrum
    switch (opts->contspec_type) {

    case nsev_cstype_BOTH:

        offset = M;
        // fall through
    case nsev_cstype_REFLECTION_COEFFICIENT:

        ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, a_vals);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix+2*(deg+1), A, V, M,result);
        CHECK_RETCODE(ret_code, leave_fun);

        phase_factor_rho = -2.0*(T[1] + eps_t*bnd_coeff);
        for (i = 0; i < M; i++) {
            xi = XI[0] + i*eps_xi;
            if (a_vals[i] == 0.0) {
                ret_code = E_DIV_BY_ZERO;
                goto leave_fun;
            }
            result[i] *= CEXP(I*xi*phase_factor_rho) / a_vals[i];
        }

        if (opts->contspec_type == nsev_cstype_REFLECTION_COEFFICIENT)
            break;
        // fall through
    case nsev_cstype_AB:

        ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, result + offset);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix+2*(deg+1), A, V, M,
            result + offset + M);
        CHECK_RETCODE(ret_code, leave_fun);

        scale = POW(2.0, W); // needed since the transfer matrix might
                                  // have been scaled by nse_fscatter

        if (bnd_coeff == 0.0) {
            phase_factor_a = T[0] + T[1];
            phase_factor_b = T[0] - T[1];
        } else {
            phase_factor_a =  2.0*D*eps_t + T[1] + T[0];
            phase_factor_b = 2.0*(D-bnd_coeff)*eps_t + T[0] - T[1];
        }

        for (i = 0; i < M; i++) {
            xi = XI[0] + i*eps_xi;       
            result[offset + i] *= scale * CEXP(I*xi*phase_factor_a);
        	result[offset + M + i] *= scale * CEXP(I*xi*phase_factor_b);
        }

        break;

    default:

        ret_code = E_INVALID_ARGUMENT(opts->contspec_type);
        goto leave_fun;
    }
    
leave_fun:
    if (a_vals != transfer_matrix + (deg+1))
        free(a_vals);
    return ret_code;
}

// Auxiliary function for filtering: We assume that bound states must have
// real part in the interval [-re_bound, re_bound].
static inline REAL re_bound(REAL eps_t, REAL map_coeff)
{
    // At least for discretizations in which the continuous-time
    // spectral parameter lam is mapped to z=exp(map_coeff*j*lam*eps_t), we
    // can only resolve the region
    // -pi/(map_coeff*eps_t)<Re(lam)<pi/(map_coeff*eps_t).
    // Numerical artefacts often occur close to the border of this
    // region, which is why we filter such bound_states
    return 0.9*PI/FABS(map_coeff * eps_t);
}

// Auxiliary function for filtering: We assume that bound states must have
// imaginary part in the interval [0, im_bound].
static inline REAL im_bound(const UINT D, COMPLEX const * const q,
    REAL const * const T)
{
    // The nonlinear Parseval relation tells us that the squared L2 norm of
    // q(t) is >= 4*(sum of the imaginary parts of the bound states). Thus,
    // any bound state with an imaginary part greater than four times the
    // squared L2 norm of q(t) can be removed. A factor of 1.5 has been
    // added to account for numerical discrepancies when computing the norm
    // numerically (e.g., truncation errors or large step sizes).
    return 1.5 * 0.25 * misc_l2norm2(D, q, T[0], T[1]);
}


// Auxiliary function: Computes the bound states from a given transfer matrix.
static inline INT tf2boundstates(
    const UINT D,
    COMPLEX const * const q,
    const UINT deg,
    COMPLEX * const transfer_matrix,
    REAL const * const T,
    const REAL eps_t,
    UINT * const K_ptr,
    COMPLEX * const bound_states,
    fnft_nsev_opts_t * const opts)
{
    REAL map_coeff;
    UINT i, K;
    REAL bounding_box[4] = { NAN };
    COMPLEX * buffer = NULL;
    INT ret_code = SUCCESS;

    map_coeff = nse_discretization_mapping_coeff(opts->discretization);
    if (map_coeff == NAN)
        return E_INVALID_ARGUMENT(opts->discretization);

    // Localize bound states ...
    switch (opts->bound_state_localization) {
        
        // ... using Newton's method
        case nsev_bsloc_NEWTON:

            K = *K_ptr;
            buffer = bound_states; // Store intermediate results directly

            // Perform Newton iterations. Initial guesses of bound-states
            // should be in the continuous-time domain.
            ret_code = refine_roots_newton(D, q, T, K, buffer,
                nse_discretization_BO, opts->niter);
            CHECK_RETCODE(ret_code, leave_fun);
            
            break;
            
        // ... using the fast eigenvaluebased root finding
        case nsev_bsloc_FAST_EIGENVALUE:
            
            K = deg;
            if (*K_ptr >= K) {
                buffer = bound_states;
            } else {
                // Store intermediate results in unused part of transfer matrix.
                // This buffer is large enough to store all deg roots of the
                // polynomial, while bound_states provided by the user might be
                // smaller. The latter only needs to store the bound states that
                // survive the filtering.
                buffer = transfer_matrix + (deg+1);
            }

            ret_code = poly_roots_fasteigen(deg, transfer_matrix, buffer);
            CHECK_RETCODE(ret_code, leave_fun);

            // Roots are returned in discrete-time domain -> coordinate
            // transform (from discrete-time to continuous-time domain).
            for (i = 0; i < K; i++)
                buffer[i] = CLOG(buffer[i]) / (map_coeff*I*eps_t);
            
            break;
            
        default:
            
            return E_INVALID_ARGUMENT(opts->bound_state_localization);
    }
    
    // Filter bound states
    if (opts->bound_state_filtering != nsev_bsfilt_NONE) {

        bounding_box[0] = -INFINITY;
        bounding_box[1] = INFINITY;
        bounding_box[2] = 0.0;
        bounding_box[3] = INFINITY;
        ret_code = misc_filter(&K, buffer, NULL, bounding_box);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = misc_merge(&K, buffer, SQRT(EPSILON));
        CHECK_RETCODE(ret_code, leave_fun);

    }
    if (opts->bound_state_filtering == nsev_bsfilt_FULL) {

        bounding_box[1] = re_bound(eps_t, map_coeff);
        bounding_box[0] = -bounding_box[1];
        bounding_box[3] = im_bound(D, q, T);
        bounding_box[2] = 0;
        ret_code = misc_filter(&K, buffer, NULL, bounding_box);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Copy result from buffer to user-supplied array (if not identical)
    if (buffer != bound_states) {
        if (*K_ptr < K) {
            WARN("Found more than *K_ptr bound states. Returning as many as possible.");
            K = *K_ptr;
        }
        memcpy(bound_states, buffer, K * sizeof(COMPLEX));
    }
    
    // Update number of bound states
    *K_ptr = K;

leave_fun:    
    return ret_code;
}

// Auxiliary function: Computes the norming constants and/or residues
// using the BO scheme
static inline INT tf2normconsts_or_residues(
    const UINT D,
    COMPLEX const * const q,
    REAL const * const T,
    const UINT K,
    COMPLEX * const transfer_matrix,
    const UINT deg,
    COMPLEX * const bound_states,
    COMPLEX * const normconsts_or_residues,
    fnft_nsev_opts_t * const opts)
{
    
    COMPLEX *a_vals = NULL, *aprime_vals = NULL;
    UINT trunc_index;
    UINT i, offset = 0;
    INT ret_code = SUCCESS;
   
    // Check inputs
    if (K == 0) // no bound states to refine
        return SUCCESS;
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    
    // Reuse no longer used parts of the transfer matrix as buffers
    a_vals = transfer_matrix + (deg+1);
    aprime_vals = transfer_matrix + 3*(deg+1);
    
    // trunc_index is the index where we will split
    // trunc_index should be integer between 0 and D-1
    // trunc_index = D corresponds to splitting based on L1-norm
    trunc_index = D;
    ret_code = nse_scatter_bound_states(D, q, T, &trunc_index, K,
        bound_states, a_vals, aprime_vals, normconsts_or_residues, nse_discretization_BO);
    CHECK_RETCODE(ret_code, leave_fun);    

    // Update to or add residues if requested
    if (opts->discspec_type != nsev_dstype_NORMING_CONSTANTS) {
        
        if (opts->discspec_type == nsev_dstype_RESIDUES) {
            offset = 0;
        } else if (opts->discspec_type == nsev_dstype_BOTH) {
            offset = K;
            memcpy(normconsts_or_residues + offset,
                    normconsts_or_residues,
                    offset*sizeof(complex double));
        } else
            return E_INVALID_ARGUMENT(opts->discspec_type);
        
        // Divide norming constants by derivatives to get residues
        for (i = 0; i < K; i++) {
            if (aprime_vals[i] == 0.0)
                return E_DIV_BY_ZERO;
            normconsts_or_residues[offset + i] /= aprime_vals[i];
        }
    }

leave_fun:
    return ret_code;
}

// Auxiliary function: Refines the bound-states using Newtons method
static inline INT refine_roots_newton(
    const UINT D,
    COMPLEX const * const q,
    REAL const * const T,
    const UINT K,
    COMPLEX * bound_states,
    nse_discretization_t discretization,
    const UINT niter)     
{
    INT ret_code = SUCCESS;
    UINT i, iter;
    COMPLEX a_val, b_val, aprime_val, error;
    REAL eprecision = EPSILON * 100;
    REAL re_bound_val, im_bound_val, eps_t, map_coeff;
    UINT trunc_index;
    trunc_index = D;
    eps_t = (T[1] - T[0])/(D - 1);
    
    // Check inputs
    if (K == 0) // no bound states to refine
        return SUCCESS;
    if (niter == 0) // no refinement requested
        return SUCCESS;
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    
    im_bound_val = im_bound(D, q, T);
    if (im_bound_val == NAN)
        return E_OTHER("Upper bound on imaginary part of bound states is NaN");
    
    map_coeff = nse_discretization_mapping_coeff(discretization);
    if (map_coeff == NAN)
        return E_INVALID_ARGUMENT(discretization);
    re_bound_val = re_bound(eps_t, discretization);
        
    // Perform iterations of Newton's method
    for (i = 0; i < K; i++) {
        iter = 0;
        do {
            // Compute a(lam) and a'(lam) at the current root
            ret_code = nse_scatter_bound_states(D, q, T, &trunc_index, 1,
                bound_states + i, &a_val, &aprime_val, &b_val, discretization);
            if (ret_code != SUCCESS)
                return E_SUBROUTINE(ret_code);

            // Perform Newton updates: lam[i] <- lam[i] - a(lam[i])/a'(lam[i])
            if (aprime_val == 0.0)
                return E_DIV_BY_ZERO;
            error = a_val / aprime_val;
            bound_states[i] -= error;
            iter++;

            if (CIMAG(bound_states[i]) > im_bound_val
                || CREAL(bound_states[i]) > re_bound_val
                || CREAL(bound_states[i]) < -re_bound_val
                || CIMAG(bound_states[i]) < 0.0)
            break;

        } while (CABS(error) > eprecision && iter < niter);       
    }
    
    return SUCCESS;
}
