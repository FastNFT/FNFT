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
 * Marius Brehler (TU Dortmund) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES


#include "fnft_nsev.h"



static fnft_nsev_opts_t default_opts = {
    .bound_state_filtering = nsev_bsfilt_FULL,
    .bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE,
    .niter = 10,
    .Dsub = 0, // auto
    .discspec_type = nsev_dstype_NORMING_CONSTANTS,
    .contspec_type = nsev_cstype_REFLECTION_COEFFICIENT,
    .normalization_flag = 1,
    .discretization = nse_discretization_2SPLIT4B,
    .richardson_extrapolation_flag = 0
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

static inline INT nsev_compute_boundstates(
        UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        const UINT deg,
        COMPLEX * const transfer_matrix,
        REAL const * const T,
        const REAL eps_t,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        fnft_nsev_opts_t * const opts);

static inline INT fnft_nsev_base(
        const UINT D,
        COMPLEX * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT M,
        COMPLEX * const contspec,
        REAL const * const XI,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        const INT kappa,
        fnft_nsev_opts_t *opts);

static inline INT nsev_compute_contspec(
        const UINT deg,
        const INT W,
        COMPLEX * const transfer_matrix,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT D,
        REAL const * const XI,
        const UINT M,
        COMPLEX * const result,
        const INT kappa,
        fnft_nsev_opts_t * const opts);

static inline INT nsev_compute_normconsts_or_residues(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT K,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        fnft_nsev_opts_t * const opts);

static inline INT nsev_refine_bound_states_newton(const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        UINT K,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        FNFT_UINT niter);

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
    COMPLEX *qsub_preprocessed = NULL;
    COMPLEX *rsub_preprocessed = NULL;
    COMPLEX *q_preprocessed = NULL;
    COMPLEX *r_preprocessed = NULL;
    UINT Dsub = 0;
    REAL Tsub[2] = {0.0 ,0.0};
    UINT first_last_index[2] = {0};
    UINT K_sub;
    COMPLEX *contspec_sub = NULL;
    COMPLEX *bound_states_sub = NULL;
    COMPLEX *normconsts_or_residues_sub = NULL;
    COMPLEX *normconsts_or_residues_reserve = NULL;
    INT bs_loc_opt = 0, ds_type_opt = 0;
    INT ret_code = SUCCESS;
    UINT i, j, upsampling_factor, D_effective, nskip_per_step;
    
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
    
    // This switch checks for incompatible bound_state_localization options
    switch (opts->discretization) {
        case nse_discretization_2SPLIT2_MODAL:
        case nse_discretization_2SPLIT1A:
        case nse_discretization_2SPLIT1B:
        case nse_discretization_2SPLIT2A:
        case nse_discretization_2SPLIT2B:
        case nse_discretization_2SPLIT2S:
        case nse_discretization_2SPLIT3S:
        case nse_discretization_2SPLIT4B:
        case nse_discretization_2SPLIT3A:
        case nse_discretization_2SPLIT3B:
        case nse_discretization_2SPLIT4A:
        case nse_discretization_2SPLIT6B:
        case nse_discretization_2SPLIT6A:
        case nse_discretization_2SPLIT8B:
        case nse_discretization_2SPLIT5A:
        case nse_discretization_2SPLIT5B:
        case nse_discretization_2SPLIT8A:
        case nse_discretization_2SPLIT7A:
        case nse_discretization_2SPLIT7B:
        case nse_discretization_4SPLIT4A:
        case nse_discretization_4SPLIT4B:
            break;
        case nse_discretization_BO:
        case nse_discretization_CF4_2:
        case nse_discretization_CF4_3:
        case nse_discretization_CF5_3:
        case nse_discretization_CF6_4:
        case nse_discretization_ES4:
        case nse_discretization_TES4:
            if (opts->bound_state_localization != nsev_bsloc_NEWTON && kappa == +1){
                ret_code = E_INVALID_ARGUMENT(opts->bound_state_localization);
                goto leave_fun;
            }
            break;
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(opts->discretization);
    }
    
    
    upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    D_effective = D * upsampling_factor; // upsampling_factor*D gives the effective number of samples
    
    // Determine step size
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    
    normconsts_or_residues_reserve = normconsts_or_residues;
    if (opts->richardson_extrapolation_flag == 1){
        ds_type_opt = opts->discspec_type;
        if (ds_type_opt == nsev_dstype_RESIDUES){
            opts->discspec_type = nsev_dstype_BOTH;
            normconsts_or_residues_reserve = malloc(*K_ptr*2 * sizeof(COMPLEX));
            if (normconsts_or_residues_reserve == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
        }
    }
    
    // Preprocess signal
    Dsub = D;
    ret_code = nse_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
            first_last_index, opts->discretization);
    CHECK_RETCODE(ret_code, leave_fun);
    
    
    if (kappa == +1 && bound_states != NULL && opts->bound_state_localization == nsev_bsloc_SUBSAMPLE_AND_REFINE) {
        // the mixed method gets special treatment
        
        // First step: Find initial guesses for the bound states using the
        // fast eigenvalue method. To bound the complexity, a subsampled
        // version of q, qsub, will be passed to the fast eigenroutine.
        Dsub = opts->Dsub;
        if (Dsub == 0) // The user wants us to determine Dsub
            Dsub = SQRT(D * LOG2(D) * LOG2(D));
        nskip_per_step = ROUND((REAL)D / Dsub);
        Dsub = ROUND((REAL)D / nskip_per_step); // actual Dsub
        
        ret_code = nse_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
                first_last_index, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        
        Tsub[0] = T[0] + first_last_index[0] * eps_t;
        Tsub[1] = T[0] + first_last_index[1] * eps_t;
        
        // Fixed bound states of qsub using the fast eigenvalue method
        opts->bound_state_localization = nsev_bsloc_FAST_EIGENVALUE;
        ret_code = fnft_nsev_base(Dsub * upsampling_factor, qsub_preprocessed, rsub_preprocessed, Tsub, 0, NULL, XI, K_ptr,
                bound_states, NULL, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
        
        // Second step: Refine the found bound states using Newton's method
        // on the full signal and compute continuous spectrum
        opts->bound_state_localization = nsev_bsloc_NEWTON;
        ret_code = fnft_nsev_base(D_effective, q_preprocessed, r_preprocessed, T, M, contspec, XI, K_ptr,
                bound_states, normconsts_or_residues_reserve, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
        
        // Restore original state of opts
        opts->bound_state_localization = nsev_bsloc_SUBSAMPLE_AND_REFINE;
    } else { 
        ret_code = fnft_nsev_base(D_effective, q_preprocessed, r_preprocessed, T, M, contspec, XI, K_ptr,
                    bound_states, normconsts_or_residues, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    
    if (opts->richardson_extrapolation_flag == 1){
        // Allocating memory
        UINT contspec_len = 0;
        if (contspec != NULL && M > 0){
            switch (opts->contspec_type) {
                case nsev_cstype_BOTH:
                    contspec_len = 3*M;
                    break;
                case nsev_cstype_REFLECTION_COEFFICIENT:
                    contspec_len = M;
                    break;
                case nsev_cstype_AB:
                    contspec_len = 2*M;
                    break;
                default:
                    ret_code = E_INVALID_ARGUMENT(opts->contspec_type);
                    goto leave_fun;
            }
            contspec_sub = malloc(contspec_len * sizeof(COMPLEX));
            if (contspec_sub == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
        }
        UINT discspec_len = 0;
        if (kappa == +1 && bound_states != NULL && *K_ptr != 0) {
            K_sub = *K_ptr;
            bound_states_sub = malloc(K_sub * sizeof(COMPLEX));
            discspec_len = K_sub;
            switch (opts->discspec_type) {
                case nsev_dstype_BOTH:
                case nsev_dstype_RESIDUES:
                    discspec_len = 2*K_sub;
                    break;
                case nsev_dstype_NORMING_CONSTANTS:
                    discspec_len = K_sub;
                    break;
                default:
                    ret_code = E_INVALID_ARGUMENT(opts->discspec_type);
                    goto leave_fun;
            }
            normconsts_or_residues_sub = malloc(discspec_len * sizeof(COMPLEX));
            if (normconsts_or_residues_sub == NULL || bound_states_sub == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }
            for (i=0; i<K_sub; i++)
                bound_states_sub[i] = bound_states[i];
        }
        // Preparing q_preprocessed
        UINT method_order;
        method_order = nse_discretization_method_order(opts->discretization);
        if (method_order == 0){
            ret_code =  E_INVALID_ARGUMENT(discretization);
            goto leave_fun;
        }
        Dsub = CEIL(D/2);        
        ret_code = nse_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
                first_last_index, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        
        Tsub[0] = T[0] + first_last_index[0]*eps_t;
        Tsub[1] = T[0] + first_last_index[1]*eps_t;
        const REAL eps_t_sub = (Tsub[1] - Tsub[0])/(Dsub - 1);
        
        // Calling fnft_nsev_base with subsampled signal
        bs_loc_opt = opts->bound_state_localization;
        opts->bound_state_localization = nsev_bsloc_NEWTON;

        ret_code = fnft_nsev_base(Dsub * upsampling_factor, qsub_preprocessed, rsub_preprocessed, Tsub, M, contspec_sub, XI, &K_sub,
                bound_states_sub, normconsts_or_residues_sub, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
        opts->bound_state_localization = bs_loc_opt;
        opts->discspec_type = ds_type_opt;
        
        // Richardson step
        REAL const scl_num = POW(eps_t_sub/eps_t,method_order);
        REAL const scl_den = scl_num - 1.0;
        REAL const dxi = (XI[1]-XI[0])/(M-1);
        if (contspec != NULL && M > 0){
            for (i=0; i<M; i++){
                if (FABS(XI[0]+dxi*i) < 0.9*PI/(2.0*eps_t_sub)){
                    for (j=0; j<contspec_len; j+=M)
                        contspec[i+j] = (scl_num*contspec[i+j] - contspec_sub[i+j])/scl_den;
                }
            }
        }
        if (kappa == +1 && bound_states != NULL && *K_ptr != 0 && K_sub != 0) {
            UINT loc = K_sub;
            REAL bs_err_thres = eps_t;
            REAL bs_err = eps_t;
            UINT K = *K_ptr;
            
            
            for (i=0; i<K; i++){
                loc = K_sub;
                bs_err_thres = eps_t;
                for (j=0; j<K_sub; j++){
                    bs_err = CABS(bound_states[i]-bound_states_sub[j])/CABS(bound_states[i]);
                    if (bs_err < bs_err_thres){
                        bs_err_thres = bs_err;
                        loc = j;
                    }
                }
                if (loc < K_sub){
                    bound_states[i] = (scl_num*bound_states[i] - bound_states_sub[loc])/scl_den;
                    if (ds_type_opt == nsev_dstype_RESIDUES || ds_type_opt == nsev_dstype_BOTH){
                        // Computing aprimes from residues and norming constants
                        normconsts_or_residues_reserve[K+i] = normconsts_or_residues_reserve[i]/normconsts_or_residues_reserve[K+i];
                        normconsts_or_residues_sub[K_sub+loc] = normconsts_or_residues_sub[loc]/normconsts_or_residues_sub[K_sub+loc];
                        // Richardson step on aprime
                        normconsts_or_residues_reserve[K+i] = (scl_num*normconsts_or_residues_reserve[K+i] - normconsts_or_residues_sub[loc+K_sub])/scl_den;
                        // Computing residue
                        normconsts_or_residues_reserve[K+i] = normconsts_or_residues_reserve[i]/normconsts_or_residues_reserve[K+i];
                    }
                }
            }
            if (ds_type_opt == nsev_dstype_RESIDUES)
                memcpy(normconsts_or_residues,normconsts_or_residues_reserve+K,K* sizeof(COMPLEX));
            else if(ds_type_opt == nsev_dstype_BOTH)
                memcpy(normconsts_or_residues,normconsts_or_residues_reserve,2*K* sizeof(COMPLEX));
            
        }
    }
    
    leave_fun:
        free(qsub_preprocessed);
        free(rsub_preprocessed);
        free(q_preprocessed);
        free(r_preprocessed);
        free(contspec_sub);
        free(bound_states_sub);
        free(normconsts_or_residues_sub);
        return ret_code;
}

// Auxiliary function: Base routine for fnft_nsev. fnft_nsev preprocesses the signals
// and calls this function with different options as needed. This prevents
// code doubling while being efficient.
static inline INT fnft_nsev_base(
        const UINT D,
        COMPLEX * const q,
        COMPLEX * r,
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
    UINT deg;
    INT W = 0, *W_ptr = NULL;
    INT ret_code = SUCCESS;
    UINT i, upsampling_factor, D_given;
    
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
        return E_INVALID_ARGUMENT(opts);
    
    
    // Determine step size
    // D is interpolated number of samples but eps_t is the step-size
    // corresponding to original number of samples.
    upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    D_given = D/upsampling_factor;
    const REAL eps_t = (T[1] - T[0])/(D_given - 1);
    
    // D should be the effective number of samples in q
    i = nse_fscatter_numel(D, opts->discretization);
    // NOTE: At this stage if i == 0 it means the discretization corresponds
    // to a slow method. Incorrect discretizations will have been checked for
    // in fnft_nsev main
    
    if (i != 0){
    //This corresponds to methods based on polynomial transfer matrix
       // Allocate memory for the transfer matrix. 
        transfer_matrix = malloc(i*sizeof(COMPLEX));
        if (transfer_matrix == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        // Compute the transfer matrix
        if (opts->normalization_flag)
            W_ptr = &W;
        ret_code = nse_fscatter(D, q, eps_t, kappa, transfer_matrix, &deg, W_ptr,
                opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
    }else{
        // These indicate to the functions to follow that the discretization
        // is a method not based on polynomial transfer matrix
        deg = 0;
        W = 0;
    }
    // Compute the continuous spectrum
    if (contspec != NULL && M > 0) {
        ret_code = nsev_compute_contspec(deg, W, transfer_matrix, q, r, T, D, XI, M,
                contspec, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    
    // Compute the discrete spectrum
    if (kappa == +1 && bound_states != NULL) {
        
        // Compute the bound states
        ret_code = nsev_compute_boundstates(D, q, r, deg, transfer_matrix, T,
                eps_t, K_ptr, bound_states, opts);
        CHECK_RETCODE(ret_code, leave_fun);
        
        
        // Norming constants and/or residues)
        if (normconsts_or_residues != NULL && *K_ptr != 0) {

              ret_code = nsev_compute_normconsts_or_residues(D, q, r, T, *K_ptr,
                    bound_states, normconsts_or_residues, opts);
            CHECK_RETCODE(ret_code, leave_fun);
            
        }
    } else if (K_ptr != NULL) {
        *K_ptr = 0;
    }
    
    leave_fun:
        free(transfer_matrix);
        return ret_code;
}

// Auxiliary function for filtering: We assume that bound states must have
// real part in the interval [-re_bound, re_bound].
static inline REAL re_bound(const REAL eps_t, const REAL map_coeff)
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

// Auxiliary function: Computes the bound states.
static inline INT nsev_compute_boundstates(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        const UINT deg,
        COMPLEX * const transfer_matrix,
        REAL const * const T,
        const REAL eps_t,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        fnft_nsev_opts_t * const opts)
{
    REAL degree1step = 0.0 , map_coeff = 2.0;
    UINT K, upsampling_factor, i, j, D_given;
    REAL bounding_box[4] = { NAN };
    COMPLEX * buffer = NULL;
    INT ret_code = SUCCESS;
    COMPLEX * q_tmp = NULL;
    nse_discretization_t discretization;
    
    degree1step = nse_discretization_degree(opts->discretization);
    // degree1step == 0 here indicates a valid slow method. Incorrect 
    // discretizations should have been caught earlier.
    upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    if (degree1step != 0)
        map_coeff = 2/(degree1step);
    D_given = D/upsampling_factor;
    // Localize bound states ...
    switch (opts->bound_state_localization) {
        
        // ... using Newton's method
        case nsev_bsloc_NEWTON:
            
            K = *K_ptr;
            buffer = bound_states; // Store intermediate results directly
            
            // Perform Newton iterations. Initial guesses of bound-states
            // should be in the continuous-time domain.
            
            // Setting 'discretization' as the base method for discretizations based on
            // splitting schemes.
            if (upsampling_factor == 1 && degree1step != 0){
                discretization = nse_discretization_BO;
            }else if(upsampling_factor == 2 && degree1step != 0){
                discretization = nse_discretization_CF4_2;
            }else
                discretization = opts->discretization;
            
            ret_code = nsev_refine_bound_states_newton(D, q, r, T, K, buffer,
                    discretization, opts->niter);
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
            ret_code = nse_z_to_lambda(K, eps_t, buffer, opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);
            
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
        bounding_box[2] = 0;
        // This step is required as q contains scaled values on a
            // non-equispaced grid
            if (upsampling_factor == 1){
                bounding_box[3] = im_bound(D_given, q, T);
            } else {
                q_tmp = malloc(D_given * sizeof(COMPLEX));
                if (q_tmp == NULL) {
                    ret_code = E_NOMEM;
                    goto leave_fun;
                }
                j = 1;
                for (i = 0; i < D_given; i++) {
                    q_tmp[i] = upsampling_factor*q[j];
                    j = j+upsampling_factor;
                }
                bounding_box[3] = im_bound(D_given, q_tmp, T);
            }
        ret_code = misc_filter(&K, buffer, NULL, bounding_box);
        CHECK_RETCODE(ret_code, leave_fun);
        
        ret_code = misc_merge(&K, buffer, SQRT(EPSILON));
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
        free(q_tmp);
        return ret_code;
}

// Auxiliary function: Computes continuous spectrum on a frequency grid
static inline INT nsev_compute_contspec(
        const UINT deg,
        const INT W,
        COMPLEX * const transfer_matrix,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT D,
        REAL const * const XI,
        const UINT M,
        COMPLEX * const result,
        const INT kappa,
        fnft_nsev_opts_t * const opts)
{
    COMPLEX *H11_vals = NULL, *H21_vals = NULL;
    COMPLEX A, V;
    REAL scale;
    REAL phase_factor_rho, phase_factor_a, phase_factor_b;
    INT ret_code = SUCCESS;
    UINT i, offset = 0, upsampling_factor, D_given;
    COMPLEX * scatter_coeffs = NULL, * xi = NULL;

    // Determine step size
    // D is interpolated number of samples but eps_t is the step-size
    // corresponding to original number of samples.
    upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    D_given = D/upsampling_factor;
    const REAL eps_t = (T[1] - T[0])/(D_given - 1);
    const REAL eps_xi = (XI[1] - XI[0])/(M - 1);   
     
    // Build xi-grid which is required for applying boundary conditions
    xi = malloc(M * sizeof(COMPLEX));
    if (xi == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < M; i++)
        xi[i] = XI[0] + eps_xi*i;

    // Allocate memory for transfer matrix values
    H11_vals = malloc(2*M * sizeof(COMPLEX));
    if (H11_vals == NULL){
        return E_NOMEM;
        goto leave_fun;}
    H21_vals = H11_vals + M;
    
    // If the discretization is a slow method then there should be no transfer_matrix
    if (deg == 0 && transfer_matrix == NULL && W == 0){
        
        // Allocate memory for call to nse_scatter_matrix
        scatter_coeffs = malloc(4 * M * sizeof(COMPLEX));
        if (scatter_coeffs == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }

        ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, M,
                xi, scatter_coeffs, opts->discretization, 0);
        CHECK_RETCODE(ret_code, leave_fun);
        
        // This is necessary because nse_scatter_matrix to ensure
        // boundary conditions can be applied using common code for slow
        // methods and polynomial transfer matrix based methods.
        for (i = 0; i < M; i++){
            H11_vals[i] = scatter_coeffs[i*4];
            H21_vals[i] = scatter_coeffs[i*4+2];
        }

    }else{
        // Prepare the use of the chirp transform. The entries of the transfer
        // matrix that correspond to a and b will be evaluated on the frequency
        // grid xi(i) = XI1 + i*eps_xi, where i=0,...,M-1. Since
        // z=exp(2.0*I*XI*eps_t/degree1step), we find that the z at which z the transfer
        // matrix has to be evaluated are given by z(i) = 1/(A * V^-i), where:
        V = eps_xi;
        ret_code = nse_lambda_to_z(1, eps_t, &V, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        A = -XI[0];
        ret_code = nse_lambda_to_z(1, eps_t, &A, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        
        ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, H11_vals);
        CHECK_RETCODE(ret_code, leave_fun);
        
        ret_code = poly_chirpz(deg, transfer_matrix+2*(deg+1), A, V, M,
                H21_vals);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    // Compute the continuous spectrum
    switch (opts->contspec_type) {
        
        case nsev_cstype_BOTH:
            
            offset = M;
            // fall through
            
        case nsev_cstype_REFLECTION_COEFFICIENT:
            
            
            ret_code = nse_phase_factor_rho(eps_t, T[1], &phase_factor_rho,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);
            
            
            for (i = 0; i < M; i++) {
                if (H11_vals[i] == 0.0){
                    return E_DIV_BY_ZERO;
                    goto leave_fun;
                }
                result[i] = H21_vals[i] * CEXP(I*xi[i]*phase_factor_rho) / H11_vals[i];
            }
            
            if (opts->contspec_type == nsev_cstype_REFLECTION_COEFFICIENT)
                break;
            // fall through
            
        case nsev_cstype_AB:
            
            scale = POW(2.0, W); // needed since the transfer matrix might
            // have been scaled by nse_fscatter. W == 0 for slow methods.
            
            ret_code = nse_phase_factor_a(eps_t, D_given, T, &phase_factor_a,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);
            
            ret_code = nse_phase_factor_b(eps_t, D_given, T, &phase_factor_b,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);
            
            
            for (i = 0; i < M; i++) {
                result[offset + i] = H11_vals[i] * scale * CEXP(I*xi[i]*phase_factor_a);
                result[offset + M + i] = H21_vals[i] * scale * CEXP(I*xi[i]*phase_factor_b);
            }
            
            break;
            
        default:
            
            ret_code = E_INVALID_ARGUMENT(opts->contspec_type);
            goto leave_fun;
    }
    
    leave_fun:
        free(H11_vals);
        free(scatter_coeffs);
        free(xi);
        return ret_code;
}



// Auxiliary function: Computes the norming constants and/or residues
// using slow scattering schemes
static inline INT nsev_compute_normconsts_or_residues(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT K,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        fnft_nsev_opts_t * const opts)
{
    
    COMPLEX *a_vals = NULL, *aprime_vals = NULL;
    UINT i, offset = 0;
    INT ret_code = SUCCESS;
    nse_discretization_t discretization;
    // Check inputs
    if (K == 0) // no bound states 
        return SUCCESS;
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    
    a_vals = malloc(K * sizeof(COMPLEX));
    aprime_vals = malloc(K * sizeof(COMPLEX));
    if (a_vals == NULL || aprime_vals == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    const UINT upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    // Setting discretization as the base method for discretizations based on
    // splitting schemes.
    REAL degree1step = nse_discretization_degree(opts->discretization);
    // degree1step == 0 here implies method not based on polynomial
    // transfer matrix.
    if (upsampling_factor == 1 && degree1step != 0)
        discretization  = nse_discretization_BO;
    else if (upsampling_factor == 2 && degree1step != 0)
        discretization  = nse_discretization_CF4_2;
    else
        discretization = opts->discretization;
    
    ret_code = nse_scatter_bound_states(D, q, r, T, K,
            bound_states, a_vals, aprime_vals, normconsts_or_residues, discretization, 0);
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
        free(a_vals);
        free(aprime_vals);
        return ret_code;
}


// Auxiliary function: Refines the bound-states using Newtons method
static inline INT nsev_refine_bound_states_newton(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT K,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        const UINT niter)
{
    INT ret_code = SUCCESS;
    UINT i, iter, upsampling_factor, D_given;
    COMPLEX a_val, b_val, aprime_val, error;
    REAL eprecision = EPSILON * 100;
    REAL re_bound_val, im_bound_val = NAN;
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
    
    upsampling_factor = nse_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(discretization);
        goto leave_fun;
    }
    D_given = D/upsampling_factor;
    const REAL eps_t = (T[1] - T[0])/(D_given - 1);
    
    im_bound_val = upsampling_factor*upsampling_factor*im_bound(D, q, T);
    if (im_bound_val == NAN){
        ret_code = E_OTHER("Upper bound on imaginary part of bound states is NaN");
        CHECK_RETCODE(ret_code, leave_fun);
    }
    re_bound_val = re_bound(eps_t, 2.0);
    
    // Perform iterations of Newton's method
    for (i = 0; i < K; i++) {
        iter = 0;
        do {
            // Compute a(lam) and a'(lam) at the current root
            ret_code = nse_scatter_bound_states(D, q, r, T, 1,
                    bound_states + i, &a_val, &aprime_val, &b_val, discretization, 1);
            if (ret_code != SUCCESS){
                ret_code = E_SUBROUTINE(ret_code);
                CHECK_RETCODE(ret_code, leave_fun);
            }
            // Perform Newton updates: lam[i] <- lam[i] - a(lam[i])/a'(lam[i])
            
//             printf("l=%1.5e+%1.5e;\t",CREAL(bound_states[i]),CIMAG(bound_states[i]));
//             printf("a=%1.5e+%1.5e;\t",CREAL(a_val),CIMAG(a_val));
//             printf("aprime=%1.5e+%1.5e;\n",CREAL(aprime_val),CIMAG(aprime_val));
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
    
    leave_fun:
        return ret_code;
}

