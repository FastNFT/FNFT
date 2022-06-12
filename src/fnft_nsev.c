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
 * Sander Wahls (TU Delft) 2017-2018, 2022.
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Marius Brehler (TU Dortmund) 2018.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_nsev.h"
#include "fnft__fft_wrapper.h"
#include "fnft__delaunay_triangulation.h"


static fnft_nsev_opts_t default_opts = {
    .bound_state_filtering = nsev_bsfilt_AUTO,
    .bounding_box[0] = -FNFT_INF,
    .bounding_box[1] = FNFT_INF,
    .bounding_box[2] = -FNFT_INF,
    .bounding_box[3] = FNFT_INF,
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
        UINT niter,
        REAL const * const bounding_box);

static inline INT nsev_refine_bound_states_muller(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT K,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        const UINT niter,
        REAL const * const bounding_box);

static inline INT nsev_localize_bound_states_pjt(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        UINT * const K_ptr,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        UINT max_evals,
        REAL const * const bounding_box);

static inline INT nsev_localize_bound_states_GRPF(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        UINT * const K_ptr,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        UINT max_evals,
        REAL const * const bounding_box);

static inline REAL re_bound(const REAL eps_t, const REAL map_coeff);

static inline REAL im_bound(const UINT D, COMPLEX const * const q,
        REAL const * const T);

static UINT FindNextNode(REAL *NodesCoordX, REAL *NodesCoordY, UINT PrevNode, 
        UINT RefNode, UINT *TempNodes, UINT NrOfTempNodes);

/**
 * Fast nonlinear Fourier transform for the nonlinear Schroedinger
 * equation with vanishing boundary conditions.
 * This function takes care of the necessary preprocessing (for example
 * signal resampling or subsampling) based on the options before calling
 * fnft_nsev_base. If the richardson_extrapolation_flag is set, this function
 * calls fnft_nsev_base with half of the samples and then performs
 * Richardson extrapolation on the spectrum.
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
    fft_wrapper_plan_t plan_fwd = fft_wrapper_safe_plan_init();
    COMPLEX *buf0 = NULL, *buf1 = NULL;
    REAL * FTq = NULL;

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
    if ( !(opts->bounding_box[0] <= opts->bounding_box[1]) //!(...) ensures error with NANs
    || !(opts->bounding_box[2] <= opts->bounding_box[3]) )
        return E_INVALID_ARGUMENT(opts->bounding_box);
    
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
            if ((opts->bound_state_localization != nsev_bsloc_NEWTON && 
                    opts->bound_state_localization != nsev_bsloc_MULLER &&
                    opts->bound_state_localization != nsev_bsloc_PJT &&
                    opts->bound_state_localization != nsev_bsloc_GRPF) && kappa == +1){
                ret_code = E_INVALID_ARGUMENT(opts->bound_state_localization);
                goto leave_fun;
            }
            break;
        default: // Unknown discretization
            return E_INVALID_ARGUMENT(opts->discretization);
    }
    
    // Some higher-order discretizations require samples on a non-equidistant grid
    // while others require derivatives which are computed in the form of
    // finite-differences. The input array q has D samples corresponding to
    // an equidistant grid. The upsampling_factor*D gives the effective
    // number of samples for the chosen discretization. The effective number
    // of samples is required for the auxiliary function calls.
    upsampling_factor = nse_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }

    D_effective = D * upsampling_factor;
    
    // Determine step size
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    
    // Set-up bounding_box based on choice of filtering
    if (opts->bound_state_filtering == nsev_bsfilt_BASIC) {
        opts->bounding_box[0] = -INFINITY;
        opts->bounding_box[1] = INFINITY;
        opts->bounding_box[2] = 0.0;
        opts->bounding_box[3] = INFINITY;
    }else if (opts->bound_state_filtering == nsev_bsfilt_AUTO) {
        REAL degree1step = 0.0 , map_coeff = 2.0;
        degree1step = nse_discretization_degree(opts->discretization);
        // degree1step == 0 here indicates a valid slow method. Incorrect
        // discretizations should have been caught earlier.
        if (degree1step != 0)
            map_coeff = 2/(degree1step);
        
        opts->bounding_box[1] = re_bound(eps_t, map_coeff);
        opts->bounding_box[0] = -opts->bounding_box[1];
        opts->bounding_box[2] = 0;
        opts->bounding_box[3] = im_bound(D, q, T);
        
        
        // Fixing search box and some other parameters
        REAL Cq = 1e-3;
        // Allocate memory
        const UINT lenmem = D * sizeof(COMPLEX);
        buf0 = fft_wrapper_malloc(lenmem);
        buf1 = fft_wrapper_malloc(lenmem);
        if (buf0 == NULL || buf1 == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        ret_code = fft_wrapper_create_plan(&plan_fwd, D, buf0, buf1, -1);
        CHECK_RETCODE(ret_code, leave_fun);
        
        // The nonlinear Fourier spectrum corresponds to the flipped and scaled
        // linear Fourier spectrum
        for (i = 0; i < D; i++)
            buf0[i] = -CONJ(q[i]);
        
        ret_code = fft_wrapper_execute_plan(plan_fwd, buf0, buf1);
        CHECK_RETCODE(ret_code, leave_fun);
        
        FTq = calloc(D, sizeof(REAL));
        if (FTq == NULL) {
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        // Computing absolute value squared and applying fftshift
        UINT D_shift = (UINT)CEIL(D/2);
        for (i = 0; i < D_shift-1; i++)
            FTq[i] = CABS(buf1[i+D_shift])*CABS(buf1[i+D_shift]);
        for (i = 0; i < D_shift; i++)
            FTq[i+D_shift-1] = CABS(buf1[i])*CABS(buf1[i]);
        
        REAL FTq_max = misc_max(D,FTq);
        
        REAL eps_xi = (2*0.9*PI/(2*eps_t))/(D-1);
        REAL tmp = 0;
        for (i = D - 1; i > 0; i--){
            tmp = FTq[i]/FTq_max;
            if (tmp > Cq){
                opts->bounding_box[1] = (0.9*PI/(2*eps_t))-eps_xi*(D-i);
                break;
            }
        }
        opts->bounding_box[1] = (0.9*PI/(2*eps_t))-eps_xi*(D-i);
        for (i = 0; i < D; i++){
            if ((FTq[i]/FTq_max)>Cq){
                opts->bounding_box[0] = -(0.9*PI/(2*eps_t))+eps_xi*(i);
                break;
            }
        }

    }else if (opts->bound_state_filtering == nsev_bsfilt_NONE){
        opts->bounding_box[0] = -INFINITY;
        opts->bounding_box[1] = INFINITY;
        opts->bounding_box[2] = -INFINITY;
        opts->bounding_box[3] = INFINITY;
    }
    
    // If Richardson extrapolation is not requested or if it is requested but
    // opts->discspec_type is nsev_dstype_NORMCONSTS,then normconsts_or_residues_reserve
    // acts as just a place holder for normconsts_or_residues. This is equivalent
    // to calling the auxiliary functions with normconsts_or_residues instead of
    // normconsts_or_residues_reserve.
    normconsts_or_residues_reserve = normconsts_or_residues;
    // If Richardson extrapolation is requested and opts->discspec_type
    // is nsev_dstype_RESIDUES, *K_ptr*2  memory is allocated
    // for normconsts_or_residues_reserve. This overwrites the previous
    // assignment of normconsts_or_residues to normconsts_or_residues_reserve.
    // Richardson extrpolation is applied on normconsts(b) and aprimes separately
    // before combining the results to get the residues=normconsts/aprimes.
    // Hence double the memory is necessary even if the user requests only residues.
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
    
    // Some higher-order discretizations require samples on a non-equidistant grid
    // while others require derivatives. The preprocessing function computes
    // the required non-equidistant samples using bandlimited interpolation and
    // the required derivatives using finite differences. The nse_discretization_preprocess_signal
    // function also performs some further discretization specific processing. 
    // Preprocessing takes care of computing things which are required by all
    // the auxiliary functions thus helping efficiency.
    Dsub = D;
    ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
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

        ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
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
        UINT method_order;
        method_order = nse_discretization_method_order(opts->discretization);
        if (method_order == 0){
            ret_code =  E_INVALID_ARGUMENT(discretization);
            goto leave_fun;
        }

        // The signal q is now subsampled(approx. half the samples) and 
        // preprocessed as required for the discretization. This is
        // required for obtaining a second approximation of the spectrum
        // which will be used for Richardson extrapolation.
        Dsub = CEIL(D/2);        
        ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
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
        
        // Richardson extrapolation of the continuous spectrum
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
        // Richardson extrapolation of the discrete spectrum
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
        fft_wrapper_destroy_plan(&plan_fwd);
        fft_wrapper_free(buf0);
        fft_wrapper_free(buf1);
        free(FTq);
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
    REAL degree1step = 0.0;
    UINT K, upsampling_factor, D_given;
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
                    discretization, opts->niter, opts->bounding_box);
            CHECK_RETCODE(ret_code, leave_fun);           
            break;
            
        // ... using Muller's method
        case nsev_bsloc_MULLER:

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

            ret_code = nsev_refine_bound_states_muller(D, q, r, T, K, buffer,
                    discretization, opts->niter, opts->bounding_box);
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
            ret_code = nse_discretization_z_to_lambda(K, eps_t, buffer, opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            break;

        // ... using the phase jump tracking based root finding
        case nsev_bsloc_PJT:
            
            K = *K_ptr;
            buffer = bound_states; // Store intermediate results directly
            // Setting 'discretization' as the base method for discretizations based on
            // splitting schemes.
            if (upsampling_factor == 1 && degree1step != 0){
                discretization = nse_discretization_BO;
            }else if(upsampling_factor == 2 && degree1step != 0){
                discretization = nse_discretization_CF4_2;
            }else
                discretization = opts->discretization;
            
            ret_code = nsev_localize_bound_states_pjt(D, q, r, T, &K, buffer,
                    discretization, 8*D_given, opts->bounding_box);
            break;
            
        // ... using the phase jump tracking based root finding
        case nsev_bsloc_GRPF:
            
            K = *K_ptr;
            buffer = bound_states; // Store intermediate results directly
            // Setting 'discretization' as the base method for discretizations based on
            // splitting schemes.
            if (upsampling_factor == 1 && degree1step != 0){
                discretization = nse_discretization_BO;
            }else if(upsampling_factor == 2 && degree1step != 0){
                discretization = nse_discretization_ES4;
            }else
                discretization = opts->discretization;
            
            ret_code = nsev_localize_bound_states_GRPF(D, q, r, T, &K, buffer,
                    discretization, 5*D_given, opts->bounding_box);
            break;
            
            
        default:

            return E_INVALID_ARGUMENT(opts->bound_state_localization);
    }
    
    // Filter bound states
    if (opts->bound_state_filtering != nsev_bsfilt_NONE) {    
        
        ret_code = misc_filter(&K, buffer, NULL, opts->bounding_box);
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
        ret_code = nse_discretization_lambda_to_z(1, eps_t, &V, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        A = -XI[0];
        ret_code = nse_discretization_lambda_to_z(1, eps_t, &A, opts->discretization);
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
            
            ret_code = nse_discretization_phase_factor_rho(eps_t, T[1], &phase_factor_rho,opts->discretization);
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
            
            // Calculating the discretization specific phase factors.
            ret_code = nse_discretization_phase_factor_a(eps_t, D_given, T, &phase_factor_a,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            ret_code = nse_discretization_phase_factor_b(eps_t, D_given, T, &phase_factor_b,opts->discretization);
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
            bound_states, a_vals, aprime_vals, normconsts_or_residues, discretization);
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
        const UINT niter,
        REAL const * const bounding_box)
{
    INT ret_code = SUCCESS;
    UINT i, iter;
    COMPLEX a_val, aprime_val, error;
    REAL eprecision = EPSILON * 100;
    
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
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
    || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);

    // Perform iterations of Newton's method
    for (i = 0; i < K; i++) {
        iter = 0;
        do {
            // Compute a(lam) and a'(lam) at the current root
            ret_code = nse_scatter_bound_states(D, q, r, T, 1,
                    bound_states + i, &a_val, &aprime_val, NULL, discretization);
            if (ret_code != SUCCESS){
                ret_code = E_SUBROUTINE(ret_code);
                CHECK_RETCODE(ret_code, leave_fun);
            }
            // Perform some checks
            if (a_val == 0.0) // we found a zero, stop here because otherwise the
                break; // next line will cause an error if it is of higher order
            if (aprime_val == 0.0)
                return E_DIV_BY_ZERO;

            // Perform Newton updates: lam[i] <- lam[i] - a(lam[i])/a'(lam[i])
            error = a_val / aprime_val;
            bound_states[i] -= error;
            iter++;
            if (CIMAG(bound_states[i]) > bounding_box[3]
                    || CREAL(bound_states[i]) > bounding_box[1]
                    || CREAL(bound_states[i]) < bounding_box[0]
                    || CIMAG(bound_states[i]) < bounding_box[2])
                break;

        } while (CABS(error) > eprecision && iter < niter);
    }
    leave_fun:
        return ret_code;
}

// Auxiliary function: Refines the bound-states using Muller's method
static inline INT nsev_refine_bound_states_muller(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        const UINT K,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        const UINT niter,
        REAL const * const bounding_box)
{
    INT ret_code = SUCCESS;
    UINT i, iter;
    COMPLEX a_val[4], lam[4], valA, valB, valC, ratio, disc, den1, den2;
    REAL eprecision = EPSILON * 100;
    
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
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
    || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);

    // Perform iterations of Muller's method
    for (i = 0; i < K; i++) {
        lam[1] = bound_states[i];
        lam[0] = bound_states[i] + bound_states[i]/115;
        lam[2] = bound_states[i] - bound_states[i]/85;
        for  (iter = 0; iter < niter; iter++){
            // Compute a(lam)
            ret_code = nse_scatter_bound_states(D, q, r, T, 3,
                    &lam[0], &a_val[0], NULL, NULL, discretization);
            if (ret_code != SUCCESS){
                ret_code = E_SUBROUTINE(ret_code);
                CHECK_RETCODE(ret_code, leave_fun);
            }           
            
            ratio = (lam[2]-lam[1])/(lam[1]-lam[0]);
            valA = ratio*a_val[2] - ratio*(1+ratio)*a_val[1] + ratio*ratio*a_val[0];
            valB = (2*ratio + 1)*a_val[2] - (1 + ratio)*(1 + ratio)*a_val[1] + ratio*ratio*a_val[0];            
            valC = (1 + ratio)*a_val[2];
            
            if ( valA != 0 ){                
                disc = valB*valB - 4*valA*valC;
                den1 = (valB + CSQRT(disc));
                den2 = (valB - CSQRT(disc));
                
                if (CABS(den1) < CABS(den2))
                    lam[3] = lam[2] - (lam[2] - lam[1])*(2*valC/den2);
                else
                    lam[3] = lam[2] - (lam[2] - lam[1])*(2*valC/den1);
                
            }else if (valB != 0)
                lam[3] = lam[2] - (lam[2] - lam[1])*(2*valC/valB);
            else{
                break;
            }
            
            ret_code = nse_scatter_bound_states(D, q, r, T, 1,
                    &lam[3], &a_val[3], NULL, NULL, discretization);
            if (ret_code != SUCCESS){
                ret_code = E_SUBROUTINE(ret_code);
                CHECK_RETCODE(ret_code, leave_fun);
            }
            
            if (CABS(lam[3] - lam[2]) < eprecision || CABS(a_val[3]) < eprecision){
                lam[2] = lam[3];
                break;
            }
            
            lam[0] = lam[1];
            lam[1] = lam[2];
            lam[2] = lam[3];
            
            a_val[0] = a_val[1];
            a_val[1] = a_val[2];
            a_val[2] = a_val[3];
            
            if (CIMAG(lam[2]) > bounding_box[3]
                    || CREAL(lam[2]) > bounding_box[1]
                    || CREAL(lam[2]) < bounding_box[0]
                    || CIMAG(lam[2]) < bounding_box[2])
                break;            
        }
        bound_states[i] = lam[2];
    }
    leave_fun:
        return ret_code;
}

// Auxiliary function: Finds bound states by tracking the phase-jumps
static inline INT nsev_localize_bound_states_pjt(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        UINT * const K_ptr,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        UINT max_evals,
        REAL const * const bounding_box)
{
    INT ret_code = SUCCESS;
    UINT i, evals = 0;
    COMPLEX *xi = NULL;
    COMPLEX *a_vals = NULL;
    COMPLEX * jump_points = NULL;
    REAL * phi_jump_points = NULL;
    
    // Check inputs
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
    || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);
    
    if (max_evals == 0)
        max_evals = 8*D;
    REAL bounding_box_local[4] = { NAN };
    UINT K = *K_ptr;
    // D is interpolated number of samples but eps_t is the step-size
    // corresponding to original number of samples.
    UINT upsampling_factor = nse_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(discretization);
        goto leave_fun;
    }
    UINT D_given = D/upsampling_factor;    
    
    if (bounding_box[0] < -5)
        bounding_box_local[0] = -5;
    else
        bounding_box_local[0] = bounding_box[0];
    
    if (bounding_box[1] > 5)
        bounding_box_local[1] = 5;
    else
        bounding_box_local[1] = bounding_box[1];
    
    if (bounding_box[2] < 0)
        bounding_box_local[2] = 0;
    else
        bounding_box_local[2] = bounding_box[2];
    
    if (LOG(0.9*DBL_MAX)/((T[1]-T[0])/2) < bounding_box[3])
        bounding_box_local[3] = LOG(0.9*DBL_MAX)/((T[1]-T[0])/2); //TODO
    else
        bounding_box_local[3] = bounding_box[3];
    
    jump_points = malloc(K * sizeof(COMPLEX));
    phi_jump_points = malloc(K * sizeof(REAL));    
    if (jump_points == NULL || phi_jump_points == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    
    REAL Tol = 1e-3;
    
    // Fixing step-size for real-axis
    UINT M = (bounding_box_local[1] - bounding_box_local[0])/Tol;
    if (M > D_given)
        M = D_given;
    
    printf("M=%ld\n",M);
    printf("bbox[0]=%f\n",bounding_box_local[0]);
    printf("bbox[1]=%f\n",bounding_box_local[1]);
    
    xi = malloc(M * sizeof(COMPLEX));
    a_vals = malloc(M * sizeof(COMPLEX));
    if (xi == NULL || a_vals == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    REAL eps_xi = (bounding_box_local[1] - bounding_box_local[0])/(M - 1);
    for (i = 0; i < M; i++)
        xi[i] = bounding_box_local[0] + i*eps_xi;
    
    ret_code = nse_scatter_bound_states(D, q, r, T, M,
            xi, a_vals, NULL, NULL, discretization);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    evals += M;
    
    K = 0;
    for (i = 1; i < M; i++){
        if ((CARG(a_vals[i])*CARG(a_vals[i-1]))<0 && FABS(CARG(a_vals[i])-CARG(a_vals[i-1]))>(1.3*PI)){
            jump_points[K] = (xi[i]+xi[i-1])/2;
            phi_jump_points[K]  = 0;
            K++;
        }
    }
    
    REAL hG = 0.5*misc_min_dist(K,jump_points);// step-size for remaining contour
    if ((bounding_box_local[3]+bounding_box_local[1]-bounding_box_local[0])*0.01 < hG)
        hG = (bounding_box_local[3]+bounding_box_local[1]-bounding_box_local[0])*0.01;
    
    // Right boundary
    M = (bounding_box_local[3] - bounding_box_local[2])/hG;
    xi = malloc(M * sizeof(COMPLEX));
    a_vals = malloc(M * sizeof(COMPLEX));
    if (xi == NULL || a_vals == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < M; i++)
        xi[i] = bounding_box_local[1] + (bounding_box_local[2]+i*hG)*I;
    
    ret_code = nse_scatter_bound_states(D, q, r, T, M,
            xi, a_vals, NULL, NULL, discretization);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    evals += M;
    
    for (i = 1; i < M; i++){
        if ((CARG(a_vals[i])*CARG(a_vals[i-1]))<0 && FABS(CARG(a_vals[i])-CARG(a_vals[i-1]))>(1.3*PI)){
            jump_points[K] = (xi[i]+xi[i-1])/2;
            phi_jump_points[K]  = PI/2;
            K++;
        }
    }
    
    // Upper boundary
    M = (bounding_box_local[1] - bounding_box_local[0])/hG; 
    xi = malloc(M * sizeof(COMPLEX));
    a_vals = malloc(M * sizeof(COMPLEX));
    if (xi == NULL || a_vals == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < M; i++)
        xi[i] = bounding_box_local[3]*I + (bounding_box_local[0]+i*hG);
    
    ret_code = nse_scatter_bound_states(D, q, r, T, M,
            xi, a_vals, NULL, NULL, discretization);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    evals += M;
    
    for (i = 1; i < M; i++){
        if ((CARG(a_vals[i])*CARG(a_vals[i-1]))<0 && FABS(CARG(a_vals[i])-CARG(a_vals[i-1]))>(1.3*PI)){
            jump_points[K] = (xi[i]+xi[i-1])/2;
            phi_jump_points[K]  = PI;
            K++;
        }
    }
    
    // Left boundary
    M = (bounding_box_local[3] - bounding_box_local[2])/hG; 
    xi = malloc(M * sizeof(COMPLEX));
    a_vals = malloc(M * sizeof(COMPLEX));
    if (xi == NULL || a_vals == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i = 0; i < M; i++)
        xi[i] = bounding_box_local[0] + (bounding_box_local[2]+i*hG)*I;
    
    ret_code = nse_scatter_bound_states(D, q, r, T, M,
            xi, a_vals, NULL, NULL, discretization);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
    evals += M;
    
    for (i = 1; i < M; i++){
        if ((CARG(a_vals[i])*CARG(a_vals[i-1]))<0 && FABS(CARG(a_vals[i])-CARG(a_vals[i-1]))>(1.3*PI)){
            jump_points[K] = (xi[i]+xi[i-1])/2;
            phi_jump_points[K]  = 3*PI/2;
            K++;
        }
    }
    
    REAL hg = (1.0/15.0)*misc_min_dist(K,jump_points);
//     if (hg < Tol/2)
//         hg = Tol/2;
    printf("hg=%f\n",hg);
   COMPLEX curr_jump_point;
   COMPLEX jump_point_neighbors[2] = {0};
   UINT quadrants[2];
   REAL angles[2] = {0}, phi_mid, m, hgc;
   for (i = 0; i < K; i++){ 
    hgc = 4*hg;
    curr_jump_point = jump_points[i];
    while (hgc>hg){
        jump_point_neighbors[0] = curr_jump_point-hgc*CEXP(I*phi_jump_points[i]);
        jump_point_neighbors[1] = curr_jump_point+hgc*CEXP(I*phi_jump_points[i]);        
        ret_code = nse_scatter_bound_states(D, q, r, T, 2,
                &jump_point_neighbors[0], a_vals, NULL, NULL, discretization);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        
        evals += 2;
        
       ret_code = misc_quadrant(2,a_vals,&quadrants[0],&angles[0]);
       if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        
        phi_mid = (PI/2)*ROUND(((angles[0]+angles[1])/2)/(PI/2));
        m = (angles[1]-angles[0])/hgc;
        jump_point_neighbors[1] = jump_point_neighbors[1]-CEXP(I*phi_jump_points[i])*(angles[1]-(phi_mid+m*0.5*hgc))/m;
        jump_point_neighbors[0] = jump_point_neighbors[0]+CEXP(I*phi_jump_points[i])*((phi_mid-m*0.5*hgc)-angles[0])/m;
        curr_jump_point = (jump_point_neighbors[0]+jump_point_neighbors[1])/2;
        hgc = hgc/2;         
    }
    if (quadrants[0] == quadrants[1]){
        ret_code = E_ASSERTION_FAILED;
        goto leave_fun;
    }
    jump_points[i] = curr_jump_point;
   }

   misc_print_buf(K,jump_points,"g");
   // Actual phase jump tracking
//    hg = 3*hg;
   REAL phic, hgs, ang_diff = PI/60.0;
   COMPLEX jump_point_neighbors_guesses[2] = {0}, new_jump_point;
   UINT quadrants_guess[2];
   UINT sustep = 0; // Number of successful steps
   UINT rstep = 0; // Number of direction change attempts.
   for (i = 0; i < K; i++){
       hgc = hg;
       hgs = hg;       
       jump_point_neighbors[0] = jump_points[i]-0.5*hgc*CEXP(I*phi_jump_points[i]);
       jump_point_neighbors[1] = jump_points[i]+0.5*hgc*CEXP(I*phi_jump_points[i]);
       
       ret_code = nse_scatter_bound_states(D, q, r, T, 2,
               &jump_point_neighbors[0], a_vals, NULL, NULL, discretization);
       if (ret_code != SUCCESS){
           ret_code = E_SUBROUTINE(ret_code);
           CHECK_RETCODE(ret_code, leave_fun);
       }
       evals += 2;
       
       ret_code = misc_quadrant(2,a_vals,&quadrants[0],&angles[0]);
       if (ret_code != SUCCESS){
           ret_code = E_SUBROUTINE(ret_code);
           CHECK_RETCODE(ret_code, leave_fun);
       }
       
       phi_mid = (PI/2)*ROUND(((angles[0]+angles[1])/2)/(PI/2));
       m = (angles[1]-angles[0])/hg;
       jump_point_neighbors[1] = jump_point_neighbors[1]-CEXP(I*phi_jump_points[i])*(angles[1]-(phi_mid+m*0.5*hgc))/m;
       jump_point_neighbors[0] = jump_point_neighbors[0]+CEXP(I*phi_jump_points[i])*((phi_mid-m*0.5*hgc)-angles[0])/m;
       curr_jump_point = (jump_point_neighbors[0]+jump_point_neighbors[1])/2;
       phic = phi_jump_points[i] + PI/2;
       if (ang_diff/m > hgc)
           hgc = ang_diff/m;
       hgs = hgc/4;
        printf("hgc=%f\n",hgc);
       while (evals < max_evals && phic != NAN){
           
           jump_point_neighbors_guesses[0] = jump_point_neighbors[0] + hgs*CEXP(I*phic);
           jump_point_neighbors_guesses[1] = jump_point_neighbors[1] + hgs*CEXP(I*phic);
           
           ret_code = nse_scatter_bound_states(D, q, r, T, 2,
                   &jump_point_neighbors_guesses[0], a_vals, NULL, NULL, discretization);
           if (ret_code != SUCCESS){
               ret_code = E_SUBROUTINE(ret_code);
               CHECK_RETCODE(ret_code, leave_fun);
           }
           evals += 2;
           
           ret_code = misc_quadrant(2,a_vals,&quadrants_guess[0],&angles[0]);
           if (ret_code != SUCCESS){
               ret_code = E_SUBROUTINE(ret_code);
               CHECK_RETCODE(ret_code, leave_fun);
           }
           
           if (quadrants[0] != quadrants_guess[0] || quadrants[1] != quadrants_guess[1]){
               sustep = 0;
               if (hgs<hgc){
                   if (rstep == 1){
                       bound_states[i] = curr_jump_point;
                       break;
                   }else{
                       rstep++;
                       if (quadrants[0] != quadrants_guess[0]){
                           jump_point_neighbors[0] = CEXP(I*phic)*hgc + jump_point_neighbors[1];
                           phic = phic - PI/2;
                       }else if (quadrants[1] != quadrants_guess[1]){
                           jump_point_neighbors[1] = CEXP(I*phic)*hgc + jump_point_neighbors[0];
                           phic = phic + PI/2;
                       }else{
                           bound_states[i] = curr_jump_point;
                           break;
                       }
                       curr_jump_point = (jump_point_neighbors[0]+jump_point_neighbors[1])/2;
                   }
               }
               else
                   hgs = hgs/2;
           }else{
               new_jump_point = jump_point_neighbors_guesses[1] - (jump_point_neighbors_guesses[1]-jump_point_neighbors_guesses[0])*(angles[1]-phi_mid)/(angles[1]-angles[0]);
               phic = CARG(new_jump_point-curr_jump_point);
               jump_point_neighbors[0] = CEXP(I*(phic+PI/2))*0.5*hgc + new_jump_point;
               jump_point_neighbors[1] = CEXP(I*(phic-PI/2))*0.5*hgc + new_jump_point;
               curr_jump_point = new_jump_point;
               
               sustep ++; // Successful step
               if (hgs < hgc*16 && sustep == 2)  // If two successful steps the double step-size
                   hgs = hgs*2;
               
           }
       }
   }
//      printf("K=%ld\n",K);
//        printf("evals=%ld\n",evals);
   *K_ptr = K;
   leave_fun:
         free(xi);
         free(a_vals);
         free(jump_points);
         free(phi_jump_points);
        return ret_code;
}

// Auxiliary function: Finds bound states by tracking the phase-jumps
static inline INT nsev_localize_bound_states_GRPF(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX * r,
        REAL const * const T,
        UINT * const K_ptr,
        COMPLEX * bound_states,
        nse_discretization_t discretization,
        UINT NodesMax,
        REAL const * const bounding_box)
{
    INT ret_code = SUCCESS;
    UINT i, j;
    
    // Check inputs
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (bounding_box == NULL)
        return E_INVALID_ARGUMENT(bounding_box);
    if ( !(bounding_box[0] <= bounding_box[1]) //!(...) ensures error with NANs
    || !(bounding_box[2] <= bounding_box[3]) )
        return E_INVALID_ARGUMENT(bounding_box);
    
    REAL bounding_box_local[4] = { NAN };
    UINT NodesShift = 0;
    UINT NrOfElements = 0;
    UINT NrOfValidEdges = 0;
    UINT NrOfCandidateEdges = 0;
    UINT NrOfCandidateNodes = 0;
    UINT NrOfCandidateElements = 0;
    UINT *ArrayOfCandidateElements = NULL;
    UINT NrOfFirstZoneCandidateElements = 0;
    UINT NrOfSecondZoneCandidateElements = 0;
    REAL XCoordOfTempEdgeNode1, XCoordOfTempEdgeNode2;
    REAL YCoordOfTempEdgeNode1, YCoordOfTempEdgeNode2;
    REAL XTempNodeCoord, YTempNodeCoord, TempEdgeLength, DistNodes;
    UINT AddCoordFlag;
    UINT *IDOfCandidateElements = NULL;
    UINT *MultiplicationOfTempEdges = NULL;
    UINT *ContourEdges1 = NULL;
    UINT *ContourEdges2 = NULL;
    UINT *NrOfNodesOfRegions = NULL, *Regions = NULL, *IndexOfNextEdge = NULL, *TempNodes = NULL, *SideFlag = NULL;
    INT *z_m = NULL;
    COMPLEX *z = NULL;
    REAL *NewNodesCoordX = NULL;
    REAL *NewNodesCoordY = NULL;
    UINT NrOfNewNodes = 0;
    UINT NrOfNodes = 0;
    REAL *NodesCoordX = NULL;
    REAL *NodesCoordY = NULL;
    COMPLEX * FuntionValues = NULL;
    COMPLEX * ComplexNodes = NULL;
    UINT * Quadrants = NULL;
    REAL * Angles = NULL;
    UINT *ValidEdges = NULL;
    UINT * Elements1 = NULL, * Elements2 = NULL, * Elements3 = NULL;
    UINT *ValidEdges1 = NULL, *ValidEdges2 = NULL;
    UINT *CandidateEdges1 = NULL, *CandidateEdges2 = NULL;
    UINT *ElementsCount = NULL, *CandidateNodes = NULL, * Temp = NULL;
    UINT *FlagOfFirstZoneCandidateElements = NULL;
    UINT *FlagOfSecondZoneCandidateElements = NULL;
    UINT *TempExtraEdges1 = NULL, *TempExtraEdges2 = NULL;
    UINT DT_built = 0;
    
    
    UINT * NodeOfTriangles1 = NULL;
    UINT * NodeOfTriangles2 = NULL;
    UINT * NodeOfTriangles3 = NULL;
    UINT * EdgesOfTriangles1 = NULL;
    UINT * EdgesOfTriangles2 = NULL;
    UINT * EdgesOfTriangles3 = NULL;
    UINT * Edges1 = NULL;
    UINT * Edges2 = NULL;
    UINT NrOfEdges = 0;
    UINT NrOfTriangles = 0;
    UINT NrOfNodesToRemove = 0;
    
    // D is interpolated number of samples but eps_t is the step-size
    // corresponding to original number of samples.
    UINT upsampling_factor = nse_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(discretization);
        goto leave_fun;
    }
    UINT D_given = D/upsampling_factor;
    if (NodesMax == 0)
        NodesMax = 4*D_given;
    
    bounding_box_local[1] = bounding_box[1];
    bounding_box_local[0] = bounding_box[0];
    
    if (bounding_box[2] < 0)
        bounding_box_local[2] = 0;
    else
        bounding_box_local[2] = bounding_box[2];
    
    if (LOG(0.9*DBL_MAX)/((T[1]-T[0])/2) < bounding_box[3])
        bounding_box_local[3] = LOG(0.9*DBL_MAX)/((T[1]-T[0])/2); //TODO
    else
        bounding_box_local[3] = bounding_box[3];
    
//     printf("bbox[0]=%f\n",bounding_box_local[0]);
//     printf("bbox[1]=%f\n",bounding_box_local[1]);
//     printf("bbox[2]=%f\n",bounding_box_local[2]);
//     printf("bbox[3]=%f\n",bounding_box_local[3]);

    ret_code =  delaunay_triangulation_rect_dom(bounding_box_local[0], bounding_box_local[1],
            bounding_box_local[2], bounding_box_local[3],
            (bounding_box_local[3] - bounding_box_local[2])/5,
            &NrOfNewNodes, &NewNodesCoordX, &NewNodesCoordY);
    CHECK_RETCODE(ret_code, leave_fun);
//     printf("NrOfNewNodes=%ld \n",NrOfNewNodes);

    fnft__delaunay_triangulation_data_t * DT_ptr = NULL;
    fnft__delaunay_triangulation_data_t DT;
    DT_ptr = &DT;
    DT_ptr->NrOfNodes = 0;
    DT_ptr->NodesMax = 0;    
    UINT status_flag = 0;
    // Setting all poniters to NULL so that they can be
    // freed if allocation fails for some pointer
    // before it.
    DT_ptr->NodeOfTriangles1 = NULL;
    DT_ptr->NodeOfTriangles2 = NULL;
    DT_ptr->NodeOfTriangles3 = NULL;
    DT_ptr->EdgesOfTriangles1 = NULL;
    DT_ptr->EdgesOfTriangles2 = NULL;
    DT_ptr->EdgesOfTriangles3 = NULL;
    DT_ptr->Edges1 = NULL;
    DT_ptr->Edges2 = NULL;
    DT_ptr->NodesCoordX = NULL;
    DT_ptr->NodesCoordY = NULL;
    DT_ptr->XCoordOfCCOfTriangles = NULL;
    DT_ptr->YCoordOfCCOfTriangles = NULL;
    DT_ptr->RadiusOfCCOfTriangles = NULL;
    DT_ptr->StatusOfTriangles = NULL;
    DT_ptr->CheckEdge = NULL;
    DT_built = 1;
    
    NodesCoordX = malloc(NodesMax * sizeof(REAL));
    NodesCoordY = malloc(NodesMax * sizeof(REAL));
    if (NodesCoordX == NULL || NodesCoordX == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    FuntionValues = malloc(NodesMax * sizeof(COMPLEX));
    ComplexNodes = malloc(NodesMax * sizeof(COMPLEX));
    Quadrants = malloc(NodesMax * sizeof(UINT));
    Angles = malloc(NodesMax * sizeof(REAL));
    if (FuntionValues == NULL || Quadrants == NULL 
            || ComplexNodes == NULL || Angles == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ValidEdges = malloc(1*sizeof(UINT));
    if (ValidEdges == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    Elements1 = malloc(1*sizeof(UINT));
    Elements2 = malloc(1*sizeof(UINT));
    Elements3 = malloc(1*sizeof(UINT));
    if (Elements1 == NULL || Elements2 == NULL|| Elements3 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ValidEdges1 = malloc(1*sizeof(UINT));
    ValidEdges2 = malloc(1*sizeof(UINT));
    if (ValidEdges1 == NULL || ValidEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ArrayOfCandidateElements = malloc(1*sizeof(UINT));        
    CandidateEdges1 = malloc(1*sizeof(UINT));
    CandidateEdges2 = malloc(1*sizeof(UINT));
    if (CandidateEdges1 == NULL || CandidateEdges2 == NULL
            || ArrayOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }    
    
    Temp = malloc(NodesMax*sizeof(UINT));
     if (Temp == NULL ){
        ret_code = E_NOMEM;
        goto leave_fun;
    }    
   
    CandidateNodes = malloc(1*sizeof(UINT));
    if (CandidateNodes == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ElementsCount = malloc(1*sizeof(UINT));
    if (ElementsCount == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    FlagOfFirstZoneCandidateElements = malloc(1*sizeof(UINT));
    FlagOfSecondZoneCandidateElements = malloc(1*sizeof(UINT));
    if (FlagOfFirstZoneCandidateElements == NULL || FlagOfSecondZoneCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    TempExtraEdges1 = malloc(1*sizeof(UINT));
    TempExtraEdges2 = malloc(1*sizeof(UINT));
    if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
//     REAL Tol = (bounding_box_local[1]-bounding_box_local[0])/D_given;
//     if ((bounding_box_local[3]-bounding_box_local[2])/D_given < Tol)
//         Tol = (bounding_box_local[3]-bounding_box_local[2])/D_given;
    REAL Tol = 0.9*PI/(T[1]-T[0]);
    ////////////////////////////////////////
    for (UINT iter=1; iter<=50; iter++){
        
//         printf("iter=%ld\n",iter);
//         printf("NrOfNewNodes=%ld \n",NrOfNewNodes);
        if (NrOfNodes+NrOfNewNodes > NodesMax)
            break;
        
        NodesShift = NrOfNodes;
        NrOfNodes = NrOfNodes+NrOfNewNodes;
        
        for (i=0; i<NrOfNewNodes; i++){
            NodesCoordX[NodesShift+i] = NewNodesCoordX[i];
            NodesCoordY[NodesShift+i] = NewNodesCoordY[i];
            ComplexNodes[i] = NewNodesCoordX[i]+I*NewNodesCoordY[i];
        }
        
        ret_code = nse_scatter_bound_states(D, q, r, T, NrOfNewNodes,
                ComplexNodes, FuntionValues+NodesShift, NULL, NULL, discretization);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        misc_quadrant(NrOfNewNodes, FuntionValues+NodesShift,
                Quadrants+NodesShift, Angles+NodesShift);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
//         printf("[");        
//         for (j=0; j<NrOfNodes; j++)
//             printf("%ld, ", Quadrants[j]);
//         printf("]\n");

        ret_code = delaunay_triangulation(DT_ptr, NodesMax, NrOfNewNodes,
                NewNodesCoordX, NewNodesCoordY, &status_flag);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        
        if (status_flag == 1){
            WARN("Duplicate point.");
            break;
        }
        if (status_flag == 2){
            WARN("Ran out of memory for triangles.");
            break;
        }
         if (status_flag == 3){
            WARN("Ran out of memory for Edges.");
            break;
        }
        
//         printf("NrOfNodes=%ld \n",DT_ptr->NrOfNodes);
//         printf("NodesMax=%ld \n",DT_ptr->NodesMax);
//         printf("NrOfTriangles=%ld \n",DT_ptr->NrOfTriangles);
//         printf("NrOfEdges=%ld \n",DT_ptr->NrOfEdges);
        
        NodeOfTriangles1 = DT_ptr->NodeOfTriangles1;
        NodeOfTriangles2 = DT_ptr->NodeOfTriangles2;
        NodeOfTriangles3 = DT_ptr->NodeOfTriangles3;
        Edges1 = DT_ptr->Edges1;
        Edges2 = DT_ptr->Edges2;
        NrOfEdges = DT_ptr->NrOfEdges;
        EdgesOfTriangles1 = DT_ptr->EdgesOfTriangles1;
        EdgesOfTriangles2 = DT_ptr->EdgesOfTriangles2;
        EdgesOfTriangles3 = DT_ptr->EdgesOfTriangles3;
        NrOfTriangles = DT_ptr->NrOfTriangles;
        NrOfNodesToRemove = DT_ptr->NrOfNodesToRemove;
        
        ValidEdges = realloc(ValidEdges,NrOfEdges*sizeof(UINT));
        if (ValidEdges == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfEdges; i++)
            ValidEdges[i] = 0;
        
        NrOfElements = NrOfTriangles;       
        for (i=0; i<NrOfTriangles; i++){
            ValidEdges[EdgesOfTriangles1[i]] = 1;
            ValidEdges[EdgesOfTriangles2[i]] = 1;
            ValidEdges[EdgesOfTriangles3[i]] = 1;
        }
        
        Elements1 = realloc(Elements1,NrOfElements*sizeof(UINT));
        Elements2 = realloc(Elements2,NrOfElements*sizeof(UINT));
        Elements3 = realloc(Elements3,NrOfElements*sizeof(UINT));
        if (Elements1 == NULL || Elements2 == NULL|| Elements3 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfTriangles; i++){
            Elements1[i] = NodeOfTriangles1[i]-NrOfNodesToRemove;
            Elements2[i] = NodeOfTriangles2[i]-NrOfNodesToRemove;
            Elements3[i] = NodeOfTriangles3[i]-NrOfNodesToRemove;
        }
                
//         for (j=0; j<NrOfElements; j++)
//             printf("[%ld,%ld,%ld]\n", Elements1[j],Elements2[j],Elements3[j]);
        
        NrOfValidEdges = 0;
        for (i=0; i<NrOfEdges; i++){
            if (ValidEdges[i] == 1)
                NrOfValidEdges = NrOfValidEdges+1;
        }
//         printf(" NrOfValidEdges=%ld \n", NrOfValidEdges);
        
        ValidEdges1 = realloc(ValidEdges1, NrOfValidEdges*sizeof(UINT));
        ValidEdges2 = realloc(ValidEdges2, NrOfValidEdges*sizeof(UINT));
        if (ValidEdges1 == NULL || ValidEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        j = 0;
        for (i=0; i<NrOfEdges; i++){
            if (ValidEdges[i] == 1){
                ValidEdges1[j] = Edges1[i]-NrOfNodesToRemove;
                ValidEdges2[j] = Edges2[i]-NrOfNodesToRemove;
//                 printf("[%ld,%ld]\n", ValidEdges1[j],ValidEdges2[j]);
                j = j+1;
            }
        }
        // CandidateEdges are edges across whom phase difference is 2
        CandidateEdges1 = realloc(CandidateEdges1, NrOfValidEdges*sizeof(UINT));
        CandidateEdges2 = realloc(CandidateEdges2, NrOfValidEdges*sizeof(UINT));
        if (CandidateEdges1 == NULL || CandidateEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        NrOfCandidateEdges = 0;
        for (i=0; i<NrOfValidEdges; i++){
            if (((Quadrants[ValidEdges1[i]]-Quadrants[ValidEdges2[i]])%4) == 2){
                CandidateEdges1[NrOfCandidateEdges] =  ValidEdges1[i];
                CandidateEdges2[NrOfCandidateEdges] =  ValidEdges2[i];
                NrOfCandidateEdges = NrOfCandidateEdges + 1;
            }
        }
        
        if  (NrOfCandidateEdges == 0){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
//         printf("NrOfCandidateEdges=%ld\n",NrOfCandidateEdges );
        
        for (i=0; i<(DT_ptr->NrOfNodes); i++)
            Temp[i] = 0;
        
        REAL MinCandidateEdgesLengths = INFINITY;
        REAL MaxCandidateEdgesLengths = 0;
        REAL CurrEdgeLength = 0;
        // Length of the CandidateEdges
        for (i=0; i<NrOfCandidateEdges; i++){
            CurrEdgeLength = SQRT((NodesCoordX[CandidateEdges2[i]]-NodesCoordX[CandidateEdges1[i]])*(NodesCoordX[CandidateEdges2[i]]-NodesCoordX[CandidateEdges1[i]])
            +(NodesCoordY[CandidateEdges2[i]]-NodesCoordY[CandidateEdges1[i]])*(NodesCoordY[CandidateEdges2[i]]-NodesCoordY[CandidateEdges1[i]]));
            if (CurrEdgeLength < MinCandidateEdgesLengths)
                MinCandidateEdgesLengths = CurrEdgeLength;
            if (CurrEdgeLength > MaxCandidateEdgesLengths)
                MaxCandidateEdgesLengths = CurrEdgeLength;
            // If a CandidateEdge is longer than Tol then its end points are marked
            if (CurrEdgeLength > Tol){
                Temp[CandidateEdges1[i]] = 1;// TODO may require to set everything to 0 first
                Temp[CandidateEdges2[i]] = 1;
            }
        }
        
        if (MaxCandidateEdgesLengths < Tol){
             printf("Tolerence achieved in interation %ld \n", iter);
            break;
        }
        
        // Collecting all the nodes which are the end points
        // Marking first and then collecting to preventing counting a node twice
        CandidateNodes = realloc(CandidateNodes, NrOfNodes*sizeof(UINT));
        if (CandidateNodes == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        NrOfCandidateNodes = 0;
        for (i=0; i<NrOfNodes; i++){
            if (Temp[i] == 1){
                CandidateNodes[NrOfCandidateNodes] = i;
                NrOfCandidateNodes = NrOfCandidateNodes+1;
             //printf("CandidateNode =%ld\n", i);
            }
        }
//         printf(" NrOfCandidateNodes =%ld\n", NrOfCandidateNodes  );

        // Listing all the triangles that the nodes in CandidateNodes are part
        // of
        NrOfCandidateElements = 0;
        ArrayOfCandidateElements= realloc(ArrayOfCandidateElements, 3*NrOfElements*sizeof(UINT));
        if (ArrayOfCandidateElements == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        ret_code =  delaunay_triangulation_vertexAttachments(NrOfElements,Elements1, Elements2, Elements3,
                NrOfCandidateNodes, CandidateNodes, &NrOfCandidateElements, ArrayOfCandidateElements);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
//         printf(" NrOfCandidateElements =%ld\n", NrOfCandidateElements  );
//  for (i=0; i<NrOfCandidateElements; i++)
//      printf(" CandidateElement =%ld\n", ArrayOfCandidateElements[i]);
//
        // Calculating the number of times a particular triangle shows up in
        // ArrayOfCandidateElements
        ElementsCount = realloc(ElementsCount, NrOfElements*sizeof(UINT));
        if (ElementsCount == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfElements; i++)
            ElementsCount[i] = 0;
        for (i=0; i<NrOfCandidateElements; i++)
            ElementsCount[ArrayOfCandidateElements[i]] = ElementsCount[ArrayOfCandidateElements[i]]+1;
                
        NrOfFirstZoneCandidateElements = 0;
        NrOfSecondZoneCandidateElements = 0;
        
        FlagOfFirstZoneCandidateElements = realloc(FlagOfFirstZoneCandidateElements, NrOfElements*sizeof(UINT));
        FlagOfSecondZoneCandidateElements = realloc(FlagOfSecondZoneCandidateElements, NrOfElements*sizeof(UINT));
        if (FlagOfFirstZoneCandidateElements == NULL || FlagOfSecondZoneCandidateElements == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfElements; i++){
            FlagOfFirstZoneCandidateElements[i] = 0;
            FlagOfSecondZoneCandidateElements[i] = 0;
        }
        
        for (i=0; i<NrOfElements; i++){
            // If a triangle shows up more than once
            if (ElementsCount[i]>1){
                FlagOfFirstZoneCandidateElements[i] = 1;
                NrOfFirstZoneCandidateElements = NrOfFirstZoneCandidateElements+1;
            }
            // If a triangle shows up exactly once
            if (ElementsCount[i] == 1){
                FlagOfSecondZoneCandidateElements[i] = 1;
                NrOfSecondZoneCandidateElements = NrOfSecondZoneCandidateElements+1;
            }
        }
//         printf("NrOfFirstZoneCandidateElements =%ld\n", NrOfFirstZoneCandidateElements);
//         printf("NrOfSecondZoneCandidateElements =%ld\n", NrOfSecondZoneCandidateElements);
//         
        // TODO possibility to used Edges1 and Edges2 directly
        // Vectorizing FirstZoneCandidateElements
        TempExtraEdges1 = realloc(TempExtraEdges1, NrOfFirstZoneCandidateElements*3*sizeof(UINT));
        TempExtraEdges2 = realloc(TempExtraEdges2, NrOfFirstZoneCandidateElements*3*sizeof(UINT));
        if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        j = 0;
        for (i=0; i<NrOfElements; i++){
            if (FlagOfFirstZoneCandidateElements[i] == 1){
                TempExtraEdges1[j*3] = Elements1[i];
                TempExtraEdges2[j*3] = Elements2[i];
                
                TempExtraEdges1[j*3+1] = Elements2[i];
                TempExtraEdges2[j*3+1] = Elements3[i];
                
                TempExtraEdges1[j*3+2] = Elements3[i];
                TempExtraEdges2[j*3+2] = Elements1[i];
                j = j+1;
            }
        }
        
        // Creating newnodes at mid-points of the edges
        NewNodesCoordX = realloc(NewNodesCoordX, 3*(NrOfFirstZoneCandidateElements+NrOfSecondZoneCandidateElements)*sizeof(REAL));
        NewNodesCoordY = realloc(NewNodesCoordY, 3*(NrOfFirstZoneCandidateElements+NrOfSecondZoneCandidateElements)*sizeof(REAL));
        NewNodesCoordX[0] = (NodesCoordX[TempExtraEdges1[0]]+NodesCoordX[TempExtraEdges2[0]])/2;
        NewNodesCoordY[0] = (NodesCoordY[TempExtraEdges1[0]]+NodesCoordY[TempExtraEdges2[0]])/2;
        NrOfNewNodes = 1;
        
        for (i=1; i<(3*NrOfFirstZoneCandidateElements); i++){
            XCoordOfTempEdgeNode1 = NodesCoordX[TempExtraEdges1[i]];
            YCoordOfTempEdgeNode1 = NodesCoordY[TempExtraEdges1[i]];
            XCoordOfTempEdgeNode2 = NodesCoordX[TempExtraEdges2[i]];
            YCoordOfTempEdgeNode2 = NodesCoordY[TempExtraEdges2[i]];
            XTempNodeCoord = (XCoordOfTempEdgeNode1 + XCoordOfTempEdgeNode2)/2;
            YTempNodeCoord = (YCoordOfTempEdgeNode1 + YCoordOfTempEdgeNode2)/2;
            TempEdgeLength = SQRT((XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)*(XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)
            +(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1)*(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1));
            if (TempEdgeLength > Tol){
                // Check if new point is very close to already added new point
                AddCoordFlag = 1;
                for (j=0; j<NrOfNewNodes; j++){
                    DistNodes = SQRT((NewNodesCoordX[j]-XTempNodeCoord)*(NewNodesCoordX[j]-XTempNodeCoord)
                    + (NewNodesCoordY[j]-YTempNodeCoord)*(NewNodesCoordY[j]-YTempNodeCoord));
                    if (DistNodes < 2*EPSILON){
                        AddCoordFlag = 0;
                        break;
                    }
                }
                // If not then add new point
                if (AddCoordFlag == 1){
                    NewNodesCoordX[NrOfNewNodes] = XTempNodeCoord;
                    NewNodesCoordY[NrOfNewNodes] = YTempNodeCoord;
                    NrOfNewNodes = NrOfNewNodes+1;
                }
            }
        }
        
        // removing the first new node if the edge is too short
        XCoordOfTempEdgeNode1 = NodesCoordX[TempExtraEdges1[0]];
        YCoordOfTempEdgeNode1 = NodesCoordY[TempExtraEdges1[0]];
        XCoordOfTempEdgeNode2 = NodesCoordX[TempExtraEdges2[0]];
        YCoordOfTempEdgeNode2 = NodesCoordY[TempExtraEdges2[0]];
        TempEdgeLength = SQRT((XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)*(XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)
        +(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1)*(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1));
        if (TempEdgeLength<Tol){
            NrOfNewNodes = NrOfNewNodes-1;
            for (i=0; i<NrOfNewNodes; i++){
                NewNodesCoordX[i] = NewNodesCoordX[i+1];
                NewNodesCoordY[i] = NewNodesCoordY[i+1];
            }
        }
//         printf("NrOfNewNodes=%ld\n", NrOfNewNodes);
//     return;
        
        //  For triangles that showed up exactly once check if they are skinny        
        REAL TempLengths[3] = {0};
        REAL XNode1Coord, XNode2Coord, XNode3Coord;
        REAL YNode1Coord, YNode2Coord, YNode3Coord;
        REAL SkinnyTriangle = 3.0;
        for (i=0; i<NrOfElements; i++){
            if (FlagOfSecondZoneCandidateElements[i] == 1){
                XNode1Coord = NodesCoordX[Elements1[i]];
                YNode1Coord = NodesCoordY[Elements1[i]];
                XNode2Coord = NodesCoordX[Elements2[i]];
                YNode2Coord = NodesCoordY[Elements2[i]];
                XNode3Coord = NodesCoordX[Elements3[i]];
                YNode3Coord = NodesCoordY[Elements3[i]];
                
                TempLengths[0] = SQRT((XNode2Coord-XNode1Coord)*(XNode2Coord-XNode1Coord)
                + (YNode2Coord- YNode1Coord)*(YNode2Coord- YNode1Coord));
                TempLengths[1] = SQRT((XNode3Coord-XNode2Coord)*(XNode3Coord-XNode2Coord)
                + (YNode3Coord- YNode2Coord)*(YNode3Coord- YNode2Coord));
                TempLengths[2] = SQRT((XNode1Coord-XNode3Coord)*(XNode1Coord-XNode3Coord)
                + (YNode1Coord- YNode3Coord)*(YNode1Coord- YNode3Coord));
                // If skinny then add midpoint of triangle as newnode
                if (misc_max(3,&TempLengths[0])/misc_min(3,&TempLengths[0]) > SkinnyTriangle){
                    NewNodesCoordX[NrOfNewNodes] = (XNode1Coord+XNode2Coord+XNode3Coord)/3;
                    NewNodesCoordY[NrOfNewNodes] = (YNode1Coord+YNode2Coord+YNode3Coord)/3;
                    NrOfNewNodes = NrOfNewNodes+1;
                }
            }
        }
//         printf("NrOfNewNodes=%ld\n", NrOfNewNodes);

    }
//     printf("NrOfCandidateEdges=%ld\n",NrOfCandidateEdges );
    
    
// Evaluation of contour edges from all candidates edges
    NrOfCandidateElements = 0;
    ArrayOfCandidateElements= realloc(ArrayOfCandidateElements, 3*NrOfElements*sizeof(UINT));
    if (ArrayOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    ret_code = delaunay_triangulation_edgeAttachments(NrOfElements,Elements1, Elements2, Elements3,
            NrOfCandidateEdges, CandidateEdges1, CandidateEdges2, &NrOfCandidateElements, ArrayOfCandidateElements);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
//    printf("NrOfCandidateElements=%ld\n",NrOfCandidateElements);
// Removing repeated mention of triangles
    ElementsCount = realloc(ElementsCount, NrOfElements*sizeof(UINT));
    if (ElementsCount == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i=0; i<NrOfElements; i++)
        ElementsCount[i] = 0;
    
    for (i=0; i<NrOfCandidateElements; i++)
        ElementsCount[ArrayOfCandidateElements[i]] = 1;
    
    IDOfCandidateElements = malloc(NrOfCandidateElements*sizeof(UINT));
    if (IDOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    NrOfCandidateElements = 0; // without reptition
    for (i=0; i<NrOfElements; i++){
        if (ElementsCount[i] == 1){
            IDOfCandidateElements[NrOfCandidateElements] = i;
            NrOfCandidateElements = NrOfCandidateElements+1;
            
        }
    }
//     printf("NrOfCandidateElements=%ld\n",NrOfCandidateElements);
//
// % Collecting all the edges of all candidate triangles
// TempEdges1=zeros(NoOfCandidateElements*3,1);
// TempEdges2=zeros(NoOfCandidateElements*3,1);
    TempExtraEdges1 = realloc(TempExtraEdges1, NrOfCandidateElements*3*sizeof(UINT));
    TempExtraEdges2 = realloc(TempExtraEdges2, NrOfCandidateElements*3*sizeof(UINT));
    if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i=0; i<NrOfCandidateElements; i++){
        TempExtraEdges1[i*3] = Elements1[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3] = Elements2[IDOfCandidateElements[i]];
        TempExtraEdges1[i*3 + 1] = Elements2[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3 + 1] = Elements3[IDOfCandidateElements[i]];
        TempExtraEdges1[i*3 + 2] = Elements3[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3 + 2] = Elements1[IDOfCandidateElements[i]];
    }
    
    // Removing edges to keep only outer contours
    
    UINT NrOfContourEdges = 0;
    MultiplicationOfTempEdges = calloc(3*NrOfCandidateElements, sizeof(UINT)); //NOTE calloc
    ContourEdges1 = malloc(3*NrOfCandidateElements*sizeof(UINT));
    ContourEdges2 = malloc(3*NrOfCandidateElements*sizeof(UINT));
    if (MultiplicationOfTempEdges == NULL || ContourEdges1 == NULL ||
            ContourEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    INT NrOfEdge;
    // Possibly simplify this
    for (i=0; i<(3*NrOfCandidateElements); i++){
        if (MultiplicationOfTempEdges[i] == 0){
            NrOfEdge = -1;
            for (j=0; j<(3*NrOfCandidateElements); j++){
                if (j != i){
                    if ((TempExtraEdges2[j] == TempExtraEdges1[i] && TempExtraEdges1[j] == TempExtraEdges2[i])
                            || (TempExtraEdges1[j] == TempExtraEdges1[i] && TempExtraEdges2[j] == TempExtraEdges2[i])){
                        NrOfEdge = j;
                        MultiplicationOfTempEdges[NrOfEdge] = 2;
                    }
                }
            }
            if (NrOfEdge == -1){
                MultiplicationOfTempEdges[i] = 1;
                ContourEdges1[NrOfContourEdges] = TempExtraEdges1[i];
                ContourEdges2[NrOfContourEdges] = TempExtraEdges2[i];
                NrOfContourEdges = NrOfContourEdges+1;
            }else
                MultiplicationOfTempEdges[i] = 2;
        }
    }
    
//     printf("NrOfContourEdges=%ld\n",NrOfContourEdges);
    
    // evaluation of the regions
    NrOfNodesOfRegions = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    Regions = calloc(5*NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    IndexOfNextEdge = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    TempNodes = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    SideFlag = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    if (NrOfNodesOfRegions == NULL || Regions == NULL || IndexOfNextEdge == NULL
            || TempNodes == NULL || SideFlag == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    
    UINT NodesCount = 1;
    UINT NrOfRegions = 1;
    Regions[NodesCount-1] = ContourEdges1[0];
    NrOfNodesOfRegions[NrOfRegions-1]++;
    UINT RefNode = ContourEdges2[0];
    UINT ContourEdgesStart = 1;
    UINT NrOfIndexOfNextEdge;
    UINT PrevNode, Index;
    while (ContourEdgesStart < NrOfContourEdges){
        NrOfIndexOfNextEdge = 0;
        for (j=ContourEdgesStart; j<NrOfContourEdges; j++){
            if (ContourEdges1[j] == RefNode){
                IndexOfNextEdge[NrOfIndexOfNextEdge] = j;
                SideFlag[NrOfIndexOfNextEdge] = 1;
                NrOfIndexOfNextEdge = NrOfIndexOfNextEdge+1;
            }else if (ContourEdges2[j] == RefNode){
                IndexOfNextEdge[NrOfIndexOfNextEdge] = j;
                SideFlag[NrOfIndexOfNextEdge] = 2;
                NrOfIndexOfNextEdge = NrOfIndexOfNextEdge+1;
            }
        }
        if (NrOfIndexOfNextEdge == 0){
            NodesCount++;
            Regions[NodesCount-1] = RefNode;
            NrOfNodesOfRegions[NrOfRegions-1]++;
            
            if (ContourEdgesStart< NrOfContourEdges){
                NrOfRegions++;
                NodesCount++;
                Regions[NodesCount-1] = ContourEdges1[ContourEdgesStart];
                NrOfNodesOfRegions[NrOfRegions-1]++;
                RefNode = ContourEdges2[ContourEdgesStart];
                ContourEdgesStart = ContourEdgesStart+1;
            }
        }else{
            if (NrOfIndexOfNextEdge > 1){
                PrevNode = Regions[NodesCount-1];
                for (i=0; i<NrOfIndexOfNextEdge; i++){
                    if (SideFlag[i] == 1)
                        TempNodes[i] = ContourEdges2[IndexOfNextEdge[NrOfIndexOfNextEdge]];
                    else if (SideFlag[i] == 2)
                        TempNodes[i] = ContourEdges1[IndexOfNextEdge[NrOfIndexOfNextEdge]];
                }
                Index = FindNextNode(NodesCoordX, NodesCoordY, PrevNode, RefNode, TempNodes, NrOfIndexOfNextEdge);
                IndexOfNextEdge[0] = IndexOfNextEdge[Index];
                SideFlag[0] = SideFlag[Index];
                NrOfIndexOfNextEdge = 1;
            }
            NodesCount++;
            Regions[NodesCount-1] = RefNode;
            NrOfNodesOfRegions[NrOfRegions-1]++;
            if (SideFlag[0] == 1)
                RefNode = ContourEdges2[IndexOfNextEdge[0]];
            else if (SideFlag[0] == 2)
                RefNode = ContourEdges1[IndexOfNextEdge[0]];
            
            NrOfContourEdges = NrOfContourEdges - 1;
            for (j=ContourEdgesStart; j<NrOfContourEdges; j++){
                if (j>=IndexOfNextEdge[0]){
                    ContourEdges1[j] = ContourEdges1[j+1];
                    ContourEdges2[j] = ContourEdges2[j+1];
                }
            }
        }
    }
    NodesCount++;
    Regions[NodesCount-1] = RefNode;
    NrOfNodesOfRegions[NrOfRegions-1]++;
    printf("NrOfRegions=%ld\n",NrOfRegions);
//     for (i=0; i<NrOfRegions; i++)
//         printf("Nodes in region %ld = %ld\n",i, NrOfNodesOfRegions[i]);
    
    
    
    z_m = calloc(NrOfRegions, sizeof(INT)); //NOTE calloc
    z = calloc(NrOfRegions, sizeof(COMPLEX)); //NOTE calloc
    if (z_m == NULL || z == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    NodesCount = 0;
    INT dQ;
    for (i=0; i<NrOfRegions; i++){
        for (j=1; j<NrOfNodesOfRegions[i]; j++){
            dQ = Quadrants[Regions[NodesCount+j]] - Quadrants[Regions[NodesCount+j-1]];
            if (dQ == 3)
                z_m[i] = z_m[i] - 1;
            else if (dQ == -3)
                z_m[i] = z_m[i] + 1;
            else if (ABS(dQ) == 2)
                z_m[i] = z_m[i] + NAN;
            else
                z_m[i] = z_m[i] + dQ;
            
        }
        z_m[i] = ABS(z_m[i])/4; //taking abs because orientation of contour not known
        
        z[i] = 0;
        for (j=0; j<NrOfNodesOfRegions[i]; j++)
            z[i] = z[i]+(NodesCoordX[Regions[NodesCount+j]] + I*NodesCoordY[Regions[NodesCount+j]]);

        z[i] = z[i]/NrOfNodesOfRegions[i];
        
        NodesCount = NodesCount+NrOfNodesOfRegions[i];
    }
    
    j = 0;
    for (i=0; i<NrOfRegions; i++){
        if (z_m[i]>0 && z_m[i] != NAN){
         bound_states[j] = z[i];
         j++;
        }
    }
    *K_ptr = j;
//     misc_print_buf(*K_ptr,bound_states,"bs");
   
    leave_fun:
        free(NewNodesCoordX);
        free(NewNodesCoordY);
        free(NodesCoordX);
        free(NodesCoordY);
        free(FuntionValues);
        free(ComplexNodes);
        free(Quadrants);
        free(Angles);
        free(ValidEdges);
        free(ValidEdges1);
        free(ValidEdges2);
        free(CandidateEdges1);
        free(CandidateEdges2);
        free(Elements1);
        free(Elements2);
        free(Elements3);
        free(FlagOfFirstZoneCandidateElements);
        free(FlagOfSecondZoneCandidateElements);
        free(ElementsCount);
        free(IDOfCandidateElements);
        free(TempExtraEdges1);
        free(TempExtraEdges2);
        free(ContourEdges1);
        free(ContourEdges2);
        free(MultiplicationOfTempEdges);
        free(SideFlag);
        free(TempNodes);
        free(IndexOfNextEdge);
        free(Regions);
        free(NrOfNodesOfRegions);
        free(z);
        free(z_m);
        free(ArrayOfCandidateElements);
        free(CandidateNodes);
        free(Temp);
        if (DT_built == 1){
            free(DT_ptr->NodeOfTriangles1);
            free(DT_ptr->NodeOfTriangles2);
            free(DT_ptr->NodeOfTriangles3);
            free(DT_ptr->EdgesOfTriangles1);
            free(DT_ptr->EdgesOfTriangles2);
            free(DT_ptr->EdgesOfTriangles3);
            free(DT_ptr->Edges1);
            free(DT_ptr->Edges2);
            free(DT_ptr->NodesCoordX);
            free(DT_ptr->NodesCoordY);
            free(DT_ptr->XCoordOfCCOfTriangles);
            free(DT_ptr->YCoordOfCCOfTriangles);
            free(DT_ptr->RadiusOfCCOfTriangles);
            free(DT_ptr->StatusOfTriangles);
            free(DT_ptr->CheckEdge);
        }
        return ret_code;

}

static UINT FindNextNode(REAL *NodesCoordX, REAL *NodesCoordY, UINT PrevNode, 
        UINT RefNode, UINT *TempNodes, UINT NrOfTempNodes)
{
// %  FindNextNode: finds the next node in the candidate region boudary
// %                process. The next one (after the reference one) is picked
// %                from the fixed set of nodes.
// %
// % INPUTS
// %
// %  NodesCoord     : nodes coordinates
// %  PrevNode       : previous node
// %  RefNode        : reference (current) node
// %  TempNodes      : set of nodes
// %
// % OUTPUTS
// %
// %  Index          : index of the next node
// %

    UINT Index = 0, j;
    REAL min_Phi = 10;
    REAL PrevNodeCoordX = NodesCoordX[PrevNode];
    REAL PrevNodeCoordY = NodesCoordY[PrevNode];
    REAL RefNodeCoordX = NodesCoordX[RefNode];
    REAL RefNodeCoordY = NodesCoordY[RefNode];
    REAL PRX, PRY, TRX, TRY, LenPR, LenTR, DotProd, Phi;
    for (j=0; j<NrOfTempNodes; j++){
        PRX = PrevNodeCoordX-RefNodeCoordX;
        PRY = PrevNodeCoordY-RefNodeCoordY;
        TRX = NodesCoordX[TempNodes[j]] - RefNodeCoordX;
        TRY = NodesCoordY[TempNodes[j]] - RefNodeCoordY;
        LenPR = SQRT(PRX*PRX + PRY*PRY);
        LenTR = SQRT(TRX*TRX + TRY*TRY);
        DotProd = PRX*TRX + PRY*TRY;
        Phi = ACOS(DotProd/(LenPR*LenTR));
        if ((PRX*TRY-PRY*TRX)<0)
            Phi = 2*PI - Phi;
        if (Phi < min_Phi){
            min_Phi = Phi;
            Index = j;
        }
    }
    
    return Index;
}
