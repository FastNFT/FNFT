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
 * Sander Wahls (TU Delft) 2017-2018, 2023.
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Marius Brehler (TU Dortmund) 2018.
 * Peter J Prins (TU Delft) 2020-2021.
 * Sander Wahls (KIT) 2023.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_kdvv.h"

static fnft_kdvv_opts_t default_opts = {
    .bound_state_localization = kdvv_bsloc_GRIDSEARCH_AND_REFINE,
    .niter = 10,
    .discspec_type = kdvv_dstype_NORMING_CONSTANTS,
    .contspec_type = kdvv_cstype_REFLECTION_COEFFICIENT,
    .normalization_flag = 1,
    .discretization = kdv_discretization_2SPLIT4B,
    .richardson_extrapolation_flag = 0,
    .grid_spacing = -1
};

/**
 * Creates a new options variable for fnft_kdvv with default settings.
 * See the header file for a detailed description.
 */
fnft_kdvv_opts_t fnft_kdvv_default_opts()
{
    return default_opts;
}

/**
 * Declare auxiliary routines used by the main routine fnft_kdvv.
 * Their bodies follow below.
 */

static inline INT kdvv_compute_boundstates(
        UINT const D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const REAL eps_t,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        fnft_kdvv_opts_t const * const opts);

static inline INT fnft_kdvv_base(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT M,
        COMPLEX * const contspec,
        REAL const * const XI,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        const INT kappa,
        fnft_kdvv_opts_t *opts);

static inline INT kdvv_compute_contspec(
        const UINT deg,
        const INT W,
        COMPLEX * const transfer_matrix,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT D,
        REAL const * const XI,
        const UINT M,
        COMPLEX * const result,
        const INT kappa,
        fnft_kdvv_opts_t const * const opts);

static inline INT kdvv_compute_normconsts_or_residues(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT K,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        fnft_kdvv_opts_t const * const opts);

static inline INT kdvv_refine_bound_states_newton(const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        UINT K,
        COMPLEX * const bound_states,
        kdv_discretization_t const discretization,
        UINT const niter,
        REAL const * const bounding_box,
        const INT normalization_flag);

/**
 * Fast nonlinear Fourier transform for the nonlinear Schroedinger
 * equation with vanishing boundary conditions.
 * This function takes care of the necessary preprocessing (for example
 * signal resampling or subsampling) based on the options before calling
 * fnft_kdvv_base. If the richardson_extrapolation_flag is set, this function
 * calls fnft_kdvv_base with half of the samples and then performs
 * Richardson extrapolation on the spectrum.
 */
INT fnft_kdvv(
        const UINT D,
        COMPLEX const * const q,
        REAL const * const T,
        const UINT M,
        COMPLEX * const contspec,
        REAL const * const XI,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        fnft_kdvv_opts_t *opts)
{
    const INT kappa = 1;
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
    fnft_kdvv_bsloc_t bs_loc_opt = 0;
    fnft_kdvv_dstype_t ds_type_opt = 0;
    INT ret_code = SUCCESS;
    UINT i, j, upsampling_factor, D_effective;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    for (i=0; i<D; i++) {
        if (CIMAG(q[i])!=0.0)
            return E_INVALID_ARGUMENT(q);
    }
    if (T == NULL || T[0] >= T[1])
        return E_INVALID_ARGUMENT(T);
    if (contspec != NULL) {
        if (XI == NULL || XI[0] >= XI[1])
            return E_INVALID_ARGUMENT(XI);
    }
    for (i=0; i<M; i++){
        if (XI[0]+i*(XI[1]-XI[0])/(M-1) == 0)
            return E_INVALID_ARGUMENT(XI);
    }

    if (opts == NULL)
        opts = &default_opts;

    if (bound_states != NULL) {
        if (K_ptr == NULL) {
            return E_INVALID_ARGUMENT(K_ptr);
        } else if (opts->bound_state_localization==fnft_kdvv_bsloc_NEWTON) {
            // This method needs user-provided initial guesses. Check those:
            for (i=0; i<*K_ptr; i++) {
                if ( CREAL(bound_states[i])!=0.0 || CIMAG(bound_states[i])<=0.0 )
                    return E_INVALID_ARGUMENT(bound_states);
            }
        }
    }

    // Some higher-order discretizations require samples on a non-equidistant grid
    // while others require derivatives which are computed in the form of
    // finite-differences. The input array q has D samples corresponding to
    // an equidistant grid. The upsampling_factor*D gives the effective
    // number of samples for the chosen discretization. The effective number
    // of samples is required for the auxiliary function calls.
    upsampling_factor = kdv_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }

    D_effective = D * upsampling_factor;

    // Determine step size
    const REAL eps_t = (T[1] - T[0])/(D - 1);

    // If Richardson extrapolation is not requested or if it is requested but
    // opts->discspec_type is kdvv_dstype_NORMCONSTS,then normconsts_or_residues_reserve
    // acts as just a place holder for normconsts_or_residues. This is equivalent
    // to calling the auxiliary functions with normconsts_or_residues instead of
    // normconsts_or_residues_reserve.
    normconsts_or_residues_reserve = normconsts_or_residues;
    // If Richardson extrapolation is requested and opts->discspec_type
    // is kdvv_dstype_RESIDUES, *K_ptr*2  memory is allocated
    // for normconsts_or_residues_reserve. This overwrites the previous
    // assignment of normconsts_or_residues to normconsts_or_residues_reserve.
    // Richardson extrpolation is applied on normconsts(b) and aprimes separately
    // before combining the results to get the residues=normconsts/aprimes.
    // Hence double the memory is necessary even if the user requests only residues.
    if (opts->richardson_extrapolation_flag == 1){
        ds_type_opt = opts->discspec_type;
        if (ds_type_opt == kdvv_dstype_RESIDUES){
            opts->discspec_type = kdvv_dstype_BOTH;
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
    // the required derivatives using finite differences. The kdv_discretization_preprocess_signal
    // function also performs some further discretization specific processing.
    // Preprocessing takes care of computing things which are required by all
    // the auxiliary functions thus helping efficiency.
    Dsub = D;
    ret_code = kdv_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
            first_last_index, opts->discretization);
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = fnft_kdvv_base(D_effective, q_preprocessed, r_preprocessed, T, M, contspec, XI, K_ptr,
                bound_states, normconsts_or_residues, kappa, opts);
    CHECK_RETCODE(ret_code, leave_fun);

    if (opts->richardson_extrapolation_flag == 1){
        // Allocating memory
        UINT contspec_len = 0;
        if (contspec != NULL && M > 0){
            switch (opts->contspec_type) {
                case kdvv_cstype_BOTH:
                    contspec_len = 3*M;
                    break;
                case kdvv_cstype_REFLECTION_COEFFICIENT:
                    contspec_len = M;
                    break;
                case kdvv_cstype_AB:
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
                case kdvv_dstype_BOTH:
                case kdvv_dstype_RESIDUES:
                    discspec_len = 2*K_sub;
                    break;
                case kdvv_dstype_NORMING_CONSTANTS:
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
        method_order = kdv_discretization_method_order(opts->discretization);
        if (method_order == 0){
            ret_code =  E_INVALID_ARGUMENT(discretization);
            goto leave_fun;
        }

        // The signal q is now subsampled(approx. half the samples) and
        // preprocessed as required for the discretization. This is
        // required for obtaining a second approximation of the spectrum
        // which will be used for Richardson extrapolation.
        Dsub = (UINT) CEIL(D/2);
        ret_code = kdv_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
                first_last_index, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        Tsub[0] = T[0] + first_last_index[0]*eps_t;
        Tsub[1] = T[0] + first_last_index[1]*eps_t;
        const REAL eps_t_sub = (Tsub[1] - Tsub[0])/(Dsub - 1);

        // Calling fnft_kdvv_base with subsampled signal
        bs_loc_opt = opts->bound_state_localization;
        opts->bound_state_localization = kdvv_bsloc_NEWTON;

        ret_code = fnft_kdvv_base(Dsub * upsampling_factor, qsub_preprocessed, rsub_preprocessed, Tsub, M, contspec_sub, XI, &K_sub,
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
                    if (ds_type_opt == kdvv_dstype_RESIDUES || ds_type_opt == kdvv_dstype_BOTH){
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
            if (ds_type_opt == kdvv_dstype_RESIDUES)
                memcpy(normconsts_or_residues,normconsts_or_residues_reserve+K,K* sizeof(COMPLEX));
            else if(ds_type_opt == kdvv_dstype_BOTH)
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
        if (normconsts_or_residues_reserve != normconsts_or_residues)
            free(normconsts_or_residues_reserve);
        return ret_code;
}

// Auxiliary function: Base routine for fnft_kdvv. fnft_kdvv preprocesses the signals
// and calls this function with different options as needed. This prevents
// code doubling while being efficient.
static inline INT fnft_kdvv_base(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT M,
        COMPLEX * const contspec,
        REAL const * const XI,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        const INT kappa,
        fnft_kdvv_opts_t *opts)
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
    upsampling_factor = kdv_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }
    D_given = D/upsampling_factor;
    const REAL eps_t = (T[1] - T[0])/(D_given - 1);

    // D should be the effective number of samples in q
    i = kdv_fscatter_numel(D, opts->discretization);
    // NOTE: At this stage if i == 0 it means the discretization corresponds
    // to a slow method. Incorrect discretizations will have been checked for
    // in fnft_kdvv main

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
            
        ret_code = kdv_fscatter(D, q, r, eps_t, kappa, transfer_matrix, &deg, W_ptr,
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
        ret_code = kdvv_compute_contspec(deg, W, transfer_matrix, q, r, T, D, XI, M,
                contspec, kappa, opts);
        CHECK_RETCODE(ret_code, leave_fun);
    }

    // Compute the discrete spectrum
    if (bound_states != NULL && *K_ptr!=0) {
        // Compute the bound states
        ret_code = kdvv_compute_boundstates(D, q, r, T,
                eps_t, K_ptr, bound_states, opts);
        CHECK_RETCODE(ret_code, leave_fun);

        // Norming constants and/or residues)
        if (normconsts_or_residues != NULL && *K_ptr != 0) {
            ret_code = kdvv_compute_normconsts_or_residues(D, q, r, T, *K_ptr,
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

// Auxiliary function: Computes the bound states.
static inline INT kdvv_compute_boundstates(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const REAL eps_t,
        UINT * const K_ptr,
        COMPLEX * const bound_states,
        fnft_kdvv_opts_t const * const opts)
{
    INT ret_code = SUCCESS;
    UINT K=*K_ptr;
    UINT upsampling_factor, i;
    COMPLEX * xi = NULL;
    COMPLEX * scatter_coeffs = NULL;
    INT * log_scaling_factors = NULL;
    INT const kappa = 1;

    upsampling_factor = kdv_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }

    // Copy the values of opts and modify it from fast to slow discretization if needed.
    fnft_kdvv_opts_t opts_slow = *opts;
    ret_code = fnft__kdv_slow_discretization(&(opts_slow.discretization));
    CHECK_RETCODE(ret_code, leave_fun);

    // Calculate the bounding box for the bound states. Therefore we need to find the maximum of the real part of q(t), skipping t-derivative samples
    REAL bound_squared = 0;

    if (upsampling_factor<=2) {
        // For the currently implemented discretizations with an upsampling factor
        // of 1 or 2, all substep sizes are real. Therefore the spectrum we calculate
        // is the exact spectrum of a staircase reconstruction from the preprocessed
        // samples. Hence the maximum of the potential samples gives an upperbound
        // on the eigenvalues.
        for (i = 0; i < D; i++)
            bound_squared = bound_squared > (REAL)CREAL(q[i]) ? bound_squared : (REAL)CREAL(q[i]);
    } else {
        // For other discretizations the above bound can be violated
        // due to interpolation of the potential. Therefore we apply a safety factor.
        switch (opts->discretization) {
            case kdv_discretization_ES4_VANILLA:
            case kdv_discretization_TES4_VANILLA:
            case kdv_discretization_ES4:
            case kdv_discretization_TES4:
                // These discretizations contain samples of the first and
                // second derivative in the preprocessed potential. We need
                // to skip those.
                for (i = 0; i < D; i+=3)
                    bound_squared = bound_squared > (REAL)CREAL(q[i]) ? bound_squared : (REAL)CREAL(q[i]);
                // Apply a safety factor, because the effective maximum may be higher due to interpolation
                break;
            default:
                // Other discretizations with an order >2 contain weighted averages
                // of interpolated potential samples in the preprocessed q.
                // A theoretical calculation of the eigenvalue bound is
                // complicated and not worth the effort. Instead, we estimate
                // the bound from the maximum of the real part of the
                // preprocessed q and proceed with fingers crossed.
                for (i = 0; i < D; i++)
                    bound_squared = bound_squared > (REAL)CREAL(q[i]) ? bound_squared : (REAL)CREAL(q[i]);
                break;
        }
        bound_squared *= 2.0; // Safety factor
    }

    if (bound_squared==0.0){
        // Non-positive potentials have no bound states
        *K_ptr=0;
        goto leave_fun;
    }

    REAL const bounding_box[4] = {0,0,0,SQRT(bound_squared)};

    // Localize bound states ...
    switch (opts->bound_state_localization) {
        case kdvv_bsloc_GRIDSEARCH_AND_REFINE:
        {
            // Set the search grid
            const REAL dist = bounding_box[3] - bounding_box[2];
            if (dist <= 0) {
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            UINT M = 0;
            if (opts->grid_spacing == -1) {
                M = 10000;
            } else if (opts->grid_spacing > 0) {
                M = CEIL(dist/opts->grid_spacing);
                if (M<2) {
                    *K_ptr = 0;
                    goto leave_fun;
                }   
            } else {
                ret_code = E_INVALID_ARGUMENT(grid_spacing);
                goto leave_fun;
            }
            COMPLEX const eps_xi = I*dist/(M - 1);
            xi = malloc(M * sizeof(COMPLEX));
            if (xi == NULL) {
                ret_code = E_NOMEM;
                goto leave_fun;
            }

            xi[0] = I*nexttoward(bounding_box[2],bounding_box[3]);
            for (i = 1; i < M-1; i++)
                xi[i] = I * CIMAG(eps_xi*i);
            xi[M-1] = I*nexttoward(bounding_box[3],bounding_box[2]);

            // Allocate memory for call to kdv_scatter_matrix
            scatter_coeffs = malloc(4 * M * sizeof(COMPLEX));
            CHECK_NOMEM(scatter_coeffs, ret_code, leave_fun);

            if (opts->normalization_flag) {
                log_scaling_factors = malloc(M * sizeof(COMPLEX));
                CHECK_NOMEM(log_scaling_factors, ret_code, leave_fun);
            }

            // Compute a(xi) on the grid on the imaginary axis
            UINT const derivative_flag = 0;
            ret_code = kdv_scatter_matrix(D, q, r, eps_t, kappa,
                                          M, xi, scatter_coeffs, 
                                          log_scaling_factors,
                                          opts_slow.discretization,
                                          derivative_flag);
            CHECK_RETCODE(ret_code, leave_fun);

            // Search for zerocrossings of a(xi)
            K=0;
            for (i = 0; i < M; i++) {
                if (scatter_coeffs[4*i]==0) {
                    if (K>=*K_ptr) { // Lucky us: Exact bound state found
                        ret_code = E_OTHER("More than *K_ptr initial guesses for bound states found. Increase *K_ptr and try again.")
                        CHECK_RETCODE(ret_code, leave_fun);
                    }
                    bound_states[K]=xi[i];
                    K++;
                } else if ( i>0 && CREAL(scatter_coeffs[4*i])*CREAL(scatter_coeffs[4*(i-1)])<0) {
                    if (K>=*K_ptr) { // Bound state between this sample and the previous
                        ret_code = E_OTHER("More than *K_ptr initial guesses for bound states found. Increase *K_ptr and try again.")
                        CHECK_RETCODE(ret_code, leave_fun);
                    }
                    // Zero crossing estimate with regula falsi:
                    REAL scl1 = 1, scl2 = 1;
                    if (opts->normalization_flag) {
                        scl1 = POW(2, log_scaling_factors[i]);
                        scl2 = POW(2, log_scaling_factors[i-1]);
                    }
                    bound_states[K] = (xi[i-1] * scl1*CREAL(scatter_coeffs[4*i]) - xi[i] * scl2*CREAL(scatter_coeffs[4*(i-1)])) / CREAL( scl1*scatter_coeffs[4*i] - scl2*scatter_coeffs[4*(i-1)]);
                    K++;
                }
            }
        }
            // fall through
        case kdvv_bsloc_NEWTON: // Refine using Newton's method

            // Perform Newton iterations. Initial guesses of bound-states
            // should be in the continuous-time domain.

            ret_code = kdvv_refine_bound_states_newton(D, q, r, T, K, bound_states,
                    opts_slow.discretization, opts_slow.niter, bounding_box, opts_slow.normalization_flag);
            CHECK_RETCODE(ret_code, leave_fun);

            break;

        default:

            return E_INVALID_ARGUMENT(opts->bound_state_localization);
    }

    // Update number of bound states
    *K_ptr = K;

leave_fun:
    free(xi);
    free(scatter_coeffs);
    free(log_scaling_factors);
    return ret_code;
}

// Auxiliary function: Computes continuous spectrum on a frequency grid
static inline INT kdvv_compute_contspec(
        const UINT deg,
        const INT W,
        COMPLEX * const transfer_matrix,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT D,
        REAL const * const XI,
        const UINT M,
        COMPLEX * const result,
        const INT kappa,
        fnft_kdvv_opts_t const * const opts)
{
    COMPLEX * H11_vals = NULL;
    COMPLEX A, V;
    REAL scale;
    REAL phase_factor_rho, phase_factor_a, phase_factor_b;
    INT ret_code = SUCCESS;
    UINT i, offset = 0, upsampling_factor, D_given;
    COMPLEX * scatter_coeffs = NULL, * xi = NULL;
    INT * log_scaling_factors = NULL;

    // Determine step size
    // D is interpolated number of samples but eps_t is the step-size
    // corresponding to original number of samples.
    upsampling_factor = kdv_discretization_upsampling_factor(opts->discretization);
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
    H11_vals = malloc(4*M * sizeof(COMPLEX));
    if (H11_vals == NULL){
        return E_NOMEM;
        goto leave_fun;}
    COMPLEX * const H12_vals = H11_vals + M;
    COMPLEX * const H21_vals = H12_vals + M;
    COMPLEX * const H22_vals = H21_vals + M;

    // If the discretization is a slow method then there should be no transfer_matrix
    if (deg == 0 && transfer_matrix == NULL && W == 0){

        // Allocate memory for call to kdv_scatter_matrix
        scatter_coeffs = malloc(4 * M * sizeof(COMPLEX));
        CHECK_NOMEM(scatter_coeffs, ret_code, leave_fun);

        if (opts->normalization_flag) {
            log_scaling_factors = malloc(M * sizeof(COMPLEX));
            CHECK_NOMEM(log_scaling_factors, ret_code, leave_fun);
        }

        ret_code = kdv_scatter_matrix(D, q, r, eps_t, kappa, M,
                xi, scatter_coeffs, log_scaling_factors,
                opts->discretization, 0);
        CHECK_RETCODE(ret_code, leave_fun);

        // This is necessary because kdv_scatter_matrix to ensure
        // boundary conditions can be applied using common code for slow
        // methods and polynomial transfer matrix based methods.
        for (i = 0; i < M; i++){
            REAL scl = 1;
            if (opts->normalization_flag)
                scl = POW(2, log_scaling_factors[i]);
            H11_vals[i] = scl*scatter_coeffs[i*4];
//            H12_vals[i] = scatter_coeffs[i*4+1]; // Not used
            H21_vals[i] = scl*scatter_coeffs[i*4+2];
//            H22_vals[i] = scatter_coeffs[i*4+3]; // Not used
        }

    }else{
        // Prepare the use of the chirp transform. The entries of the transfer
        // matrix that correspond to a and b will be evaluated on the frequency
        // grid xi(i) = XI1 + i*eps_xi, where i=0,...,M-1. Since
        // z=exp(2.0*I*XI*eps_t/degree1step), we find that the z at which z the transfer
        // matrix has to be evaluated are given by z(i) = 1/(A * V^-i), where:
        V = eps_xi;
        ret_code = kdv_discretization_lambda_to_z(1, eps_t, &V, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);
        A = -XI[0];
        ret_code = kdv_discretization_lambda_to_z(1, eps_t, &A, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, H11_vals);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix+1*(deg+1), A, V, M, H12_vals);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix+2*(deg+1), A, V, M, H21_vals);
        CHECK_RETCODE(ret_code, leave_fun);

        ret_code = poly_chirpz(deg, transfer_matrix+3*(deg+1), A, V, M, H22_vals);
        CHECK_RETCODE(ret_code, leave_fun);

        // NOTE: At this point the matrix [Hnm_vals] is NOT the change of state matrix in AKNS basis. All its values differ by a scalar factor exp(-I*xi*D) due to dividing the Laurent polynomials by the most negative power to obtain ordinary polynomials. This will finally be corrected by the phase_factor_a and phase_factor_b below.
        // Change the basis of the change of state matrix to the S-basis.
        COMPLEX Tmx[2][2], Mmx[2][2], Hmx[2*2];
        for (i=0; i<M; i++) { // Loop over lambda samples
            // Copy the values of the change of state matrix in AKNS basis to a 2 dimensional array
            for (UINT n=0; n<4; n++){
                Hmx[n] = H11_vals[n*M+i]; // This exploints the adjacency of H11_vals ... H22_vals in memory.
            }

            // Fetch the change of basis matrix from S to the chosen discretization basis
            ret_code = kdv_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],xi[i],0,eps_t,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            // Right-multiply the change of state matrix by the change of basis matrix from S to the chosen discretization basis
            misc_matrix_mult(2,2,2,&Hmx[0],&Tmx[0][0],&Mmx[0][0]);

            // Fetch the change of basis matrix from the chosen discretization basis to S
            ret_code = kdv_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],xi[i],0,eps_t,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            // Left-multiply the change of state matrix by the change of basis matrix from the chosen discretization basis to S
            misc_matrix_mult(2,2,2,&Tmx[0][0],&Mmx[0][0],&Hmx[0]);

            // Copy to the result
            for (UINT n=0; n<4; n++){
                H11_vals[n*M+i] = Hmx[n]; // This exploints the adjacency of H11_vals ... H22_vals in memory.
            }
        }
    }
    // Compute the continuous spectrum
    switch (opts->contspec_type) {

        case kdvv_cstype_BOTH:

            offset = M;

        // fall through
        case kdvv_cstype_REFLECTION_COEFFICIENT:

            ret_code = kdv_discretization_phase_factor_rho(eps_t, T[1], &phase_factor_rho,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            for (i = 0; i < M; i++) {
                if (H11_vals[i] == 0.0){
                    return E_DIV_BY_ZERO;
                    goto leave_fun;
                }
                result[i] = H21_vals[i] * CEXP(I*xi[i]*phase_factor_rho) / H11_vals[i];
            }

            if (opts->contspec_type == kdvv_cstype_REFLECTION_COEFFICIENT)
                break;
            // fall through

        case kdvv_cstype_AB:

            scale = POW(2, W); // needed since the transfer matrix might
            // have been scaled by kdv_fscatter. W == 0 for slow methods.

            // Calculating the discretization specific phase factors.
            ret_code = kdv_discretization_phase_factor_a(eps_t, D_given, T, &phase_factor_a,opts->discretization);
            CHECK_RETCODE(ret_code, leave_fun);

            ret_code = kdv_discretization_phase_factor_b(eps_t, D_given, T, &phase_factor_b,opts->discretization);
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
        free(log_scaling_factors);
        return ret_code;
}

// Auxiliary function: Computes the norming constants and/or residues
// using slow scattering schemes
static inline INT kdvv_compute_normconsts_or_residues(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT K,
        COMPLEX * const bound_states,
        COMPLEX * const normconsts_or_residues,
        fnft_kdvv_opts_t const * const opts)
{
    COMPLEX *a_vals = NULL, *aprime_vals = NULL;
    UINT i, offset = 0;
    INT ret_code = SUCCESS;

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

    const UINT upsampling_factor = kdv_discretization_upsampling_factor(opts->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto leave_fun;
    }

    // Copy the values of opts and modify it from fast to slow discretization if needed.
    fnft_kdvv_opts_t opts_slow = *opts;
    ret_code=fnft__kdv_slow_discretization(&(opts_slow.discretization));
    CHECK_RETCODE(ret_code, leave_fun);

    ret_code = kdv_scatter_bound_states(D, q, r, T, K, bound_states, a_vals, aprime_vals, normconsts_or_residues, opts_slow.discretization, 0, opts->normalization_flag);
    CHECK_RETCODE(ret_code, leave_fun);

    // Update to or add residues if requested
    if (opts_slow.discspec_type != kdvv_dstype_NORMING_CONSTANTS) {

        if (opts_slow.discspec_type == kdvv_dstype_RESIDUES) {
            offset = 0;
        } else if (opts_slow.discspec_type == kdvv_dstype_BOTH) {
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
static inline INT kdvv_refine_bound_states_newton(
        const UINT D,
        COMPLEX const * const q,
        COMPLEX const * const r,
        REAL const * const T,
        const UINT K,
        COMPLEX * const bound_states,
        kdv_discretization_t const discretization,
        const UINT niter,
        REAL const * const bounding_box,
        const INT normalization_flag)
{
    INT ret_code = SUCCESS;
    UINT i, iter;
    COMPLEX a_val, b_val, aprime_val, error;
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
            ret_code = kdv_scatter_bound_states(D, q, r, T, 1,
                    bound_states + i, &a_val, &aprime_val, &b_val, discretization, 1, normalization_flag);
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
            error = CREAL(a_val) / (I * CIMAG(aprime_val));
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
