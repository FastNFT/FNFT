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
 * Sander Wahls (TU Delft) 2017-2018, 2020-2021.
 * Marius Brehler (TU Dortmund) 2018.
 * Shrinivas Chimmalgi (TU Delft) 2019-2020.
 * Peter J Prins (TU Delft) 2020.
 * Sander Wahls (KIT) 2023.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_nsep.h"

static fnft_nsep_opts_t default_opts = {
    .localization = fnft_nsep_loc_SUBSAMPLE_AND_REFINE,
    .filtering = fnft_nsep_filt_AUTO,
    .max_evals = 50,
    .bounding_box[0] = -FNFT_INF,
    .bounding_box[1] = FNFT_INF,
    .bounding_box[2] = -FNFT_INF,
    .bounding_box[3] = FNFT_INF,
    .normalization_flag = 1,
    .discretization = nse_discretization_2SPLIT2A,
    .floquet_range = {-1, 1},
    .points_per_spine = 2,
    .Dsub = 0, // => the algorithm chooses Dsub automatically
    .tol = -1 // negative tol => the algorithm chooses tol automatically
};

static const UINT oversampling_factor = 32;

// Returns an options object for the main routine with default settings.
// See the header file for a detailed description.
fnft_nsep_opts_t fnft_nsep_default_opts()
{
    return default_opts;
}

// Declaration of auxiliary routines. The function bodies are below.
static inline INT refine_mainspec(
        const UINT D, COMPLEX const * const q, COMPLEX * r,
        const REAL eps_t, const UINT K,
        COMPLEX * const mainspec, const UINT max_evals,
        REAL rhs, const REAL tol, INT kappa, nse_discretization_t discretization,
        const REAL * bounding_box, const INT normalization_flag);

static inline INT refine_auxspec(
        const UINT D, COMPLEX const * const q, COMPLEX * r,
        const REAL eps_t, const UINT K,
        COMPLEX * const auxspec, const UINT max_evals, const REAL tol,
        const INT kappa, nse_discretization_t discretization,
        const REAL * bounding_box, const INT normalization_flag);

static inline INT subsample_and_refine(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr, INT skip_real_flag,
        INT warn_flags[2]);

static inline INT newton(
        const UINT D, COMPLEX const * const q, REAL const * const T,
        UINT * const K_ptr, COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr);

static inline INT gridsearch(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr,
        INT warn_flags[2]);

static inline void update_bounding_box_if_auto(const REAL eps_t,
        REAL map_coeff, fnft_nsep_opts_t * const opts_ptr);

// Main routine.
INT fnft_nsep(const UINT D, COMPLEX const * const q,
        REAL const * const T, REAL const phase_shift, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr)
{
    INT ret_code = SUCCESS;
    UINT i, K1, K2, M1, M2;
    INT warn_flags[2] = { 0, 0 }; // 0 = no warning about too many points so
    // far, 1st val is for main spec, 2nd for aux
    COMPLEX *q_preprocessed = NULL;

    // Check inputs
    // D (and Dsub in line 493) only need to be even to ensure that deg in line 279 and 537 are also even for all discretizations.
    // deg is required to be even so that its use in lines 321, 355 and 587 make sense.
    // Due to the periodic boundary conditions, Dsub has to be an even factor of D. For some choices like D=514, this means Dsub can only be 514 as D=514=2*257.
    // This makes the idea of subsampling unintuitive. Hence at this stage D is restricted to powers of two for which subsampling still makes sense.
    if (D < 2 || (D & (D-1)) != 0 )
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL || T[0] >= T[1])
        return E_INVALID_ARGUMENT(T);
    if (abs(kappa) != 1)
        return E_INVALID_ARGUMENT(kappa);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, Pass sheet_indices="NULL".);
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->filtering != fnft_nsep_filt_NONE && main_spec == NULL && aux_spec != NULL)
        return E_INVALID_ARGUMENT(main_spec. Filtering of the auxiliary spectrum is not possible if the main spectrum is not computed.);

    REAL Lam_shift = phase_shift/(-2*(T[1] - T[0]));
    const REAL eps_t = (T[1] - T[0])/D;
    q_preprocessed = malloc(D * sizeof(COMPLEX));
    if (q_preprocessed == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    // Removing phase rotation along q
    for (i=0; i < D; i++)
        q_preprocessed[i] = q[i]*CEXP(2*I*Lam_shift*(T[0]+eps_t*i));

    // Bounding box needs to be shifted to ensure it works properly for
    // quasi-periodic signals
    if (opts_ptr->filtering == fnft_nsep_filt_MANUAL){
        opts_ptr->bounding_box[0] -= Lam_shift;
        opts_ptr->bounding_box[1] -= Lam_shift;
    }

    switch (opts_ptr->localization) {

        case fnft_nsep_loc_MIXED:

            // Store lengths of user-provided arrays in temp variables

            K1 = *K_ptr;
            M1 = *M_ptr;

            // Compute non-real points in the spectra via subsample & refine

            if (kappa == +1) {
                ret_code = subsample_and_refine(D, q_preprocessed, T, &K1, main_spec,
                        &M1, aux_spec, sheet_indices, kappa, opts_ptr,
                        1 /*skip_real_flag*/, warn_flags);
                CHECK_RETCODE(ret_code, leave_fun);
            } else { // no non-real main spec in the defocusing case, pass NULL
                K1 = 0;
                ret_code = subsample_and_refine(D, q_preprocessed, T, &K1, NULL,
                        &M1, aux_spec, sheet_indices, kappa, opts_ptr,
                        1 /*skip_real_flag*/, warn_flags);
                CHECK_RETCODE(ret_code, leave_fun);
            }

            // Store lengths of still unused parts of user-provided arrays in temp
            // variables
            if (K1 > *K_ptr || M1 > *M_ptr)
                return E_ASSERTION_FAILED;
            K2 = *K_ptr - K1;
            M2 = *M_ptr - M1;

            // Compute real points in the spectra via gridsearch

            ret_code = gridsearch(D, q_preprocessed, T, &K2, main_spec+K1,
                    &M2, aux_spec+M1, sheet_indices, kappa, opts_ptr, warn_flags);
            CHECK_RETCODE(ret_code, leave_fun);

            // Signal numbers of found main and aux spec points to the user

            *K_ptr = K1 + K2;
            *M_ptr = M1 + M2;

            break;

        case fnft_nsep_loc_SUBSAMPLE_AND_REFINE:

            ret_code = subsample_and_refine(D, q_preprocessed, T, K_ptr, main_spec,
                    M_ptr, aux_spec, sheet_indices, kappa, opts_ptr,
                    0/*skip_real_flag*/, warn_flags);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        case fnft_nsep_loc_NEWTON:

            ret_code = newton(D, q_preprocessed, T, K_ptr, main_spec, M_ptr, aux_spec,
                              sheet_indices, kappa, opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        case fnft_nsep_loc_GRIDSEARCH:

            ret_code = gridsearch(D, q_preprocessed, T, K_ptr, main_spec,
                    M_ptr, aux_spec, sheet_indices, kappa, opts_ptr, warn_flags);
            CHECK_RETCODE(ret_code, leave_fun);
            break;

        default:

            return E_INVALID_ARGUMENT(opts_ptr->discretization);
    }
    if (main_spec != NULL) {
    for (i=0; i < *K_ptr; i++)
        main_spec[i] += Lam_shift;
    }
    if (aux_spec != NULL) {
    for (i=0; i < *M_ptr; i++)
        aux_spec[i] += Lam_shift;
    }

    // Restoring changed limits
    if (opts_ptr->filtering == fnft_nsep_filt_MANUAL){
        opts_ptr->bounding_box[0] += Lam_shift;
        opts_ptr->bounding_box[1] += Lam_shift;
    }
        leave_fun:
            free(q_preprocessed);
            return ret_code;

}

// Auxiliary function: Uses grid-search to find main and auxiliary spectrum
// on the real-line.
static inline INT gridsearch(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr, INT warn_flags[2])
{
    COMPLEX * transfer_matrix = NULL;
    COMPLEX * p = NULL;
    COMPLEX * roots = NULL;
    REAL degree1step, map_coeff;
    REAL PHI[2] = { 0.0, 2.0*PI };
    UINT deg;
    INT W = 0, *W_ptr = NULL;
    UINT K, K_filtered;
    UINT M = 0;
    UINT i, upsampling_factor, D_effective;// upsampling_factor*D gives the effective number of samples
    INT ret_code = SUCCESS;
    COMPLEX *q_preprocessed = NULL;
    COMPLEX *r_preprocessed = NULL;
    UINT Dsub = 0;
    UINT first_last_index[2] = {0};

    // Check inputs
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, "Pass NULL");

    upsampling_factor = nse_discretization_upsampling_factor(opts_ptr->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto release_mem;
    }
    D_effective = D * upsampling_factor;

    // Determine step size
    const REAL eps_t = (T[1] - T[0])/D;

    Dsub = D;
    ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
            first_last_index, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Allocate memory for the transfer matrix
    i = nse_fscatter_numel(D_effective, opts_ptr->discretization);
    if (i == 0) { // since Dsub>=2, this means unknown discretization
        ret_code = E_INVALID_ARGUMENT(opts_ptr->discretization);
        goto release_mem;
    }
    transfer_matrix = malloc(i*sizeof(COMPLEX));
    if (transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Compute the transfer matrix
    if (opts_ptr->normalization_flag)
        W_ptr = &W;
    ret_code = nse_fscatter(D_effective, q_preprocessed, r_preprocessed, eps_t, transfer_matrix, &deg,
            W_ptr, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Will be required later for coordinate transforms
    degree1step = nse_discretization_degree(opts_ptr->discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);
    map_coeff = 2/degree1step;
    update_bounding_box_if_auto(eps_t, map_coeff, opts_ptr);
    PHI[0] = map_coeff*eps_t*opts_ptr->bounding_box[0];
    PHI[1] = map_coeff*eps_t*opts_ptr->bounding_box[1];
    if (PHI[0] > PHI[1]) {
        REAL tmp = PHI[0];
        PHI[0] = PHI[1];
        PHI[1] = tmp;
    }

    roots = malloc(oversampling_factor*deg*sizeof(COMPLEX));
    if (roots == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Compute main spectrum if desired
    if (main_spec != NULL) {

        // The main spectrum is given by the z that solve Delta(z)=+/-2,
        // where Delta(z)=trace{monodromy matrix(z)}is the Floquet discriminant

        // Allocate memory for the polynomial p(z) approx z^{D/2} Delta(z)+/-2
        // and its roots
        p = malloc((deg + 1)*sizeof(COMPLEX));
        if (p == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }

        // First, determine p(z) for the positive sign (+)
        for (i=0; i<=deg; i++)
            p[i] = transfer_matrix[i] + CONJ(transfer_matrix[deg-i]);
        p[deg/2] += 2.0 * POW(2.0, -W); // the pow arises because
        // nse_fscatter rescales

        // Find the roots of p(z)
        K = oversampling_factor*deg;
        ret_code = poly_roots_fftgridsearch(deg, p, &K, PHI, roots);
        CHECK_RETCODE(ret_code, release_mem);

        if (K > deg) {
            ret_code = E_OTHER("Found more roots than memory is available.");
            goto release_mem;
        }

        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_discretization_z_to_lambda(K, eps_t, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        // Filter the roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        // Copy to user-provided array
        if (K > *K_ptr) {
            if (warn_flags[0] == 0) {
                WARN("Found more than *K_ptr main spectrum points. Returning as many as possible.");
                warn_flags[0] = 1;
            }
            K = *K_ptr;
        }
        memcpy(main_spec, roots, K * sizeof(COMPLEX));

        // Second, determine p(z) for the negative sign (-)
        p[deg/2] -= 4.0 * POW(2.0, -W);


        // Find the roots of the new p(z)
        K_filtered = oversampling_factor*deg;
        ret_code = poly_roots_fftgridsearch(deg, p, &K_filtered, PHI,
                roots);
        CHECK_RETCODE(ret_code, release_mem);

        if (K_filtered > deg) {
            ret_code = E_OTHER("Found more roots than memory is available.");
            goto release_mem;
        }

        // Coordinate transform of the new roots
        ret_code = nse_discretization_z_to_lambda(K_filtered, eps_t, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        // Filter the new roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K_filtered, roots, NULL,
                    opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        // Copy to user-provided array
        if (K + K_filtered > *K_ptr) {
            if (warn_flags[0] == 0) {
                WARN("Found more than *K_ptr main spectrum points. Returning as many as possible.");
                warn_flags[0] = 1;
            }
            if (K < *K_ptr)
                K_filtered = *K_ptr - K;
            else
                K_filtered = 0;
        }
        memcpy(main_spec + K, roots, K_filtered * sizeof(COMPLEX));

        // Set number of points in the main spectrum
        K += K_filtered;
    }

    // Compute auxiliary spectrum (real line only)
    if (aux_spec != NULL) {

        M = oversampling_factor*deg;
        ret_code = poly_roots_fftgridsearch(deg, transfer_matrix+(deg+1), &M,
                PHI, roots);
        CHECK_RETCODE(ret_code, release_mem);

        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_discretization_z_to_lambda(M, eps_t, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        // Filter the roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&M, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        // Copy to user-provided array
        if (M > *M_ptr) {
            if (warn_flags[1] == 0) {
                WARN("Found more than *M_ptr aux spectrum points. Returning as many as possible.");
                warn_flags[1] = 1;
            }
            M = *M_ptr;
        }
        memcpy(aux_spec, roots, M * sizeof(COMPLEX));
    }

    // Update number of main spectrum points and auxiliary spectrum points
    *K_ptr = K;
    *M_ptr = M;

    release_mem:
        free(transfer_matrix);
        free(p);
        free(roots);
        free(q_preprocessed);
        free(r_preprocessed);

        return ret_code;
}

// Auxiliary function: Finds initial guesses for the main and auxiliary spectrum
// using an eigenmethod on a subsampled version of the signal and then uses
// Newton's method to refine them.
static inline INT subsample_and_refine(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr, INT skip_real_flag,
        INT warn_flags[2])
{
    COMPLEX * transfer_matrix = NULL;
    COMPLEX * p = NULL;
    COMPLEX * roots = NULL;
    REAL degree1step, map_coeff;
    REAL tol_im;
    REAL refine_tol;
    UINT deg;
    UINT Dsub = 0;

    INT W = 0, *W_ptr = NULL;
    UINT K = 0, K_new = 0;
    UINT M = 0;
    UINT i, upsampling_factor, D_effective;// upsampling_factor*D gives the effective number of samples
    INT ret_code = SUCCESS;
    COMPLEX *qsub_preprocessed = NULL;
    COMPLEX *rsub_preprocessed = NULL;
    COMPLEX *q_preprocessed = NULL;
    COMPLEX *r_preprocessed = NULL;
    UINT first_last_index[2];
    UINT nskip_per_step;
    nse_discretization_t nse_discretization = 0;

    // To suppress unused parameter warnings.
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, "Pass NULL");

    upsampling_factor = nse_discretization_upsampling_factor(opts_ptr->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto release_mem;
    }
    D_effective = D * upsampling_factor;

    // Determine step size
    const REAL eps_t = (T[1] - T[0])/D;

    // Create the signal required for refinement of the initial guesses.
    Dsub = D;
    ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
            first_last_index, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Create a subsampled/resampled version of q for computing initial guesses.
    Dsub = opts_ptr->Dsub;
    if (Dsub == 0) // users wants Dsub to be chosen automatically
        Dsub = (UINT) POW(2.0, CEIL( 0.5 * LOG2(D * LOG2(D) * LOG2(D)) ));
    else
        Dsub = (UINT) POW(2.0, ROUND(LOG2(Dsub)));

    ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &qsub_preprocessed, &rsub_preprocessed,
            first_last_index, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    nskip_per_step = D/Dsub;
    if ( first_last_index[0] != 0 || first_last_index[1]+nskip_per_step != D )
        return E_ASSERTION_FAILED; // Correct update of T for general
                                   // downsampling still needs to be
                                 // implemented

    if (upsampling_factor == 2) {
        nse_discretization = nse_discretization_CF4_2;
    } else if (upsampling_factor == 1) {
        nse_discretization = nse_discretization_BO;
    }

    // Determine the tolerance for the refinement steps
    if (opts_ptr->tol < 0)
        refine_tol = SQRT(EPSILON);
    else
        refine_tol = opts_ptr->tol;

    // Allocate memory for the transfer matrix
    i = nse_fscatter_numel(Dsub*upsampling_factor, opts_ptr->discretization);
    if (i == 0) { // since Dsub>=2, this means unknown discretization
        ret_code = E_INVALID_ARGUMENT(opts_ptr->discretization);
        goto release_mem;
    }
    transfer_matrix = malloc(i*sizeof(COMPLEX));
    if (transfer_matrix == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Determine step size
    const REAL eps_t_sub = nskip_per_step*eps_t;

    // Compute the transfer matrix
    if (opts_ptr->normalization_flag)
        W_ptr = &W;
    ret_code = nse_fscatter(Dsub*upsampling_factor, qsub_preprocessed, rsub_preprocessed, eps_t_sub, transfer_matrix, &deg,
            W_ptr, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    // Will be required later for coordinate transforms and filtering
    degree1step = nse_discretization_degree(opts_ptr->discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);
    map_coeff = 2/degree1step;
    update_bounding_box_if_auto(eps_t_sub, map_coeff, opts_ptr);
    tol_im = opts_ptr->bounding_box[3] - opts_ptr->bounding_box[2];
    tol_im /= oversampling_factor*(D - 1);

    // Use unused part of transfer matrix as buffer for the roots
    roots = transfer_matrix + 3*(deg+1);

    // Compute main spectrum if desired
    if (main_spec != NULL) {

        // The main spectrum is given by the z that solve Delta(z)=+/-2,
        // where Delta(z)=trace{monodromy matrix(z)}is the Floquet discriminant

        // Allocate memory for the polynomial p(z) approx z^{D/2} Delta(z)+/-2
        p = malloc((deg + 1)*sizeof(COMPLEX));
        if (p == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }

        // First, determine p(z) for the positive sign (+)
        for (i=0; i<=deg; i++)
            p[i] = transfer_matrix[i] + CONJ(transfer_matrix[deg-i]);

        // Below, we'll search for solutions of Delta(z)=rhs, where rhs
        // iterates through a grid of nvals values between rhs_0 and rhs_1.
        // Normally, we use rhs_0=-2, rhs_1=2 and points_per_spine=2, which
        // provides us the main spectrum (= endpoints of spines). Using higher
        // values for points_per_spine allows us to compute spines.
        const REAL rhs_0 = opts_ptr->floquet_range[0];
        const REAL rhs_1 = opts_ptr->floquet_range[1];
        const UINT nvals = opts_ptr->points_per_spine;
        REAL rhs_step = rhs_1 - rhs_0;
        if (nvals>1)
            rhs_step /= nvals - 1;
        const COMPLEX p_center_coeff = p[deg/2];

        for (UINT nval=0; nval<nvals; nval++) {

            const REAL rhs = 2.0*(rhs_0 + nval*rhs_step);

            p[deg/2] = p_center_coeff - rhs * POW(2.0, -W); // the pow arises
            // because nse_fscatter rescales

            // Find the roots of p(z)-rhs
            ret_code = poly_roots_fasteigen(deg, p, roots);
            CHECK_RETCODE(ret_code, release_mem);

            // Coordinate transform (from discrete-time to continuous-time domain)
            ret_code = nse_discretization_z_to_lambda(deg, eps_t_sub, roots,
                    opts_ptr->discretization);
            CHECK_RETCODE(ret_code, release_mem);

            // Filter the roots
            K_new = deg;
            if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
                ret_code = misc_filter(&K_new, roots, NULL,
                        opts_ptr->bounding_box);
                CHECK_RETCODE(ret_code, release_mem);
            }
            if (skip_real_flag != 0) {
                ret_code = misc_filter_nonreal(&K_new, roots, tol_im);
                CHECK_RETCODE(ret_code, release_mem);
            }

            // Refine the remaining roots
            ret_code = refine_mainspec(D_effective, q_preprocessed, r_preprocessed, eps_t, K_new, roots,
                                       opts_ptr->max_evals, -rhs, refine_tol, kappa, nse_discretization,
                                       opts_ptr->bounding_box, opts_ptr->normalization_flag);
            CHECK_RETCODE(ret_code, release_mem);

            // Filter the refined roots
            if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
                ret_code = misc_filter(&K_new, roots, NULL,
                        opts_ptr->bounding_box);
                CHECK_RETCODE(ret_code, release_mem);
            }
            if (skip_real_flag != 0) {
                ret_code = misc_filter_nonreal(&K_new, roots, tol_im);
                CHECK_RETCODE(ret_code, release_mem);
            }


            // Copy to user-provided array
            if (K+K_new > *K_ptr) {
                if (warn_flags[0] == 0) {
                    WARN("Found more than *K_ptr main spectrum points. Returning as many as possible.");
                    warn_flags[0] = 1;
                }
                K_new = *K_ptr-K;
            }

            memcpy(main_spec+K, roots, K_new * sizeof(COMPLEX));
            K += K_new;
            if (warn_flags[0] == 1)
                break; // user-provided array for main spectrum is full
        }

    }

    // Compute aux spectrum if desired
    if (aux_spec != NULL) {
        ret_code = poly_roots_fasteigen(deg, transfer_matrix + (deg + 1),
                roots);


        // Set number of points in the aux spectrum
        M = deg;

        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_discretization_z_to_lambda(M, eps_t_sub, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        // Filter the roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&M, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        // Refine the roots
        ret_code = refine_auxspec(D_effective, q_preprocessed, r_preprocessed, eps_t, M, roots,
                                  opts_ptr->max_evals, refine_tol, kappa, nse_discretization,
                                  opts_ptr->bounding_box, opts_ptr->normalization_flag);
        CHECK_RETCODE(ret_code, release_mem);


        // Filter the refined roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&M, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        if (skip_real_flag != 0) {
            ret_code = misc_filter_nonreal(&M, roots, tol_im);
            CHECK_RETCODE(ret_code, release_mem);
        }

        // Copy to user-provided array
        if (M > *M_ptr) {
            if (warn_flags[1] == 0) {
                WARN("Found more than *M_ptr aux spectrum points. Returning as many as possible.");
                warn_flags[1] = 1;
            }
            M = *M_ptr;
        }
        memcpy(aux_spec, roots, M * sizeof(COMPLEX));

    }

    // Update number of main spectrum points and auxiliary spectrum points
    *K_ptr = K;
    *M_ptr = M;

    release_mem:
        free(transfer_matrix);
        free(p);
        free(q_preprocessed);
        free(r_preprocessed);
        free(qsub_preprocessed);
        free(rsub_preprocessed);
        return ret_code;
}

// Auxiliary function: Refines user-provided initial guesses for the main and
// auxiliary spectrum using Newton's method.
static inline INT newton(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr)
{
    COMPLEX * initial_guesses = NULL;
    REAL degree1step, map_coeff;
    REAL tol_im;
    REAL refine_tol;
    UINT Dsub = 0;
    UINT K = 0;
    UINT M = 0;
    COMPLEX *q_preprocessed = NULL;
    COMPLEX *r_preprocessed = NULL;
    UINT first_last_index[2];
    INT ret_code = SUCCESS;

    // To suppress unused parameter warnings.
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, "Pass NULL");

    if (K_ptr != NULL)
        K = *K_ptr;

    if (M_ptr != NULL)
        M = *M_ptr;

    const UINT upsampling_factor = nse_discretization_upsampling_factor(opts_ptr->discretization);
    if (upsampling_factor == 0) {
        ret_code = E_INVALID_ARGUMENT(opts->discretization);
        goto release_mem;
    }
    const UINT D_effective = D * upsampling_factor;

    // Determine step size
    const REAL eps_t = (T[1] - T[0])/D;

    // Create the signal required for refinement of the initial guesses.
    Dsub = D;
    ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
            first_last_index, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);

    if (Dsub != D) {
        ret_code = E_ASSERTION_FAILED;
        goto release_mem;
    }

    // Determine the tolerance for the refinement steps
    if (opts_ptr->tol < 0)
        refine_tol = SQRT(EPSILON);
    else
        refine_tol = opts_ptr->tol;

    // Will be required later for coordinate transforms and filtering
    degree1step = nse_discretization_degree(opts_ptr->discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);
    map_coeff = 2/degree1step;
    update_bounding_box_if_auto(eps_t, map_coeff, opts_ptr);
    tol_im = opts_ptr->bounding_box[1] - opts_ptr->bounding_box[0];
    tol_im /= oversampling_factor*(D - 1);

    // Refine main spectrum if desired
    if (main_spec != NULL) {

        // The main spectrum is given by the z that solve Delta(z)=+/-2,
        // where Delta(z)=trace{monodromy matrix(z)} is the Floquet discriminant

        // Below, we'll search for solutions of Delta(z)=rhs, where rhs
        // iterates through a grid of nvals values between rhs_0 and rhs_1.
        // Normally, we use rhs_0=-2, rhs_1=2 and points_per_spine=2, which
        // provides us the main spectrum (= endpoints of spines). Using higher
        // values for points_per_spine allows us to compute spines.
        const REAL rhs_0 = opts_ptr->floquet_range[0];
        const REAL rhs_1 = opts_ptr->floquet_range[1];
        const UINT nvals = opts_ptr->points_per_spine;
        REAL rhs_step = rhs_1 - rhs_0;
        if (nvals>1)
            rhs_step /= nvals - 1;
        initial_guesses = main_spec;

        for (UINT nval=0; nval<nvals; nval++) {

            const REAL rhs = 2.0*(rhs_0 + nval*rhs_step);

            // Refine the initial guesses
            ret_code = refine_mainspec(D_effective, q_preprocessed, r_preprocessed, eps_t,
                                       K, initial_guesses, opts_ptr->max_evals,
                                       -rhs, refine_tol, kappa,
                                       opts_ptr->discretization, opts_ptr->bounding_box,
                                       opts_ptr->normalization_flag);
            CHECK_RETCODE(ret_code, release_mem);

            initial_guesses += K;
        }

        K *= opts_ptr->points_per_spine;

        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K, main_spec, NULL,
                                   opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        *K_ptr = K;
    }

    // Refine aux spectrum if desired
    if (aux_spec != NULL) {

        // Refine the roots
        ret_code = refine_auxspec(D_effective, q_preprocessed, r_preprocessed, eps_t, M,
                                  aux_spec, opts_ptr->max_evals, refine_tol,
                                  kappa, opts_ptr->discretization, opts_ptr->bounding_box,
                                  opts_ptr->normalization_flag);
        CHECK_RETCODE(ret_code, release_mem);

        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&M, aux_spec, NULL,
                                   opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }

        *M_ptr = M;
    }

release_mem:
    free(q_preprocessed);
    free(r_preprocessed);
    return ret_code;
}

// Auxiliary function: Checks if a complex number is in a given box
static INT in_box(const COMPLEX z, const REAL * box)
{
    return (CREAL(z) >= box[0]) && (CREAL(z) <= box[1]) &&
           (CIMAG(z) >= box[2]) && (CIMAG(z) <= box[3]);
}

// Auxiliary function: Uses Newton's method for roots of higher order to refine the main spectrum.
static inline INT refine_mainspec(
        const UINT D, COMPLEX const * const q, COMPLEX * r,
        const REAL eps_t, const UINT K,
        COMPLEX * const mainspec,
        const UINT max_evals, const REAL rhs, const REAL tol,
        const INT kappa, nse_discretization_t discretization,
        const REAL * bounding_box, const INT normalization_flag)
{
    UINT k;
    COMPLEX M[8];
    COMPLEX lam, f, f_prime, incr, tmp, next_f, next_f_prime, best_incr;
    REAL cur_abs, min_abs;
    UINT nevals;
    INT ret_code;
    UINT m, best_m;
    const UINT max_m = 2; // roots might be single or double
    INT warn_flag = 0;
    INT W = 0;
    INT * const W_ptr = (normalization_flag) ? &W : NULL;

    if (max_evals == 0)
        return SUCCESS;

    for (k=0; k<K; k++) { // Iterate over the provided main spectrum estimates.

        // Initialization. Computes the monodromy matrix at the current main
        // spectrum estimate lam and determines the value of
        // f=a(lam)+a~(lam)+rhs as well as of f' = df/dlam. (The main spectrum
        // consists of the roots of f for rhs=+/- 2.0.)

        ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1,
                &mainspec[k], M, W_ptr, discretization, 1);
        if (ret_code != SUCCESS)
            return E_SUBROUTINE(ret_code);

        REAL scl = (normalization_flag) ? POW(2, W) : 1;
        next_f = scl*(M[0] + M[3]) + rhs; // f = a(lam) + atil(lam) + rhs
        next_f_prime = scl*(M[4] + M[7]); // f' = a'(lam) + atil'(lam)

        // Iteratively refine the current main spectrum point by applying
        // Newton's method for higher order roots: next_x=x-m*f/f', where
        // m is the order of the root. Since we do not know m, several values
        // are tested in a line search-like procedure. Per iteration, max_m
        // values of m are tested, leading to max_m monodromoy mat evaluations.
        for (nevals=1; nevals<=max_evals;) {

            // The current values of f and f' at lam = mainspec[k]
            f = next_f;
            f_prime = next_f_prime;
            if (f_prime == 0.0)
                return E_DIV_BY_ZERO;
            incr = f / f_prime;

	        // Test different increments to deal with the many higher order
            // roots (Newton's method for a root of order m is x<-x-m*f/f').
            min_abs = INFINITY;
            best_m = 1;
            for (m=1; m<=max_m; m++) {
                lam = mainspec[k] - m*incr;
                ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1, &lam, M, W_ptr, discretization, 1);
                if (ret_code != SUCCESS)
                    return E_SUBROUTINE(ret_code);
                scl = (normalization_flag) ? POW(2, W) : 1;
                nevals++;
                tmp = scl*(M[0] + M[3]) + rhs;
                cur_abs = CABS(tmp);
                // keep this m if the new value of |f| would be lower than the
                // ones for the previously tested values of m
                if ( cur_abs < min_abs ) {
                    min_abs = cur_abs;
                    best_m = m;
                    next_f = tmp;
                    next_f_prime = scl*(M[4] + M[7]);
                    if ( cur_abs < tol )
                        break;
                }
            }
            best_incr = best_m*incr;
            mainspec[k] -= best_incr; // Newton step

	        // printf("MAIN k=%zu, nevals=%zu: lam=%e+%ej, |f|=%e, |f_prime|=%e, |incr|=%e, best_m=%d\n", k, nevals, CREAL(mainspec[k]), CIMAG(mainspec[k]), CABS(f), CABS(f_prime),CABS(incr),best_m);

            if ( CABS(best_incr) <= tol || CABS(f) <= tol ) {
                // We already know f and f_prime at the new mainspec[k], so
                // let's use that for a final first-order Newton step.
                if (next_f_prime == 0.0)
                    return E_DIV_BY_ZERO;
                mainspec[k] -= next_f / next_f_prime;
                break;
            }

            // stop refining the current value if it left the bounding box
            if (!in_box(mainspec[k], bounding_box))
                break;
        }

	    if (CABS(best_incr) > tol && CABS(f) > tol && in_box(mainspec[k], bounding_box))
		    warn_flag = 1;

        //printf("==> used %zu evaluations\n", nevals);
    }

    if (warn_flag)
    	WARN("Refinement of one or more main spectrum (or spine) points did not converge. Try manual filtering, or increase max_evals or tol.")

    return SUCCESS;
}

// Auxiliary function: Uses Newton's method to refine the aux spectrum.
static inline INT refine_auxspec(
        const UINT D, COMPLEX const * const q, COMPLEX * r,
        const REAL eps_t, const UINT K,
        COMPLEX * const auxspec, const UINT max_evals, const REAL tol,
        const INT kappa, nse_discretization_t discretization,
        const REAL * bounding_box, const INT normalization_flag)
{
    UINT k, nevals;
    COMPLEX M[8];
    COMPLEX f, f_prime, incr;
    INT warn_flag = 0;
    INT ret_code;
    INT W = 0;
    INT * const W_ptr = (normalization_flag) ? &W : NULL;

    if (max_evals == 0)
        return SUCCESS;

    for (k=0; k<K; k++) {

        for (nevals=0; nevals<max_evals;) {

            ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1,
                    &auxspec[k], M, W_ptr, discretization, 1);
            if (ret_code != SUCCESS)
                return E_SUBROUTINE(ret_code);
            nevals++;
            // We should in principle multiply all entries of M with scl
            // due to the normalization, but we don't because we only
            // care about the ratio of f and f_prime, so that scl cancels out
            f = M[1]; // f = b(lam)
            f_prime = M[5]; // f' = b'(lam)
            if (f_prime == 0.0)
                return E_DIV_BY_ZERO;

            //printf("AUX k=%zu, iter=%zu: lam=%g+%gj, |f|=%g, |f_prime|=%g\n", k, iter, CREAL(auxspec[k]), CIMAG(auxspec[k]), CABS(f), CABS(f_prime));

            incr = f / f_prime;
            auxspec[k] -= incr;
            if ( CABS(incr) <= tol || CABS(f) <= tol ) // Intentionally put after the previous line
                // => we use the already known values for f and f_prime for a
                // last Newton step even if already |f|<tol.
                break;

            // stop refining the current value if it left the bounding box
            if (!in_box(auxspec[k], bounding_box))
                break;
        }

    	if (CABS(incr) > tol && CABS(f) > tol && in_box(auxspec[k], bounding_box))
    	    warn_flag = 1;
    }

    if (warn_flag)
    	WARN("Refinement of one or more auxiliary spectrum points did not converge. Try manual filtering, or increase max_evals or tol.")

    return SUCCESS;
}

static inline void update_bounding_box_if_auto(const REAL eps_t,
        REAL map_coeff, fnft_nsep_opts_t * const opts_ptr)
{
    // The bounding box used for filtering is computed. The format of the box
    // is
    //
    //  bounding_box=[lower_bound_re, upper_bound_re, lower_bound_im,
    //      upper_bound_im].
    //
    // Real part:
    //
    // Since the coordinate transform
    //
    //  z = exp(1j*map_coeff*eps_t*lam)
    //
    // is periodic for real lam with period 2*pi/map_coeff, we can resolve
    //
    //  -pi/|map_coeff|/eps_t <= Re(lam) <= pi/|map_coeff|/eps_t.
    //
    // We enforce |Re(lam)|<=0.9*pi/|map_coeff|/eps_t.

    // Deal with "slow" discretizations for which no coordinate transform
    // of the form above is used. By setting map_coeff=1, we fall back to
    // the conventional Nyquist frequency in this case.
    if (map_coeff == INFINITY || map_coeff == 0 || map_coeff != map_coeff /* NaN */)
        map_coeff = 1.0;

    if (opts_ptr->filtering == fnft_nsep_filt_AUTO) {
        opts_ptr->bounding_box[1] = 0.9*PI / (FABS(map_coeff) * eps_t);
        opts_ptr->bounding_box[0] = -opts_ptr->bounding_box[1];
        opts_ptr->bounding_box[3] = -LOG(0.1) / (FABS(map_coeff) * eps_t);
        opts_ptr->bounding_box[2] = -opts_ptr->bounding_box[3];
    }
}
