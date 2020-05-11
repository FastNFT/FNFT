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
 * Marius Brehler (TU Dortmund) 2018.
 * Shrinivas Chimmalgi (TU Delft) 2019-2020.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__misc.h" // for misc_filter
#include "fnft_nsep.h"
#include "fnft__nse_discretization.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_roots_fftgridsearch.h"
#include "fnft__nse_scatter.h"
#include "fnft__nse_fscatter.h"
#include <string.h> // for memcpy
#ifdef DEBUG
#include <stdio.h> // for printf
#endif

static fnft_nsep_opts_t default_opts = {
    .localization = fnft_nsep_loc_MIXED,
    .filtering = fnft_nsep_filt_AUTO,
    .max_evals = 20,
    .bounding_box[0] = -FNFT_INF,
    .bounding_box[1] = FNFT_INF,
    .bounding_box[2] = -FNFT_INF,
    .bounding_box[3] = FNFT_INF,
    .normalization_flag = 1,
    .discretization = nse_discretization_2SPLIT2A
};

static const UINT oversampling_factor = 32; 

// Returns an options object for the main routine with default settings.
// See the header file for a detailed description.
fnft_nsep_opts_t fnft_nsep_default_opts()
{
    return default_opts;
}

// Auxiliary routines.
static inline INT refine_mainspec(
        const UINT D, COMPLEX const * const q,
        const REAL eps_t, const UINT K,
        COMPLEX * const mainspec, const UINT max_evals,
        REAL rhs, INT kappa, nse_discretization_t discretization);
static inline INT refine_auxspec(
        const UINT D, COMPLEX const * const q,
        const REAL eps_t, const UINT K,
        COMPLEX * const auxspec, const UINT max_evals,
        const INT kappa, nse_discretization_t discretization);
static inline INT subsample_and_refine(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr, INT skip_real_flag,
        INT warn_flags[2]);
static inline INT gridsearch(const UINT D,
        COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr,
        INT warn_flags[2]);
static inline void update_bounding_box_if_auto(const REAL eps_t,
        const REAL map_coeff, fnft_nsep_opts_t * const opts_ptr);

// Main routine.
INT fnft_nsep(const UINT D, COMPLEX const * const q,
        REAL const * const T, UINT * const K_ptr,
        COMPLEX * const main_spec, UINT * const M_ptr,
        COMPLEX * const aux_spec, REAL * const sheet_indices,
        const INT kappa, fnft_nsep_opts_t * opts_ptr)
{
    INT ret_code = SUCCESS;
    UINT i, K1, K2, M1, M2, D_effective;
    INT warn_flags[2] = { 0, 0 }; // 0 = no warning about too many points so
    // far, 1st val is for main spec, 2nd for aux
    COMPLEX *q_effective = NULL;
    
    
    // Check inputs
    if (D < 2 || D%2 == 0)
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
        
    COMPLEX m_q = q[D-1]/q[0];
  
    if (CABS(1-CABS(m_q)) > 1e-2)
        WARN("Signal q may not be (quasi-)periodic. abs(1-abs(q[D-1]/q[0]))>1e-2.");

    REAL Lam_shift = CARG(m_q)/(-2*(T[1] - T[0]));
    const REAL eps_t = (T[1] - T[0])/(D - 1);
    D_effective = D-1;
    q_effective = malloc(D_effective * sizeof(COMPLEX));
    if (q_effective == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    // Removing phase rotation along q
    for (i=0; i < D_effective; i++)
         q_effective[i] = q[i]*CEXP(2*I*Lam_shift*(T[0]+eps_t*i));
    // Last sample is dropped as it is the beginning of next period


        switch (opts_ptr->localization) {
            
            case fnft_nsep_loc_MIXED:
                
                // Store lengths of user-provided arrays in temp variables
                
                K1 = *K_ptr;
                M1 = *M_ptr;
                
                // Compute non-real points in the spectra via subsample & refine
                
                if (kappa == +1) {
                    ret_code = subsample_and_refine(D_effective, q_effective, T, &K1, main_spec,
                            &M1, aux_spec, sheet_indices, kappa, opts_ptr,
                            1 /*skip_real_flag*/, warn_flags);
                    CHECK_RETCODE(ret_code, leave_fun);
                } else { // no non-real main spec in the defocusing case, pass NULL
                    K1 = 0;
                    ret_code = subsample_and_refine(D_effective, q_effective, T, &K1, NULL,
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
                
                ret_code = gridsearch(D_effective, q_effective, T, &K2, main_spec+K1,
                        &M2, aux_spec+M1, sheet_indices, kappa, opts_ptr, warn_flags);
                CHECK_RETCODE(ret_code, leave_fun);
                
                // Signal numbers of found main and aux spec points to the user
                
                *K_ptr = K1 + K2;
                *M_ptr = M1 + M2;
                
                break;
                
            case fnft_nsep_loc_SUBSAMPLE_AND_REFINE:
                
                ret_code = subsample_and_refine(D_effective, q_effective, T, K_ptr, main_spec,
                        M_ptr, aux_spec, sheet_indices, kappa, opts_ptr,
                        0/*skip_real_flag*/, warn_flags);
                CHECK_RETCODE(ret_code, leave_fun);
                break;
                
            case fnft_nsep_loc_GRIDSEARCH:
                
                ret_code = gridsearch(D_effective, q_effective, T, K_ptr, main_spec,
                        M_ptr, aux_spec, sheet_indices, kappa, opts_ptr, warn_flags);
                CHECK_RETCODE(ret_code, leave_fun);
                break;
                
            default:
                
                return E_INVALID_ARGUMENT(opts_ptr->discretization);
        }
    for (i=0; i < *K_ptr; i++)
         main_spec[i] += Lam_shift;
    for (i=0; i < *M_ptr; i++)
         aux_spec[i] += Lam_shift;
  leave_fun:
        free(q_effective);
        return ret_code;
            
}

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
    UINT M;
    UINT i, j, D_scale, D_effective;// D_scale*D gives the effective number of samples
    INT ret_code = SUCCESS;
    COMPLEX *q_effective = NULL;
    COMPLEX *q_1 = NULL;
    COMPLEX *q_2 = NULL;
    // Check inputs
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, "Pass NULL");
    
    D_scale = nse_discretization_D_scale(opts_ptr->discretization);
    D_effective = D * D_scale;
    
    // Determine step size
    const REAL eps_t = (T[1] - T[0])/D;
    const REAL scl_factor = SQRT(3.0)/6.0;


    if (D_scale == 2) {
        q_1 = malloc(D * sizeof(COMPLEX));
        q_2 = malloc(D * sizeof(COMPLEX));
        q_effective = malloc(D_effective * sizeof(COMPLEX));
        if (q_1 == NULL || q_2 == NULL || q_effective == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }
    }

    q_effective = malloc(D_effective * sizeof(COMPLEX));
    if (q_effective == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Create a resampled version of q if higher-order method is requested
    if (D_scale == 2) {
        ret_code = misc_resample(D, eps_t, q, -eps_t*scl_factor, q_1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = misc_resample(D, eps_t, q, eps_t*scl_factor, q_2);
        CHECK_RETCODE(ret_code, release_mem);

        j = 0;
        for (i=0; i < D; i++) {
            q_effective[j] = (q_1[i]+q_2[i])/4.0 - (q_2[i]-q_1[i])*scl_factor;
            j = j+2;
        }
        j = 1;
        for (i=0; i < D; i++) {
            q_effective[j] = (q_1[i]+q_2[i])/4.0 + (q_2[i]-q_1[i])*scl_factor;
            j = j+2;
        }

    } else if (D_scale == 1) {
        memcpy(q_effective, q, D * sizeof(COMPLEX));
    }

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
    ret_code = nse_fscatter(D_effective, q_effective, eps_t, kappa, transfer_matrix, &deg,
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
            p[i] = transfer_matrix[i] + transfer_matrix[3*(deg+1)+i]; //CONJ(transfer_matrix[deg-i]);
        
        p[deg/2] += 2.0 * POW(2.0, -W); // the pow arises because

        // Find the roots of p(z)
        K = oversampling_factor*deg;
        ret_code = poly_roots_fftgridsearch(deg, p, &K, PHI, roots);
        CHECK_RETCODE(ret_code, release_mem);
        if (K > deg) {
            ret_code = E_OTHER("Found more roots than memory is available.");
            goto release_mem;
        }

        
        misc_print_buf(K,roots,"zi");
        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_z_to_lambda(K, eps_t, roots, opts_ptr->discretization);
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
        ret_code = nse_z_to_lambda(K_filtered, eps_t, roots, opts_ptr->discretization);
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
        ret_code = nse_z_to_lambda(M, eps_t, roots, opts_ptr->discretization);
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
        free(q_1);
        free(q_2);
        free(q_effective);
        
        return ret_code;
}

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
    COMPLEX * qsub = NULL;
    REAL degree1step, map_coeff;
    REAL tol_im;
    UINT deg;
    UINT Dsub = 0;
    INT W = 0, *W_ptr = NULL;
    UINT K = 0, K_filtered = 0;
    UINT M = 0;
    UINT i, j, D_scale, D_effective;// D_scale*D gives the effective number of samples
    INT ret_code = SUCCESS;
    COMPLEX *qsub_1 = NULL;
    COMPLEX *qsub_2 = NULL;
    COMPLEX *qsub_effective = NULL;
    COMPLEX *q_effective = NULL;
    COMPLEX *q_1 = NULL;
    COMPLEX *q_2 = NULL;
    UINT first_last_index[2];
    UINT nskip_per_step;
    nse_discretization_t nse_discretization = 0;
    // To suppress unused parameter warnings.
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, "Pass NULL");
    
    D_scale = nse_discretization_D_scale(opts_ptr->discretization);
    D_effective = D * D_scale;
    
    // Determine step size
    const REAL eps_t = (T[1] - T[0])/D;
    const REAL scl_factor = SQRT(3.0)/6.0;
    // Create a subsampled/resampled version of q for computing initial guesses. (The
    // refinement will be carried out based on the original signal.)
    if (D_scale == 2) {
        q_1 = malloc(D * sizeof(COMPLEX));
        q_2 = malloc(D * sizeof(COMPLEX));
        if (q_1 == NULL || q_2 == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }
    }

    Dsub = ROUND(SQRT(D * LOG2(D) * LOG2(D)));
    nskip_per_step = ROUND((REAL)D / Dsub);
    Dsub = ROUND((REAL)D / nskip_per_step); // actual Dsub
    if (Dsub%2 == 1)
        Dsub --;

    qsub_effective = malloc(Dsub * D_scale * sizeof(COMPLEX));
    q_effective = malloc(D_effective * sizeof(COMPLEX));
    if (qsub_effective == NULL || q_effective == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }   
    
    if (D_scale == 2) {
        ret_code = misc_resample(D, eps_t, q, -eps_t*scl_factor*nskip_per_step, q_1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = misc_resample(D, eps_t, q, eps_t*scl_factor*nskip_per_step, q_2);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = misc_downsample(D, q_1, &Dsub, &qsub_1, first_last_index);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = misc_downsample(D, q_2, &Dsub, &qsub_2, first_last_index);
        CHECK_RETCODE(ret_code, release_mem);

        j = 0;
        for (i=0; i < Dsub; i++) {
            qsub_effective[j] = (qsub_1[i]+qsub_2[i])/4.0 - (qsub_2[i]-qsub_1[i])*scl_factor;
            j = j+2;
        }
        j = 1;
        for (i=0; i < Dsub; i++) {
            qsub_effective[j] = (qsub_1[i]+qsub_2[i])/4.0 + (qsub_2[i]-qsub_1[i])*scl_factor;
            j = j+2;
        }

        ret_code = misc_resample(D, eps_t, q, -eps_t*scl_factor, q_1);
        CHECK_RETCODE(ret_code, release_mem);
        ret_code = misc_resample(D, eps_t, q, eps_t*scl_factor, q_2);
        CHECK_RETCODE(ret_code, release_mem);

        j = 0;
        for (i=0; i < D; i++) {
            q_effective[j] = (q_1[i]+q_2[i])/4.0 - (q_2[i]-q_1[i])*scl_factor;
            j = j+2;
        }
        j = 1;
        for (i=0; i < D; i++) {
            q_effective[j] = (q_1[i]+q_2[i])/4.0 + (q_2[i]-q_1[i])*scl_factor;
            j = j+2;
        }

        nse_discretization = nse_discretization_CF4_2;
    } else if (D_scale == 1) {
        ret_code = misc_downsample(D, q, &Dsub, &qsub_effective, first_last_index);
        CHECK_RETCODE(ret_code, release_mem);
        memcpy(q_effective, q, D * sizeof(COMPLEX));
        nse_discretization = nse_discretization_BO;
    }
    

   
    // Allocate memory for the transfer matrix
    i = nse_fscatter_numel(Dsub*D_scale, opts_ptr->discretization);
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
    ret_code = nse_fscatter(Dsub*D_scale, qsub_effective, eps_t_sub, kappa, transfer_matrix, &deg,
            W_ptr, opts_ptr->discretization);
    CHECK_RETCODE(ret_code, release_mem);
    
    // Will be required later for coordinate transforms and filtering
    degree1step = nse_discretization_degree(opts_ptr->discretization);
    if (degree1step == NAN)
        return E_INVALID_ARGUMENT(opts_ptr->discretization);
    map_coeff = 2/degree1step;
    update_bounding_box_if_auto(eps_t_sub, map_coeff, opts_ptr);
    tol_im = opts_ptr->bounding_box[1] - opts_ptr->bounding_box[0];
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
        
        p[deg/2] += 2.0 * POW(2.0, -W); // the pow arises because

        
        // Find the roots of p(z)
        ret_code = poly_roots_fasteigen(deg, p, roots);
        CHECK_RETCODE(ret_code, release_mem);
        
        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_z_to_lambda(deg, eps_t_sub, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        
        // Filter the roots
        K = deg;
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        if (skip_real_flag != 0) {
            ret_code = misc_filter_nonreal(&K, roots, tol_im);
            CHECK_RETCODE(ret_code, release_mem);
        }
        // Refine the remaining roots
        ret_code = refine_mainspec(D_effective, q_effective, eps_t, K, roots,
                opts_ptr->max_evals, +2.0, kappa, nse_discretization);
        CHECK_RETCODE(ret_code, release_mem);
        
        // Filter the refined roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        if (skip_real_flag != 0) {
            ret_code = misc_filter_nonreal(&K, roots, tol_im);
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
        ret_code = poly_roots_fasteigen(deg, p, roots);
        CHECK_RETCODE(ret_code, release_mem);
        
        // Coordinate transform of the new roots
        ret_code = nse_z_to_lambda(deg, eps_t_sub, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);
        
        // Filter the new roots
        K_filtered = deg;
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K_filtered, roots, NULL,
                    opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        if (skip_real_flag != 0) {
            ret_code = misc_filter_nonreal(&K_filtered, roots, tol_im);
            CHECK_RETCODE(ret_code, release_mem);
        }
        // Refine the remaining new roots
        ret_code = refine_mainspec(D_effective, q_effective, eps_t, K_filtered,
                roots, opts_ptr->max_evals, -2.0, kappa, nse_discretization);
        
        // Filter the refined new roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&K_filtered, roots, NULL,
                    opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        if (skip_real_flag != 0) {
            ret_code = misc_filter_nonreal(&K_filtered, roots, tol_im);
            CHECK_RETCODE(ret_code, release_mem);
        }
        
        // Copy to user-provided array
        if (K_filtered + K > *K_ptr) {
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
    
    // Compute aux spectrum if desired
    if (aux_spec != NULL) {
        ret_code = poly_roots_fasteigen(deg, transfer_matrix + (deg + 1),
                roots);
        CHECK_RETCODE(ret_code, release_mem);

        // Set number of points in the aux spectrum
        M = deg;

        // Coordinate transform (from discrete-time to continuous-time domain)
        ret_code = nse_z_to_lambda(M, eps_t_sub, roots, opts_ptr->discretization);
        CHECK_RETCODE(ret_code, release_mem);

        // Filter the roots
        if (opts_ptr->filtering != fnft_nsep_filt_NONE) {
            ret_code = misc_filter(&M, roots, NULL, opts_ptr->bounding_box);
            CHECK_RETCODE(ret_code, release_mem);
        }
        
        // Refine the roots
        ret_code = refine_auxspec(D_effective, q_effective, eps_t, M, roots,
                opts_ptr->max_evals, kappa, nse_discretization);
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
        free(qsub);
        free(qsub_effective);
        free(q_1);
        free(q_2);
        return ret_code;
}

// Uses Newton's method for roots of higher order to refine the main spectrum.
static inline INT refine_mainspec(
        const UINT D, COMPLEX const * const q,
        const REAL eps_t, const UINT K,
        COMPLEX * const mainspec,
        const UINT max_evals, const REAL rhs, const INT kappa, nse_discretization_t discretization)
{
    UINT k;
    COMPLEX M[8];
    COMPLEX lam, f, f_prime, incr, tmp, next_f, next_f_prime;
    REAL cur_abs, min_abs;
    UINT nevals;
    INT ret_code;
    UINT m, best_m;
    const UINT max_m = 4; // led to the lowest number of function evals in
    // an example
    COMPLEX * r = NULL;
    
    for (k=0; k<K; k++) { // Iterate over the provided main spectrum estimates.
        
        // Initilization. Computes the monodromy matrix at the current main
        // spectrum estimate lam and determines the value of
        // f=a(lam)+a~(lam)+rhs as well as of f' = df/dlam. (The main spectrum
        // consists of the roots of f for rhs=+/- 2.0.)
        
        ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1,
                &mainspec[k], M, discretization);
        if (ret_code != SUCCESS)
            return E_SUBROUTINE(ret_code);
        next_f = M[0] + M[3] + rhs; // f = a(lam) + atil(lam) + rhs
        next_f_prime = 2.0*( M[4] + M[7] ); // f' = 2*[a'(lam) + atil'(lam)]
        
        // Iteratively refine the current main spectrum poINT by applying
        // Newton's method for higher order roots: next_x=x-m*f/f', where
        // m is the order of the root. Since we do not know m, several values
        // are tested in a line search-like procedure. Per iterion, max_m
        // values of m are tested, leading to max_m monodromoy mat evaluations.
        for (nevals=1; nevals<=max_evals; nevals+=max_m) {
            
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
                ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1, &lam, M, discretization);
                if (ret_code != SUCCESS)
                    return E_SUBROUTINE(ret_code);
                tmp = M[0] + M[3] + rhs;
                cur_abs = CABS(tmp);
                // keep this m if the new value of |f| would be lower than the
                // ones for the previously tested values of m
                if ( cur_abs < min_abs ) {
                    min_abs = cur_abs;
                    best_m = m;
                    next_f = tmp;
                    next_f_prime = 2.0*( M[4] + M[7] );
                }
            }
            
            //printf("MAIN k=%zu, nevals=%zu: lam=%g+%gj, |f|=%g, |f_prime|=%g, |incr|=%g\n", k, nevals, CREAL(mainspec[k]), CIMAG(mainspec[k]), CABS(f), CABS(f_prime),CABS(incr));
            if ( min_abs >= CABS(f) ) // stop if no improvement
                break;
            
            mainspec[k] -= best_m*incr;
        };
        //printf("==> used %zu evaluations\n", nevals);
    }
    free(r);
    return SUCCESS;
}

// Uses Newton's method to refine the aux spectrum.
static inline INT refine_auxspec(
        const UINT D, COMPLEX const * const q,
        const REAL eps_t, const UINT K,
        COMPLEX * const auxspec, const UINT max_evals,
        const INT kappa, nse_discretization_t discretization)
{
    UINT k, nevals;
    COMPLEX M[8];
    COMPLEX f, f_prime, prev_f;
    INT ret_code;
    COMPLEX * r = NULL;
    for (k=0; k<K; k++) {
        
        prev_f = NAN;
        for (nevals=0; nevals<max_evals; nevals++) {
            
            ret_code = nse_scatter_matrix(D, q, r, eps_t, kappa, 1,
                    &auxspec[k], M, discretization);
            if (ret_code != SUCCESS)
                return E_SUBROUTINE(ret_code);
            
            f = M[1]; // f = b(lam)
            f_prime = M[5]; // f' = b'(lam)
            if (f_prime == 0.0)
                return E_DIV_BY_ZERO;
            
            if ( CABS(f) >= CABS(prev_f) ) // stop if no improvement
                break;
            
            // printf("AUX k=%zu, iter=%zu: lam=%g+%gj, |f|=%g, |f_prime|=%g\n", k, iter, CREAL(auxspec[k]), CIMAG(auxspec[k]), CABS(f), CABS(f_prime));
            
            auxspec[k] -= f / f_prime;
            prev_f = f;
        }
    }
    free(r);
    return SUCCESS;
}

static inline void update_bounding_box_if_auto(const REAL eps_t,
        const REAL map_coeff, fnft_nsep_opts_t * const opts_ptr)
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
    
    if (opts_ptr->filtering == fnft_nsep_filt_AUTO) {
        opts_ptr->bounding_box[1] = 0.9*PI / (FABS(map_coeff) * eps_t);
        opts_ptr->bounding_box[0] = -opts_ptr->bounding_box[1];
        opts_ptr->bounding_box[3] = -LOG(0.1) / (FABS(map_coeff) * eps_t);
        opts_ptr->bounding_box[2] = -opts_ptr->bounding_box[3];
    }
}

