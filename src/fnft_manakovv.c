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
 * Sander Wahls (TU Delft) 2017-2018, 2020.
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Marius Brehler (TU Dortmund) 2018.
 * Peter J Prins (TU Delft) 2020.
 * Lianne de Vries (TU Delft student) 2021.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_manakovv.h"

static fnft_manakovv_opts_t default_opts = {
	.bound_state_filtering = manakovv_bsfilt_FULL,
	.bound_state_localization = manakovv_bsloc_SUBSAMPLE_AND_REFINE,
	.niter = 10,
	.Dsub = 0, // auto
	.discspec_type = manakovv_dstype_NORMING_CONSTANTS,
	.contspec_type = manakovv_cstype_REFLECTION_COEFFICIENT,
	.normalization_flag = 1,
	.discretization = manakov_discretization_2SPLIT4B,
	.richardson_extrapolation_flag = 0
};

/**
 * Creates a new options variable for fnft_manakovv with default settings.
 * See the header file for a detailed description.
 */
fnft_manakovv_opts_t fnft_manakovv_default_opts()
{
	return default_opts;
}

/**
 * Returns the maximum number of bound states that can be detected by
 * fnft_manakovv. See header file for details.
 */
UINT fnft_manakovv_max_K(const UINT D, fnft_manakovv_opts_t const* const opts)
{
	if (opts != NULL)
		return manakov_discretization_degree(opts->discretization) * D;
	else
		return manakov_discretization_degree(default_opts.discretization) * D;
}

/**
 * Declare auxiliary routines used by the main routine fnft_manakovv.
 * Their bodies follow below.
 */

 static inline INT manakovv_compute_boundstates(
		 const UINT deg,
		 COMPLEX * const transfer_matrix,
		 const REAL eps_t,
		 UINT * const K_ptr,
		 COMPLEX * const bound_states,
		 fnft_manakovv_opts_t * const opts);
		 

static inline INT fnft_manakovv_base(
	const UINT D,
	COMPLEX const* const q1,
	COMPLEX const* const q2,
	REAL const* const T,
	const UINT M,
	COMPLEX* const contspec,
	REAL const* const XI,
	UINT* const K_ptr,
	COMPLEX* const bound_states,
	const INT kappa,
	fnft_manakovv_opts_t* opts);

static inline INT manakovv_compute_contspec(
	const UINT deg,
	const INT W,
	COMPLEX* const transfer_matrix,
	COMPLEX const* const q1,
	COMPLEX const* const q2,
	REAL const* const T,
	const UINT D,
	REAL const* const XI,
	const UINT M,
	COMPLEX* const result,
	const INT kappa,
	fnft_manakovv_opts_t const* const opts);

/**
 * Fast nonlinear Fourier transform for the Manakov
 * equation with vanishing boundary conditions.
 * This function takes care of the necessary preprocessing (for example
 * signal resampling or subsampling) based on the options before calling
 * fnft_manakovv_base. If the richardson_extrapolation_flag is set, this function
 * calls fnft_manakovv_base with half of the samples and then performs
 * Richardson extrapolation on the spectrum.
 */
INT fnft_manakovv(
	const UINT D,
	COMPLEX const* const q1,
	COMPLEX const* const q2,
	REAL const* const T,
	const UINT M,
	COMPLEX* const contspec,
	REAL const* const XI,
	UINT* const K_ptr,
	COMPLEX* const bound_states,
	COMPLEX* const normconsts_or_residues,
	const INT kappa,
	fnft_manakovv_opts_t* opts)
{
	COMPLEX* q1sub_preprocessed = NULL;
	COMPLEX* q2sub_preprocessed = NULL;
	COMPLEX* q1_preprocessed = NULL;
	COMPLEX* q2_preprocessed = NULL;
	UINT Dsub = 0;
	REAL Tsub[2] = { 0.0 ,0.0 };
	UINT first_last_index[2] = { 0 };
	UINT K_sub;
	COMPLEX* contspec_sub = NULL;
	COMPLEX* normconsts_or_residues_reserve = NULL;
	fnft_manakovv_bsloc_t bs_loc_opt = 0;
	fnft_manakovv_dstype_t ds_type_opt = 0;
	INT ret_code = SUCCESS;
	UINT i, j, upsampling_factor;

	// Check inputs
	if (D < 2)
		return E_INVALID_ARGUMENT(D);
	if (q1 == NULL)
		return E_INVALID_ARGUMENT(q1);
	if (q2 == NULL)
		return E_INVALID_ARGUMENT(q2);
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
        case manakov_discretization_2SPLIT3A:
        case manakov_discretization_2SPLIT3B:
        case manakov_discretization_2SPLIT4A:
		case manakov_discretization_2SPLIT4B:
        case manakov_discretization_2SPLIT6B:
        case manakov_discretization_4SPLIT4A:
        case manakov_discretization_4SPLIT4B:
		case manakov_discretization_4SPLIT6B:
		case manakov_discretization_FTES4_4A:
		case manakov_discretization_FTES4_suzuki:
            break;
        case manakov_discretization_BO:
		case manakov_discretization_CF4_2:
            if (opts->bound_state_localization != manakovv_bsloc_NEWTON &&
                    kappa == +1 && bound_states != NULL){
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
	upsampling_factor = manakov_discretization_upsampling_factor(opts->discretization);
	if (upsampling_factor == 0) {
		ret_code = E_INVALID_ARGUMENT(opts->discretization);
		goto leave_fun;
	}

	// Determine step size
	const REAL eps_t = (T[1] - T[0]) / (D - 1);

	if (opts->richardson_extrapolation_flag == 1) {
		ds_type_opt = opts->discspec_type;
		if (ds_type_opt == manakovv_dstype_RESIDUES) {
			opts->discspec_type = manakovv_dstype_BOTH;
			normconsts_or_residues_reserve = malloc(*K_ptr * 2 * sizeof(COMPLEX));
			if (normconsts_or_residues_reserve == NULL) {
				ret_code = E_NOMEM;
				goto leave_fun;
			}
		}
	}

	// Some higher-order discretizations require samples on a non-equidistant grid
	// while others require derivatives. The preprocessing function computes
	// the required non-equidistant samples using bandlimited interpolation and
	// the required derivatives using finite differences. The manakov_discretization_preprocess_signal
	// function also performs some further discretization specific processing.
	// Preprocessing takes care of computing things which are required by all
	// the auxiliary functions thus helping efficiency.
	Dsub = D;
	ret_code = manakov_discretization_preprocess_signal(D, q1, q2, eps_t, &Dsub, &q1_preprocessed, &q2_preprocessed,
			first_last_index, opts->discretization);
	CHECK_RETCODE(ret_code, leave_fun);

	if (contspec != NULL){
		ret_code = fnft_manakovv_base(D*upsampling_factor, q1_preprocessed, q2_preprocessed, T, M, NULL, XI, K_ptr,
			bound_states, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);

		ret_code = fnft_manakovv_base(D*upsampling_factor, q1_preprocessed, q2_preprocessed, T, M, contspec, XI, K_ptr,
			NULL, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);
	}

	if (kappa == +1 && bound_states != NULL){
		// computing bound states.
		// use only the fast eigenvalue method to localize the bound states

		Dsub = opts->Dsub;
        if (Dsub == 0) // The user wants us to determine Dsub
            Dsub = SQRT(D * LOG2(D) * LOG2(D));
        UINT nskip_per_step = ROUND((REAL)D / Dsub);
        Dsub = ROUND((REAL)D / nskip_per_step); // actual Dsub

		ret_code = manakov_discretization_preprocess_signal(D, q1, q2, eps_t, &Dsub, &q1sub_preprocessed, &q2sub_preprocessed,
                first_last_index, opts->discretization);
        CHECK_RETCODE(ret_code, leave_fun);

		ret_code = fnft_manakovv_base(Dsub*upsampling_factor, q1sub_preprocessed, q2sub_preprocessed, T, M, NULL, XI, K_ptr,
			bound_states, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);
	}

	if (opts->richardson_extrapolation_flag == 1) {
		// Allocating memory
		UINT contspec_len = 0;

		if (contspec != NULL && M > 0) {
			switch (opts->contspec_type) {
			case manakovv_cstype_BOTH:
				contspec_len = 5 * M;
				break;
			case manakovv_cstype_REFLECTION_COEFFICIENT:
				contspec_len = 2 * M;
				break;
			case manakovv_cstype_AB:
				contspec_len = 3 * M;
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
		UINT method_order;
		method_order = manakov_discretization_method_order(opts->discretization);
		if (method_order == 0) {
			ret_code = E_INVALID_ARGUMENT(discretization);
			goto leave_fun;
		}

		// The signal q is now subsampled(approx. half the samples) and
		// preprocessed as required for the discretization. This is
		// required for obtaining a second approximation of the spectrum
		// which will be used for Richardson extrapolation.
		Dsub = (UINT)CEIL(D / 2);
        if (q1sub_preprocessed != NULL)
            free(q1sub_preprocessed);
        if (q2sub_preprocessed != NULL)
            free(q2sub_preprocessed);
		ret_code = manakov_discretization_preprocess_signal(D, q1, q2, eps_t, &Dsub, &q1sub_preprocessed, &q2sub_preprocessed,
			first_last_index, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		Tsub[0] = T[0] + first_last_index[0] * eps_t;
		Tsub[1] = T[0] + first_last_index[1] * eps_t;
		const REAL eps_t_sub = (Tsub[1] - Tsub[0]) / (Dsub - 1);

		// Calling fnft_manakovv_base with subsampled signal
		bs_loc_opt = opts->bound_state_localization;		// TODO: maybe remove this and next line? Check.
		opts->bound_state_localization = manakovv_bsloc_NEWTON;
		ret_code = fnft_manakovv_base(Dsub*upsampling_factor, q1sub_preprocessed, q2sub_preprocessed, Tsub, M, contspec_sub, XI, &K_sub,
			NULL, kappa, opts);	// bound_states already calculated. Passing NULL for boundstates to avoid recalculating
		CHECK_RETCODE(ret_code, leave_fun);
		opts->bound_state_localization = bs_loc_opt;
		opts->discspec_type = ds_type_opt;

		// Richardson extrapolation of the continuous spectrum
		REAL const scl_num = POW(eps_t_sub / eps_t, method_order);
		REAL const scl_den = scl_num - 1.0;
		REAL const dxi = (XI[1] - XI[0]) / (M - 1);
		if (contspec != NULL && M > 0) {
			for (i = 0; i < M; i++) {
				if (FABS(XI[0] + dxi * i) < 0.9 * PI / (2.0 * eps_t_sub)) {
					for (j = 0; j < contspec_len; j += M)
						contspec[i + j] = (scl_num * contspec[i + j] - contspec_sub[i + j]) / scl_den;
				}
			}
		}
	}

leave_fun:
	free(q1_preprocessed);
	free(q2_preprocessed);
    free(q1sub_preprocessed);
	free(q2sub_preprocessed);
	free(contspec_sub);
	if (normconsts_or_residues_reserve != normconsts_or_residues)
		free(normconsts_or_residues_reserve);
	return ret_code;
}

// Auxiliary function: Base routine for fnft_manakovv. fnft_manakovv preprocesses the signals
// and calls this function with different options as needed. This prevents
// code doubling while being efficient.
static inline INT fnft_manakovv_base(
	const UINT D,
	COMPLEX const* const q1,
	COMPLEX const* const q2,
	REAL const* const T,
	const UINT M,
	COMPLEX* const contspec,
	REAL const* const XI,
	UINT* const K_ptr,
	COMPLEX* const bound_states,
	const INT kappa,
	fnft_manakovv_opts_t* opts)
{
	COMPLEX* transfer_matrix = NULL;
	UINT deg;
	INT W = 0, *W_ptr = NULL;
	INT ret_code = SUCCESS;
	UINT i, upsampling_factor, D_given;

	// Check inputs
	if (D < 2)
		return E_INVALID_ARGUMENT(D);
	if (q1 == NULL)
		return E_INVALID_ARGUMENT(q1);
	if (q2 == NULL)
		return E_INVALID_ARGUMENT(q2);
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
	upsampling_factor = manakov_discretization_upsampling_factor(opts->discretization);
	if (upsampling_factor == 0) {
		ret_code = E_INVALID_ARGUMENT(opts->discretization);
		goto leave_fun;
	}
	D_given = D / upsampling_factor;
	const REAL eps_t = (T[1] - T[0]) / (D_given - 1);

	i = manakov_fscatter_numel(D, opts->discretization);
	// NOTE: At this stage if i == 0 it means the discretization corresponds
	// to a slow method. Incorrect discretizations will have been checked for
	// in fnft_manakovv main

	if (i != 0) {
		//This corresponds to methods based on polynomial transfer matrix
		   // Allocate memory for the transfer matrix.
		transfer_matrix = malloc(i * sizeof(COMPLEX));  // original code
		if (transfer_matrix == NULL) {
			ret_code = E_NOMEM;
			goto leave_fun;
		}

		// Compute the transfer matrix
		if (opts->normalization_flag)
			W_ptr = &W;
		ret_code = manakov_fscatter(D, q1, q2, kappa, eps_t, transfer_matrix, &deg, W_ptr,
			opts->discretization);  // NOTE: changed D to D_given
		CHECK_RETCODE(ret_code, leave_fun);
	}
	else {
		// These indicate to the functions to follow that the discretization
		// is a method not based on polynomial transfer matrix
		deg = 0;
		W = 0;
	}

		// Compute the continuous spectrum
	if (contspec != NULL && M > 0) {
		ret_code = manakovv_compute_contspec(deg, W, transfer_matrix, q1, q2, T, D, XI, M,
			contspec, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);
	}


		// Compute the discrete spectrum
		if (kappa == +1 && bound_states != NULL) {

			// Compute the bound states
			ret_code = manakovv_compute_boundstates(deg, transfer_matrix,
					eps_t, K_ptr, bound_states, opts);
			CHECK_RETCODE(ret_code, leave_fun);

		}
		

leave_fun:
	free(transfer_matrix);
	return ret_code;
}

static inline INT manakovv_compute_boundstates(
		const UINT deg,
		COMPLEX * const transfer_matrix,
		const REAL eps_t,
		UINT * const K_ptr,
		COMPLEX * const bound_states,
		fnft_manakovv_opts_t * const opts)
{
	UINT K;
	REAL bounding_box[4] = { NAN };
	COMPLEX * buffer = NULL;
	INT ret_code = SUCCESS;
	// Set-up bounding_box, only eigenvalues in upper half complex plane
		bounding_box[0] = -INFINITY;
		bounding_box[1] = INFINITY;
		bounding_box[2] = 0.0;
		bounding_box[3] = INFINITY;

	// Localize bound states ...
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
			ret_code = manakov_discretization_z_to_lambda(K, eps_t, buffer, opts->discretization);
			CHECK_RETCODE(ret_code, leave_fun);

	// Filter bound states
	if (opts->bound_state_filtering != manakovv_bsfilt_NONE) {

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
	return ret_code;
}


// Auxiliary function: Computes continuous spectrum on a frequency grid
static inline INT manakovv_compute_contspec(
	const UINT deg,
	const INT W,
	COMPLEX* const transfer_matrix,
	COMPLEX const* const q1,
	COMPLEX const* const q2,
	REAL const* const T,
	const UINT D,
	REAL const* const XI,
	const UINT M,       // Desired number of xi for which we determine the cont spec
	COMPLEX* const result,
	const INT kappa,
	fnft_manakovv_opts_t const* const opts)
{
	COMPLEX* H11_vals = NULL, * H21_vals = NULL, * H31_vals = NULL;
	COMPLEX A, V;
	REAL scale;
	REAL phase_factor_rho, phase_factor_a, phase_factor_b;
	INT ret_code = SUCCESS;
	UINT i, offset = 0, upsampling_factor, D_given;
	COMPLEX* scatter_coeffs = NULL, * xi = NULL;

	// Determine step size
	// D is interpolated number of samples but eps_t is the step-size
	// corresponding to original number of samples.
	upsampling_factor = manakov_discretization_upsampling_factor(opts->discretization);
	if (upsampling_factor == 0) {
		ret_code = E_INVALID_ARGUMENT(opts->discretization);
		goto leave_fun;
	}
	D_given = D / upsampling_factor;
	const REAL eps_t = (T[1] - T[0]) / (D_given - 1);
	const REAL eps_xi = (XI[1] - XI[0]) / (M - 1);


	// Build xi-grid which is required for applying boundary conditions
	xi = malloc(M * sizeof(COMPLEX));
	if (xi == NULL) {
		ret_code = E_NOMEM;
		goto leave_fun;
	}
	for (i = 0; i < M; i++)
		xi[i] = XI[0] + eps_xi * i;

	// Allocate memory for transfer matrix values
	H11_vals = malloc(3 * M * sizeof(COMPLEX));
	if (H11_vals == NULL) {
		return E_NOMEM;
		goto leave_fun;
	}
	H21_vals = H11_vals + M;
	H31_vals = H21_vals + M;

	// If the discretization is a slow method then there should be no transfer_matrix
	if (deg == 0 && transfer_matrix == NULL && W == 0) {
		//Code for slow methods
		// Allocate memory for call to nse_scatter_matrix
		scatter_coeffs = malloc(9 * M * sizeof(COMPLEX));
		if (scatter_coeffs == NULL) {
			ret_code = E_NOMEM;
			goto leave_fun;
		}

		ret_code = manakov_scatter_matrix(D, q1, q2, eps_t, M, xi, kappa,
			scatter_coeffs, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		for (i = 0; i < M; i++) {
			H11_vals[i] = scatter_coeffs[i * 9];
			H21_vals[i] = scatter_coeffs[i * 9 + 3];
			H31_vals[i] = scatter_coeffs[i * 9 + 6];
		}

	}
	else {
		// Prepare the use of the chirp transform. The entries of the transfer
		// matrix that correspond to a and b will be evaluated on the frequency
		// grid xi(i) = XI1 + i*eps_xi, where i=0,...,M-1. Since
		// z=exp(2.0*I*XI*eps_t/degree1step), we find that the z at which z the transfer
		// matrix has to be evaluated are given by z(i) = 1/(A * V^-i), where:
		V = eps_xi;
		ret_code = manakov_discretization_lambda_to_z(1, eps_t, &V, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);
		A = -XI[0];
		ret_code = manakov_discretization_lambda_to_z(1, eps_t, &A, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		ret_code = poly_chirpz(deg, transfer_matrix, A, V, M, H11_vals);
		CHECK_RETCODE(ret_code, leave_fun);

		ret_code = poly_chirpz(deg, transfer_matrix + 3 * (deg + 1), A, V, M, H21_vals);
		CHECK_RETCODE(ret_code, leave_fun);

		ret_code = poly_chirpz(deg, transfer_matrix + 6 * (deg + 1), A, V, M, H31_vals);
		CHECK_RETCODE(ret_code, leave_fun);
	}

	// Compute the continuous spectrum
	switch (opts->contspec_type) {

	case manakovv_cstype_BOTH:

		offset = 2 * M;       // If we want both the reflection coefficient and the a an b coefficients, we put
							// the reflection coefficient in result[0 : 2*M-1], so we have offset 2*M s.t. a=result[2*M : 3*M-1],
							// b1 = result[3*M : 4*M-1], b2 = result[4*M : 5*M-1]

	// fall through
	case manakovv_cstype_REFLECTION_COEFFICIENT:

		ret_code = manakov_discretization_phase_factor_rho(eps_t, T[1], &phase_factor_rho, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		for (i = 0; i < M; i++) {
			if (H11_vals[i] == 0.0) {
				return E_DIV_BY_ZERO;
				goto leave_fun;
			}
			result[i] = H21_vals[i] * CEXP(I * xi[i] * phase_factor_rho) / H11_vals[i];
			result[M + i] = H31_vals[i] * CEXP(I * xi[i] * phase_factor_rho) / H11_vals[i];
		}

		if (opts->contspec_type == manakovv_cstype_REFLECTION_COEFFICIENT)
			break;
		// fall through

	case manakovv_cstype_AB:

		scale = POW(2.0, W); // needed since the transfer matrix might
		// have been scaled by manakov_fscatter. W == 0 for slow methods.

		// Calculating the discretization specific phase factors.
		ret_code = manakov_discretization_phase_factor_a(eps_t, D_given, T, &phase_factor_a, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		ret_code = manakov_discretization_phase_factor_b(eps_t, D_given, T, &phase_factor_b, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		for (i = 0; i < M; i++) {
			result[offset + i] = H11_vals[i] * scale * CEXP(I * xi[i] * phase_factor_a);
			result[offset + M + i] = H21_vals[i] * scale * CEXP(I * xi[i] * phase_factor_b);
			result[offset + 2 * M + i] = H31_vals[i] * scale * CEXP(I * xi[i] * phase_factor_b);
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
