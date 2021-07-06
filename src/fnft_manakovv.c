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
 // NOTE: this is degree * number of samples. Makes sense, this gives the degree of the resulting polynomial, so the max. amount of roots.
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
		 UINT const D,
		 COMPLEX const * const q1,
		 COMPLEX const * const q2,
		 const UINT deg,
		 COMPLEX * const transfer_matrix,
		 REAL const * const T,
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
	COMPLEX* const normconsts_or_residues,
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

/*
static inline INT manakovv_compute_normconsts_or_residues(
		const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
		REAL const * const T,
		const UINT K,
		COMPLEX * const bound_states,
		COMPLEX * const normconsts_or_residues,
		fnft_manakovv_opts_t const * const opts);

static inline INT manakovv_refine_bound_states_newton(const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
		REAL const * const T,
		UINT K,
		COMPLEX * const bound_states,
		akns_discretization_t const discretization,
		UINT const niter,
		REAL const * const bounding_box);
*/
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
	UINT Dsub = 0;
	REAL Tsub[2] = { 0.0 ,0.0 };
	UINT first_last_index[2] = { 0 };
	UINT K_sub;
	COMPLEX* contspec_sub = NULL;
	COMPLEX* bound_states_sub = NULL;
	COMPLEX* normconsts_or_residues_sub = NULL;
	COMPLEX* normconsts_or_residues_reserve = NULL;
	fnft_manakovv_bsloc_t bs_loc_opt = 0;
	fnft_manakovv_dstype_t ds_type_opt = 0;
	INT ret_code = SUCCESS;
	UINT i, j, upsampling_factor, D_effective, nskip_per_step;

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


	/*    // This switch checks for incompatible bound_state_localization options
		switch (opts->discretization) {
			case akns_discretization_2SPLIT2_MODAL:
			case akns_discretization_2SPLIT1A:
			case akns_discretization_2SPLIT1B:
			case akns_discretization_2SPLIT2A:
			case akns_discretization_2SPLIT2B:
			case akns_discretization_2SPLIT2S:
			case akns_discretization_2SPLIT3S:
			case akns_discretization_2SPLIT4B:
			case akns_discretization_2SPLIT3A:
			case akns_discretization_2SPLIT3B:
			case akns_discretization_2SPLIT4A:
			case akns_discretization_2SPLIT6B:
			case akns_discretization_2SPLIT6A:
			case akns_discretization_2SPLIT8B:
			case akns_discretization_2SPLIT5A:
			case akns_discretization_2SPLIT5B:
			case akns_discretization_2SPLIT8A:
			case akns_discretization_2SPLIT7A:
			case akns_discretization_2SPLIT7B:
			case akns_discretization_4SPLIT4A:
			case akns_discretization_4SPLIT4B:
				break;
			case akns_discretization_BO:
			case akns_discretization_CF4_2:
			case akns_discretization_CF4_3:
			case akns_discretization_CF5_3:
			case akns_discretization_CF6_4:
			case akns_discretization_ES4:
			case akns_discretization_TES4:
				if (opts->bound_state_localization != manakovv_bsloc_NEWTON &&
						kappa == +1 && bound_states != NULL){
					ret_code = E_INVALID_ARGUMENT(opts->bound_state_localization);
					goto leave_fun;
				}
				break;
			default: // Unknown discretization
				return E_INVALID_ARGUMENT(opts->discretization);
		}*/

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
	// the required derivatives using finite differences. The nse_discretization_preprocess_signal
	// function also performs some further discretization specific processing.
	// Preprocessing takes care of computing things which are required by all
	// the auxiliary functions thus helping efficiency.
	/*
	Dsub = D;
	ret_code = nse_discretization_preprocess_signal(D, q, eps_t, kappa, &Dsub, &q_preprocessed, &r_preprocessed,
			first_last_index, opts->discretization);
	CHECK_RETCODE(ret_code, leave_fun);

	*/

	Dsub = D;
	ret_code = manakov_discretization_preprocess_signal(D, q1, q2, eps_t, kappa, &Dsub, &q1sub_preprocessed, &q2sub_preprocessed,
			first_last_index, opts->discretization);
	CHECK_RETCODE(ret_code, leave_fun);

	// if kappa=-1 there are no discrete eigenvalues
	if (kappa == -1){
			WARN("No discrete eigenvalues exist for kappa=-1. Not calculating bound states");
	}
		// use only the fast eigenvalue method to localize the bound states
		ret_code = fnft_manakovv_base(D*upsampling_factor, q1sub_preprocessed, q2sub_preprocessed, T, M, contspec, XI, K_ptr,
			bound_states, normconsts_or_residues, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);
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
		/* code relating to disc. spec.
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
		*/
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
		printf("D for preprocess_signal = %d\n", D);
		printf("Dsub for preprocess_signal = %d\n", Dsub);
		ret_code = manakov_discretization_preprocess_signal(D, q1, q2, eps_t, kappa, &Dsub, &q1sub_preprocessed, &q2sub_preprocessed,
			first_last_index, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		Tsub[0] = T[0] + first_last_index[0] * eps_t;
		Tsub[1] = T[0] + first_last_index[1] * eps_t;
		const REAL eps_t_sub = (Tsub[1] - Tsub[0]) / (Dsub - 1);

		// Calling fnft_manakovv_base with subsampled signal
		bs_loc_opt = opts->bound_state_localization;
		opts->bound_state_localization = manakovv_bsloc_NEWTON;
		ret_code = fnft_manakovv_base(Dsub*upsampling_factor, q1sub_preprocessed, q2sub_preprocessed, Tsub, M, contspec_sub, XI, &K_sub,
			NULL, normconsts_or_residues_sub, kappa, opts);	// bound_states already calculated. Passing NULL for boundstates to avoid recalculating
		misc_print_buf(100, contspec_sub, "contspec_sub");
		CHECK_RETCODE(ret_code, leave_fun);
		opts->bound_state_localization = bs_loc_opt;
		opts->discspec_type = ds_type_opt;

		// Richardson extrapolation of the continuous spectrum
		REAL const scl_num = POW(eps_t_sub / eps_t, method_order);     // TODO: order of method is known, but maybe it is a good idea to calculate it here and verify if it matches?
		REAL const scl_den = scl_num - 1.0;
		REAL const dxi = (XI[1] - XI[0]) / (M - 1);
		if (contspec != NULL && M > 0) {
			for (i = 0; i < M; i++) {
				if (FABS(XI[0] + dxi * i) < 0.9 * PI / (2.0 * eps_t_sub)) {    // NOTE: What does this line do? (Nyquist-Shannon sampling theorem, outside these bounds we get aliasing)
					for (j = 0; j < contspec_len; j += M)
						contspec[i + j] = (scl_num * contspec[i + j] - contspec_sub[i + j]) / scl_den;
				}
			} // NOTE: Richardson extrapolation is already applied to all values in array M. So as long as we put both b coefficients (or reflection coefs) in M, we do not need exta code to take care of that
		}

		/* Code for disc. spec.
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
		*/
	}

leave_fun:      // TODO: check which variables we need to free
//       free(q1sub_preprocessed);
//        free(q2sub_preprocessed);
//        free(contspec_sub);
//        free(bound_states_sub);
//        free(normconsts_or_residues_sub);
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
	COMPLEX* const normconsts_or_residues,
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

	// D should be the effective number of samples in q
	i = manakov_fscatter_numel(D, opts->discretization);
	// NOTE: At this stage if i == 0 it means the discretization corresponds
	// to a slow method. Incorrect discretizations will have been checked for
	// in fnft_nsev main

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
		/*            printf("inputs to manakov f_scatter\n");
					printf("D_given = %d\n",D_given);
					misc_print_buf(16,q1,"first_16_entries_q1");
					misc_print_buf(16,q2,"first_16_entries_q2");
					printf("kappa = %d\n",kappa);
					printf("eps_t = %f\n", eps_t);
					misc_print_buf(16,transfer_matrix,"first_16_entries_transfer_matrix");
					*/
		misc_print_buf(100,q1,"q1_new");
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

	//        misc_print_buf(10, transfer_matrix, "tm before manakovv_compute_contspec");
		// Compute the continuous spectrum
	if (contspec != NULL && M > 0) {
		ret_code = manakovv_compute_contspec(deg, W, transfer_matrix, q1, q2, T, D, XI, M,
			contspec, kappa, opts);
		CHECK_RETCODE(ret_code, leave_fun);
	}

	//        misc_print_buf(10, transfer_matrix, "tm after manakovv_compute_contspec");

		// Compute the discrete spectrum
		if (kappa == +1 && bound_states != NULL) {

			// Compute the bound states
			ret_code = manakovv_compute_boundstates(D, q1, q2, deg, transfer_matrix, T,
					eps_t, K_ptr, bound_states, opts);
			CHECK_RETCODE(ret_code, leave_fun);

		/*	// Norming constants and/or residues)
			if (normconsts_or_residues != NULL && *K_ptr != 0) {
				ret_code = nsev_compute_normconsts_or_residues(D, q, r, T, *K_ptr,
						bound_states, normconsts_or_residues, opts);
				CHECK_RETCODE(ret_code, leave_fun);
			}
		} else if (K_ptr != NULL) {
			*K_ptr = 0;*/
		}
		

leave_fun:
	free(transfer_matrix);
	return ret_code;
}

static inline INT manakovv_compute_boundstates(
		const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
		const UINT deg,
		COMPLEX * const transfer_matrix,
		REAL const * const T,
		const REAL eps_t,
		UINT * const K_ptr,
		COMPLEX * const bound_states,
		fnft_manakovv_opts_t * const opts)
{
	REAL degree1step = 0.0 , map_coeff = 2.0;
	UINT K, upsampling_factor, i, j, D_given;
	REAL bounding_box[4] = { NAN };
	COMPLEX * buffer = NULL;
	INT ret_code = SUCCESS;
	manakov_discretization_t discretization;
/*
	degree1step = manakov_discretization_degree(opts->discretization);
	// degree1step == 0 here indicates a valid slow method. Incorrect
	// discretizations should have been caught earlier.
	upsampling_factor = manakov_discretization_upsampling_factor(opts->discretization);
	if (upsampling_factor == 0) {
		ret_code = E_INVALID_ARGUMENT(opts->discretization);
		goto leave_fun;
	}
	if (degree1step != 0)
		map_coeff = 2/(degree1step);
	D_given = D/upsampling_factor;
*/
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




/* Code for disc. spec.
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
		REAL const * const T, REAL const eps_t)
{
	// The nonlinear Parseval relation tells us that the squared L2 norm of
	// q(t) is >= 4*(sum of the imaginary parts of the bound states). Thus,
	// any bound state with an imaginary part greater than four times the
	// squared L2 norm of q(t) can be removed. A factor of 1.5 has been
	// added to account for numerical discrepancies when computing the norm
	// numerically (e.g., truncation errors or large step sizes).
	return 1.5 * 0.25 * misc_l2norm2(D, q, T[0]-eps_t/2, T[1]+eps_t/2);
}


// Auxiliary function: Computes the bound states.
static inline INT nsev_compute_boundstates(
		const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
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

	// Set-up bounding_box based on choice of filtering
	if (opts->bound_state_filtering == nsev_bsfilt_BASIC) {
		bounding_box[0] = -INFINITY;
		bounding_box[1] = INFINITY;
		bounding_box[2] = 0.0;
		bounding_box[3] = INFINITY;
	}else if (opts->bound_state_filtering == nsev_bsfilt_FULL) {
		bounding_box[1] = re_bound(eps_t, map_coeff);
		bounding_box[0] = -bounding_box[1];
		bounding_box[2] = 0;
		// This step is required as q contains scaled values on a
		// non-equispaced grid
		switch(opts->discretization) {
			case fnft_nse_discretization_ES4:
			case fnft_nse_discretization_TES4:
			{
				// For these discretizations we need to skip the time-derivative samples for determining the L2-norm of q
				COMPLEX * const q_tmp = malloc(D_given * sizeof(COMPLEX));
				if (q_tmp == NULL) {
					ret_code = E_NOMEM;
					CHECK_RETCODE(ret_code, leave_fun);
				}
				for (i = 0, j = 0; i < D_given; i++, j+=3)
					q_tmp[i] = q[j];
				bounding_box[3] = im_bound(D_given, q_tmp, T, eps_t);
				free(q_tmp);
			}
				break;
			default:
				bounding_box[3] = im_bound(D, q, T, eps_t);
		}
	}else{
		bounding_box[0] = -INFINITY;
		bounding_box[1] = INFINITY;
		bounding_box[2] = -INFINITY;
		bounding_box[3] = INFINITY;
	}

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
					discretization, opts->niter, bounding_box);
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

		default:

			return E_INVALID_ARGUMENT(opts->bound_state_localization);
	}

	// Filter bound states
	if (opts->bound_state_filtering != nsev_bsfilt_NONE) {

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
*/

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
	//    misc_print_buf(deg*D,transfer_matrix,"transfer_matrix_in_contspec");
	COMPLEX* H11_vals = NULL, * H21_vals = NULL, * H31_vals = NULL;
	COMPLEX A, V;
	REAL scale;
	REAL phase_factor_rho, phase_factor_a, phase_factor_b1, phase_factor_b2;
	INT ret_code = SUCCESS;
	UINT i, offset = 0, upsampling_factor, D_given;
	COMPLEX* scatter_coeffs = NULL, * xi = NULL;
	// defining different pointers to the different entries of transfer_matrix
	COMPLEX* tm11_ptr, * tm12_ptr, * tm13_ptr, * tm21_ptr, * tm22_ptr, * tm23_ptr, * tm31_ptr, * tm32_ptr, * tm33_ptr;

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

		printf("Inputs to manakov_scatter_matrix\n");
		printf("D=%d\n", D);
		misc_print_buf(D, q1, "q1");
		misc_print_buf(D, q2, "q2");
		printf("eps_t=%f\n", eps_t);
		printf("M=%d\n", M);
		misc_print_buf(M, xi, "xi");
		printf("kappa=%d\n", kappa);
		ret_code = manakov_scatter_matrix(D, q1, q2, eps_t, M, xi, kappa,
			scatter_coeffs, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);
		printf("scatter_coeffs after manakov_scatter_matrix\n");
		misc_print_buf(9 * M, scatter_coeffs, "scatter_coeffs");

		// This is necessary because nse_scatter_matrix to ensure
		// boundary conditions can be applied using common code for slow
		// methods and polynomial transfer matrix based methods.
		for (i = 0; i < M; i++) {
			H11_vals[i] = scatter_coeffs[i * 9];
			H21_vals[i] = scatter_coeffs[i * 9 + 3];
			H31_vals[i] = scatter_coeffs[i * 9 + 6];
		}
		// Check if we have the right values stored in Hij:
		misc_print_buf(M, H11_vals, "H11");
		misc_print_buf(M, H21_vals, "H21");
		misc_print_buf(M, H31_vals, "H31");

/*		if (opts->discretization == manakov_discretization_CF4_2){
		for (i = 0; i < M; i++)
		xi[i] = (XI[0] + eps_xi * i)/2;
		}*/

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


	// TODO: remove this
	// printing Hij_vals to compare with matlab
	/*
	for (UINT i=0; i<M; i++){
		printf("H11[%d] = %f + i%f\n",i, creal(H11_vals[i]), cimag(H11_vals[i]));
	}
	for (UINT i=0; i<M; i++){
		printf("H21[%d] = %f + i%f\n",i, creal(H21_vals[i]), cimag(H21_vals[i]));
	}
	for (UINT i=0; i<M; i++){
		printf("H31[%d] = %f + i%f\n",i, creal(H31_vals[i]), cimag(H31_vals[i]));
	}

	//misc_print_buf(3*M, contspec+2*M, "ab_num");
	misc_print_buf(M, H11_vals, "H11");
	misc_print_buf(M, H21_vals, "H21");
	misc_print_buf(M, H31_vals, "H31");
	*/

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

		ret_code = manakov_discretization_phase_factor_b(eps_t, D_given, T, &phase_factor_b1, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);

		// TODO: remove this second call, phase factor b2= phase factor b1
		ret_code = manakov_discretization_phase_factor_b(eps_t, D_given, T, &phase_factor_b2, opts->discretization);
		CHECK_RETCODE(ret_code, leave_fun);
		printf("phase_factor a = %f\n", phase_factor_a);
		printf("phase_factor b = %f\n", phase_factor_b1);

		for (i = 0; i < M; i++) {
			result[offset + i] = H11_vals[i] * scale * CEXP(I * xi[i] * phase_factor_a);
			result[offset + M + i] = H21_vals[i] * scale * CEXP(I * xi[i] * phase_factor_b1);
			result[offset + 2 * M + i] = H31_vals[i] * scale * CEXP(I * xi[i] * phase_factor_b2);
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


/* Code for slow methods
// Auxiliary function: Computes the norming constants and/or residues
// using slow scattering schemes
static inline INT nsev_compute_normconsts_or_residues(
		const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
		REAL const * const T,
		const UINT K,
		COMPLEX * const bound_states,
		COMPLEX * const normconsts_or_residues,
		fnft_nsev_opts_t const * const opts)
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
*/

/* Code for disc. spec.
// Auxiliary function: Refines the bound-states using Newtons method
static inline INT nsev_refine_bound_states_newton(
		const UINT D,
		COMPLEX const * const q,
		COMPLEX const * const r,
		REAL const * const T,
		const UINT K,
		COMPLEX * const bound_states,
		nse_discretization_t const discretization,
		const UINT niter,
		REAL const * const bounding_box)
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
			ret_code = nse_scatter_bound_states(D, q, r, T, 1,
					bound_states + i, &a_val, &aprime_val, &b_val, discretization, 1);
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
*/