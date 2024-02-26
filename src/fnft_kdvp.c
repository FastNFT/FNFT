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
 * Sander Wahls (KIT) 2023.
 */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft_kdvp.h"
#include "fnft__kdv_scatter.h"

static fnft_kdvp_opts_t default_opts = {
    .mainspec_type = kdvp_mstype_EDGEPOINTS_AND_SIGNS,
    .normalization_flag = 1,
    .keep_degenerate_flag = 0,
    .discretization = kdv_discretization_BO,
    .grid_spacing = 0,
    .niter = 100,
    .tol = -1.0 // auto
};

/**
 * Creates a new options variable for fnft_kdvp with default settings.
 * See the header file for a detailed description.
 */
fnft_kdvp_opts_t fnft_kdvp_default_opts()
{
    return default_opts;
}

static inline COMPLEX E_to_lambda(REAL E)
{
    if (E < 0)
        return I*SQRT(-E);
    else
        return SQRT(E);
}

static inline REAL get_sign(REAL x)
{
    if (x > 0)
        return 1;
    else if (x < 0)
        return -1;
    return 0;
}

/**
 * Represents a real number A in the form
 *
 *  A = sign*2^log2A,
 *
 * where sign=-1,0,1 and log2A is real (or possibly NaN if X=0).
 */
typedef struct 
{
    REAL sign;
    REAL log2A;
} LogNumber;

/**
 * Check if A=num.sign*2^num.log2A satisfies A<=1.
 */
static inline int lt1(LogNumber num)
{
    if (num.sign == -1 || num.sign == 0)
        return 1;
    // now we have sign == 1
    return num.log2A <= 0;
}

/**
 * Same, but for <=-1.
 */
static inline int ltm1(LogNumber num)
{
    return num.sign == -1 && num.log2A >= 0;
}

/**
 * Same, but for A*2^W>=1.
 */
static inline int gt1(LogNumber num)
{
    return num.sign == 1 && num.log2A >= 0;
}

/**
 * Same, but for >=-1.
*/
static inline int gtm1(LogNumber num)
{
    if (num.sign == 0 || num.sign == 1)
        return 1;
    // now we have sign == -1
    return num.log2A <= 0;
}

/**
 * Computes the quantities Delta(E) and alpha_21(E) as given e.g. the Eqs. 3.11 and 3.12 of Osborne,
 * Math. Comput. Simul. 37 (1994), 431-450, https://doi.org/10.1016/0378-4754(94)00029-8
 */
static inline INT compute_DEL_and_al21(const UINT D, COMPLEX * const q, COMPLEX const * const r,
        const REAL eps_t, const REAL Ei, REAL * const DEL, LogNumber * const DEL_LN, REAL * al21,
        REAL * al21n, INT * W_ptr, fnft_kdvp_opts_t * opts_ptr)
{
    COMPLEX scatter_matrix[4];
    COMPLEX lambda;
    INT ret_code = SUCCESS;

    lambda = E_to_lambda(Ei);
    ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
        &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 0/*derivative_flag*/);
    CHECK_RETCODE(ret_code, leave_fun);

    /* Compute the Floquet discriminant Delta(E) */

    REAL trace = 0.5*CREAL(scatter_matrix[0] + scatter_matrix[3]);
    if (opts_ptr->normalization_flag) {
        *DEL = POW(2, *W_ptr)*trace;
        *DEL_LN = (LogNumber){.sign = get_sign(trace), .log2A = LOG2(FABS(trace)) + *W_ptr};
    } else {
        *DEL = trace;
        *DEL_LN = (LogNumber){.sign = get_sign(trace), .log2A = LOG2(FABS(trace))};
    }

    /* Compute al21(E) */

    const COMPLEX sumof = scatter_matrix[0] - scatter_matrix[1] + scatter_matrix[2] - scatter_matrix[3];
    // The straight-forward implementation would be
    //   *al21n = CREAL(-I/(2*lambda)*sumof);
    //   *al21 = POW(2, *W_ptr)**al21n;
    // (Note that the signs in sumof differ from the ones in the paper of Osborne. The
    // reason should be that we use the slightly different S-basis instead of the basis
    // of Osborne. See Appendix C-C of Prins and Wahls, IEEE Access, Jul. 2019.)
    // The variant below pulls the I/lambda term out of al21n.
    if (CIMAG(lambda) == 0) // Note: by construction, for some x>=0, either lambda=x or lambda=I*x
        *al21n = CIMAG(sumof);
    else
        *al21n = -CREAL(sumof);
    if (opts_ptr->normalization_flag)
        *al21 = POW(2, *W_ptr)* *al21n;
    else
        *al21 = *al21n;
    if (CIMAG(lambda) == 0)
        *al21 /= 2*CREAL(lambda);
    else
        *al21 /= 2*CIMAG(lambda);
leave_fun:
    return ret_code;
}

/**
 * Computes the sheet index as given e.g. in Eq. 3.13 of Osborne, Math. Comput. Simul. 37 (1994),
 * 431-450, https://doi.org/10.1016/0378-4754(94)00029-8
 */
static inline INT compute_sheet_index(const UINT D, COMPLEX * const q, COMPLEX const * const r,
        const REAL eps_t, const REAL Ei, REAL * const sheet_index, INT * W_ptr, fnft_kdvp_opts_t * opts_ptr)
{
    COMPLEX scatter_matrix[4];
    COMPLEX lambda;
    INT ret_code = SUCCESS;

    lambda = E_to_lambda(Ei);
    ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
        &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 0/*derivative_flag*/);
    CHECK_RETCODE(ret_code, leave_fun);

    *sheet_index = get_sign(scatter_matrix[1] + scatter_matrix[2]);

leave_fun:
    return ret_code;
}

static inline INT refine_mainspec(const UINT D, COMPLEX * const q, COMPLEX const * const r,
        const REAL eps_t, REAL E_left, REAL E_right, REAL * const E_ptr, const REAL s, INT * W_ptr, fnft_kdvp_opts_t * opts_ptr)
{
    COMPLEX scatter_matrix[8];
    INT ret_code = SUCCESS;
    const UINT niter = opts_ptr->niter;
    const REAL tol = opts_ptr->tol<0 ? 100*EPSILON : opts_ptr->tol;

    const REAL E0 = *E_ptr;
    REAL E = E0;
    for (UINT i=0; i<niter; i++) {
        COMPLEX lambda = E_to_lambda(E);
        if (lambda == 0) // Can't scatter at zero => stop
            break;

        ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
            &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 1/*derivative_flag*/);
        CHECK_RETCODE(ret_code, leave_fun);

        const INT W = W_ptr == NULL ? 0 : *W_ptr;
        const COMPLEX DELn = 0.5*(scatter_matrix[0] + scatter_matrix[3]);
        const COMPLEX DELn_prime = 0.5*(scatter_matrix[4] + scatter_matrix[7]);
        const COMPLEX DELn_ms = DELn - POW(2, -W)*s;
        const COMPLEX incr = DELn_ms/DELn_prime;

        if (incr != incr) // Can't compute increment => stop
            break;
        lambda -= incr;
        E = CREAL(lambda*lambda);

        if (!(E>=E_left) || !(E<=E_right)) { // Newton's method diverges => stop
            E = E0;
            break;
        }
        if (CABS(incr)<=tol || CABS(DELn_ms)<=tol) // We converged => stop
            break;
    }
    *E_ptr = E;
leave_fun:
    return ret_code;
}

static inline INT refine_auxspec(const UINT D, COMPLEX * const q, COMPLEX const * const r,
        const REAL eps_t, REAL E_left, REAL E_right, REAL * const mu_ptr, INT * W_ptr, fnft_kdvp_opts_t * opts_ptr)
{
    COMPLEX scatter_matrix[8];
    INT ret_code = SUCCESS;
    const UINT niter = opts_ptr->niter;
    const REAL tol = opts_ptr->tol<0 ? 100*EPSILON : opts_ptr->tol;

    const REAL mu0 = *mu_ptr;
    REAL mu = mu0;
    for (UINT i=0; i<niter; i++) {
        COMPLEX lambda = E_to_lambda(mu);
        if (lambda == 0) // Can't scatter at zero => stop
            break;

        ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
            &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 1/*derivative_flag*/);
        CHECK_RETCODE(ret_code, leave_fun);

        const COMPLEX sumof = scatter_matrix[0] - scatter_matrix[1] + scatter_matrix[2] - scatter_matrix[3];
        const COMPLEX sumof_prime = scatter_matrix[4] - scatter_matrix[5] + scatter_matrix[6] - scatter_matrix[7];
        //const COMPLEX incr = (-I/(2*lambda))*sumof/(I/(2*lambda*lambda)*sumof + -I/(2*lambda)*sumof_prime);
        const COMPLEX incr = -sumof/(sumof/lambda - sumof_prime); // equivalent to the line above

        if (incr != incr) // Can't compute increment => stop
            break;
        lambda -= incr;        
        mu = CREAL(lambda*lambda);

        if (!(mu>=E_left) || !(mu<=E_right)) { // Newton's method seems to diverge => stop
            mu = mu0;
            break;
        }
        if (CABS(incr)<=tol || CABS(sumof)<=tol) // We converged => stop
            break;
    }
    *mu_ptr = mu;
leave_fun:
    return ret_code;
}

/**
 * Nonlinear Fourier transform for the Korteweg-de Vries
 * equation with periodic boundary conditions.
 */
INT fnft_kdvp(  const UINT D,
                COMPLEX * const q,
                REAL const * const T,
                REAL * const E,
                UINT * const K_ptr,
                REAL * const main_spec, 
                UINT * const M_ptr,
                REAL * const aux_spec,
                REAL * const sheet_indices,
                fnft_kdvp_opts_t * opts_ptr)
{
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (T == NULL || T[0] >= T[1])
        return E_INVALID_ARGUMENT(T);
    if (E == NULL || E[0] >= E[1])
        return E_INVALID_ARGUMENT(E);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (main_spec == NULL)
        return E_INVALID_ARGUMENT(main_spec);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);   
    if (aux_spec == NULL)
        return E_INVALID_ARGUMENT(aux_spec); 
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->mainspec_type != kdvp_mstype_FLOQUET && opts_ptr->grid_spacing <= 0)
        return E_INVALID_ARGUMENT(opts_ptr->grid_spacing);
    if (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET && *K_ptr != *M_ptr)
        return E_INVALID_ARGUMENT(*K_ptr != *M_ptr for mainspec_type FLOQUET);
    if (opts_ptr->discretization != kdv_discretization_BO && opts_ptr->discretization != kdv_discretization_BO_VANILLA)
        return E_NOT_YET_IMPLEMENTED(opts_ptr->discretization, Use the BO or BO_VANILLA discretizations);

    const INT aux_spec_flag = (aux_spec == NULL || M_ptr == NULL || *M_ptr == 0) ? 0 : 1;
    const INT sheet_indices_flag = (sheet_indices == NULL || !aux_spec_flag) ? 0 : 1;

    COMPLEX * r = NULL;
    REAL Ei = 0, Ei_prev = 0;
    REAL DEL = 0, DEL_prev = 0;
    REAL al21 = 0, al21_prev = 0;
    LogNumber DEL_prev_LN = {0}, DEL_LN = {0};
    INT W = 0;
    REAL al21n = 0, al21n_prev = 0;
    const INT normalization_flag = opts_ptr->normalization_flag;
    const INT keep_degenerate_flag = opts_ptr->keep_degenerate_flag;
    INT * const W_ptr = (normalization_flag) ? &W : NULL;
    UINT K = 0, M = 0, N_bands = 0;
    UINT i = 0;
    INT ret_code = SUCCESS;

    r = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(r, ret_code, leave_fun);

    for (i=0; i<D; i++)
        r[i] = -1;

    const UINT L = (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET) ? *K_ptr : CEIL((E[1] - E[0]) / opts_ptr->grid_spacing);
    const REAL eps_t = (T[1] - T[0])/D;
    const REAL eps_E = (E[1] - E[0])/(L - 1);

    // Note that kdv_scatter_matrix currently fails at E=0 because a transformation matrix becomes singular at this point.
    // We therefore just skip a grid point in the E-domain if it is close to zero, which effectivly increases the
    // E-step size by a factor of two for a moment. For the FLOQUET mode, we just put NANs (not a problem for plotting).

    const REAL close_to_zero = sqrt(EPSILON) < 0.49*eps_E ? sqrt(EPSILON) : 0.49*eps_E; // threshold to identify when we
                                                                                        // are problematically close to E=0
    
    if (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET) { // implies auxspec_type == ALPHA21
                                                          
        for (i=0; i<L; i++) {
            Ei = E[0] + i*eps_E;
            if (FABS(Ei) >= close_to_zero) {
                ret_code = compute_DEL_and_al21(D, q, r, eps_t, Ei, &main_spec[i], &DEL_LN, &aux_spec[i], &al21n, W_ptr, opts_ptr);
                CHECK_RETCODE(ret_code, leave_fun);
            } else {
                main_spec[i] = NAN;
                aux_spec[i] = NAN;
            }
        }

    } else {

        UINT i0 = 0;
        Ei_prev = E[0];
        if (FABS(Ei_prev) < close_to_zero) {
            Ei_prev = E[0] + eps_E;
            i0 = 1;
        }
        ret_code = compute_DEL_and_al21(D, q, r, eps_t, Ei_prev, &DEL_prev, &DEL_prev_LN, &al21_prev, &al21n_prev, W_ptr, opts_ptr);
        CHECK_RETCODE(ret_code, leave_fun);

        for (i=i0; i<L; i++) {
            Ei = E[0] + i*eps_E;
            if (FABS(Ei) < close_to_zero)
                continue;

            ret_code = compute_DEL_and_al21(D, q, r, eps_t, Ei, &DEL, &DEL_LN, &al21, &al21n, W_ptr, opts_ptr);
            CHECK_RETCODE(ret_code, leave_fun);

            /* Auxiliary spectrum */

            // Note that the normalized al21 always has a zero at E=0 because
            // it does not have the I/k term (it's in al21). However, since we are
            // skipping E=0 anyways, this does not create a spurious aux spectrum point.
            const INT found_auxspec_flag = aux_spec_flag && al21n_prev*al21n < 0;
            if (found_auxspec_flag) { // zero-crossing detected
                if (M>=*M_ptr) {
                    ret_code = E_OTHER("Found more than *M_ptr auxiliary spectrum points found. Increase *M_ptr and try again.")
                    goto leave_fun;
                }

                // regula falsi for al21
                aux_spec[M] = (Ei_prev*al21 - Ei*al21_prev)/(al21 - al21_prev);
               
                ret_code = refine_auxspec(D, q, r, eps_t, Ei_prev, Ei, &aux_spec[M], W_ptr, opts_ptr);
                CHECK_RETCODE(ret_code, leave_fun);

                if (sheet_indices_flag) {
                    ret_code = compute_sheet_index(D, q, r, eps_t, aux_spec[M], &sheet_indices[M], W_ptr, opts_ptr);
                    CHECK_RETCODE(ret_code, leave_fun);
                }


                M++;
            }
 
            /* Main spectrum */

            // Check if Delta(Ei) crossed +/-1.
            REAL s = 0;
            if (!normalization_flag) { // compare in the linear domain
                if (DEL_prev<=-1 && DEL>=1)
                    s = 2; // two edge points with delta increasing
                else if (DEL_prev>=1 && DEL<=-1)
                    s = -2; // two edge points with delta decreasing
                else if ((DEL_prev<=1 && DEL>=1) || (DEL_prev>=1 && DEL<=1))
                    s = 1; // one edge point s.t. delta=1
                else if ((DEL_prev>=-1 && DEL<=-1) || (DEL_prev<=-1 && DEL>=-1))
                    s = -1; // one edge point s.t. delta=-1
            } else { // compare in the logarithmic domain
                if (ltm1(DEL_prev_LN) && gt1(DEL_LN))
                    s = 2; // two edge points with delta increasing
                else if (gt1(DEL_prev_LN) && ltm1(DEL_LN))              
                    s = -2; // two edge points with delta decreasing
                else if ((lt1(DEL_prev_LN) && gt1(DEL_LN)) || (gt1(DEL_prev_LN) && lt1(DEL_LN)))
                    s = 1; // one edge point s.t. delta=1
                else if ((gtm1(DEL_prev_LN) && ltm1(DEL_LN)) || (ltm1(DEL_prev_LN) && gtm1(DEL_LN)))
                    s = -1; // one edge point s.t. delta=-1
            }

            // Process zero-crossings
            if (s == 2 || s == -2) { // two edge points between two grid points
                if (K+1>=*K_ptr) {
                    ret_code = E_OTHER("Found more than *K_ptr main spectrum points. Increase *K_ptr and try again.")
                    goto leave_fun;
                }

                // regula falsi for DEL-j with j=-1,1; the "-j" terms in the denominator cancel
                for (INT j=-1; j<=1; j+=2) {
                    main_spec[2*K] = (Ei_prev*(DEL-j) - Ei*(DEL_prev-j))/(DEL - DEL_prev);
                    main_spec[2*K+1] = j;

                    ret_code = refine_mainspec(D, q, r, eps_t, Ei_prev, Ei, &main_spec[2*K], j, W_ptr, opts_ptr);
                    CHECK_RETCODE(ret_code, leave_fun);

                    K++;
                }

                // fix ordering if necessary
                if (s == -2) {
                    REAL tmp = main_spec[2*(K-2)];
                    main_spec[2*(K-2)] = main_spec[2*(K-1)];
                    main_spec[2*(K-1)] = tmp;
                    tmp = main_spec[2*(K-2)+1];
                    main_spec[2*(K-2)+1] = main_spec[2*(K-1)+1];
                    main_spec[2*(K-1)+1] = tmp;
                }

            } else if (s == 1 || s == -1) { // one edge point between two grid points
                if (K>=*K_ptr) {
                    ret_code = E_OTHER("More than *K_ptr initial guesses for bound states found. Increase *K_ptr and try again.")
                    goto leave_fun;
                }

                // regula falsi for DEL-s, the "-s" terms in the denominator cancel
                main_spec[2*K] = (Ei_prev*(DEL-s) - Ei*(DEL_prev-s))/(DEL - DEL_prev);
                main_spec[2*K+1] = s;

                ret_code = refine_mainspec(D, q, r, eps_t, Ei_prev, Ei, &main_spec[2*K], s, W_ptr, opts_ptr);
                CHECK_RETCODE(ret_code, leave_fun);

                K++;
            }

            // If no zero-crossings for Delta(E) +/- 1 were found, there could still be a pair of co-located
            // degenerate eigenvalues at which Delta(E) just touches, but does not cross, +/- 1. Such cases
            // are indicated by a co-located auxiliary spectrum point at seems to be not in an open band.
            // This is what we check next.
            
            INT in_open_band_flag;
            if (!normalization_flag)
                in_open_band_flag = !((DEL_prev >= -1 && DEL_prev <= 1) && (DEL >= -1 && DEL <= 1));
            else
                in_open_band_flag = !((gtm1(DEL_prev_LN) && lt1(DEL_prev_LN)) && (gtm1(DEL_LN) && lt1(DEL_LN)));

            if (!in_open_band_flag && found_auxspec_flag) { // looks as if we missed a pair of degenerate
                                                            // main spectrum points
                if (keep_degenerate_flag) {
                    if (K+1>=*K_ptr) {
                        ret_code = E_OTHER("More than *K_ptr initial guesses for bound states found. Increase *K_ptr and try again.")
                        goto leave_fun;
                    }
                
                    if (!normalization_flag)
                        s = get_sign(DEL);
                    else
                        s = DEL_LN.sign;
                
                    main_spec[2*K] = aux_spec[M-1];
                    main_spec[2*K+1] = s;
                    K++;

                    main_spec[2*K] = aux_spec[M-1];
                    main_spec[2*K+1] = s;
                    K++;
                } else {
                    M--;
                }
            }

            Ei_prev = Ei;
            DEL_prev = DEL;
            DEL_prev_LN = DEL_LN;
            al21_prev = al21;
            al21n_prev = al21n;
        }
        *M_ptr = M;

        /* Change representation of the main spectrum if requested by the user. */

        switch (opts_ptr->mainspec_type) {
        case kdvp_mstype_EDGEPOINTS_AND_SIGNS:
            *K_ptr = K; // no changes, this is what we computed above
            break;

        case kdvp_mstype_OPENBANDS:
            for (i=1; i<K; i++) {
                if (main_spec[2*i-1] == main_spec[2*i+1]) { // same signs s => open band
                    const REAL left_edge = main_spec[2*i-2];
                    const REAL right_edge = main_spec[2*i];
                    main_spec[2*N_bands] = left_edge;
                    main_spec[2*N_bands+1] = right_edge;
                    N_bands++;
                    i++; // required when a tiny band on the opposite side is missed
                }
            }

            *K_ptr = N_bands;
            break;

        case kdvp_mstype_AMPLITUDES_MODULI_FREQS:
            if (K < 3) { // nothing to do
                *K_ptr = 0;
                break;
            }

            // replace main spectrum points E_i with amplitudes and signs with moduli
            UINT cnt = 0;
            UINT i_ref = 0;
            REAL E_ref = FNFT_INF; // reference level
            INT in_radiation = 0; // flag to signal that we reached the radiation
                                  // part of the nonlinear spectrum
            for (i=1; i<K-1; i+=2) {
                const REAL E_2i = main_spec[2*i];
                const REAL E_2ip1 = main_spec[2*(i+1)];
                const REAL E_2im1 = main_spec[2*(i-1)];
                const REAL modulus = (E_2ip1 - E_2i)/(E_2ip1 - E_2im1);
                main_spec[3*cnt+1] = modulus;
                if (modulus >= 0.99 && !in_radiation) { // soliton
                    i_ref = cnt;
                    E_ref = E_2ip1;
                    main_spec[3*cnt] = E_2i; // reference level is added later
                } else { // radiation
                    main_spec[3*cnt] = 0.5*(E_2ip1 - E_2i);
                    in_radiation = 1;
                }
                main_spec[3*cnt+2] = 0.5*(E_2ip1 + E_2i); // needed for nonlinear frequencies (below)
                cnt++;
            }

            // finalize the soliton amplitudes
            if (E_ref < FNFT_INF) {
                for (i=0; i<=i_ref; i++)
                    main_spec[3*i] = 2*(E_ref - main_spec[3*i]);
            }

            // finalize the nonlinear frequencies, see Eqs. A.2 and A.3 in Bruehl et al,
            // Wave Motion 111 (2022), https://doi.org/10.1016/j.wavemoti.2022.102905
            for (i=0; i<cnt; i++) {
                const REAL E_bar = main_spec[3*i+2] - E_ref;
                if (E_bar < 0)
                    main_spec[3*i+2] = -SQRT(-E_bar);
                else
                    main_spec[3*i+2] = SQRT(E_bar);
            }

            *K_ptr = cnt;
            break;

        default:
            ret_code = E_INVALID_ARGUMENT(opts_ptr->mainspec_type);
            goto leave_fun;
        }
    }

leave_fun:
    free(r);
    return ret_code;
}
