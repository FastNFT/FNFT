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
    .discretization = kdv_discretization_BO,
    .grid_spacing = 0
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
 * Nonlinear Fourier transform for the Korteweg-de Vries
 * equation with periodic boundary conditions.
 */
INT fnft_kdvp(  const UINT D,
                COMPLEX const * const q,
                REAL const * const T,
                REAL const * const E,
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
    if (E == NULL || E[0] >= E[1] || E[0] == 0 || E[1] == 0)
        return E_INVALID_ARGUMENT(E);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (main_spec == NULL)
        return E_INVALID_ARGUMENT(main_spec);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);   
    if (aux_spec == NULL)
        return E_INVALID_ARGUMENT(aux_spec); 
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, Pass sheet_indices="NULL");
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->mainspec_type != kdvp_mstype_FLOQUET && opts_ptr->grid_spacing <= 0)
        return E_INVALID_ARGUMENT(opts_ptr->grid_spacing);
    if (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET && *K_ptr != *M_ptr)
        return E_INVALID_ARGUMENT(*K_ptr != *M_ptr for mainspec_type FLOQUET);
    if (opts_ptr->discretization != kdv_discretization_BO)
        return E_NOT_YET_IMPLEMENTED(opts_ptr->discretization, Use the BO discretization);

    COMPLEX scatter_matrix[4] = {0};
    COMPLEX * r = NULL;
    COMPLEX lambda = 0;
    REAL Ei = 0, Ei_prev = 0;
    REAL DEL = 0, DEL_prev = 0;
    REAL al21 = 0, al21_prev = 0;
    LogNumber DEL_prev_LN, DEL_LN;
    INT W = 0, W_prev = 0;
    REAL al21n = 0, al21n_prev = 0;
    const INT normalization_flag = opts_ptr->normalization_flag;
    INT * const W_ptr = (normalization_flag) ? &W : NULL;
    UINT N_bands = 0;
    UINT i = 0;
    INT ret_code = SUCCESS;

    r = malloc(D * sizeof(COMPLEX));
    CHECK_NOMEM(r, ret_code, leave_fun);

    for (i=0; i<D; i++)
        r[i] = -1;

    const UINT L = (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET) ? *K_ptr : CEIL((E[1] - E[0]) / opts_ptr->grid_spacing);
    const REAL eps_t = (T[1] - T[0])/D;
    const REAL eps_E = (E[1] - E[0])/(L - 1);

    if (opts_ptr->mainspec_type == kdvp_mstype_FLOQUET) { // implies auxspec_type == ALPHA21
        for (i=0; i<L; i++) {
            Ei = E[0] + i*eps_E;
            lambda = E_to_lambda(Ei);
            if (lambda == 0) {
                main_spec[i] = NAN; // currently kdv_scatter_matrix fails for lambda==0
                continue;
            }
            ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
                    &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 0/*derivative_flag*/);
            CHECK_RETCODE(ret_code, leave_fun);

            const REAL scl = POW(2, W);
            DEL = 0.5*scl*CREAL(scatter_matrix[0] + scatter_matrix[3]);
            main_spec[i] = DEL;

            //al21n = CREAL(-I/(2*lambda)*(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]));
            if (Ei < 0)
                al21n = CREAL(-(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]));
            else
                al21n = CIMAG(-(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]));
            if (Ei == 0)
                misc_print_buf(4, scatter_matrix, "scatter_matrix");
            al21 = scl*al21n;
            aux_spec[i] = al21;
        }
    } else {

        // Commmon step for all other mainspec_types: Find the main spectrum, i.e., the +/-1 crossings
        // E_i of the Floquet discriminant Delta(E). For efficiency reasons, the auxiliary spectrum is
        // computed here as well.
        
        // Using normalization is a bit more complicated here because we only have Delta(E)=A*2^W in a
        // normalized form, not Delta(E)+1 and Delta(E)-1. The approach below is to utilize a logarithmic
        // number system, in which we represent Delta(E) = sign*2^log2(A)*2^W = sign*2^(log2(A)+w), for
        // the comparisons Delta(E) >= +/-1 and Delta(E) <= +/-1.

        Ei_prev = E[0];
        lambda = E_to_lambda(Ei_prev);
        ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
            &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 0/*derivative_flag*/);
        CHECK_RETCODE(ret_code, leave_fun);

        REAL tmp = 0.5*CREAL(scatter_matrix[0] + scatter_matrix[3]);
        DEL_prev = POW(2,W)*tmp;    // required because kdv_scatter_matrix normalizes
                                    // the scatter matrix (unless normalization_flag is false)
        if (normalization_flag)
            DEL_prev_LN = (LogNumber){.sign = get_sign(tmp), .log2A = LOG2(FABS(tmp)) + W};

        al21n_prev = CREAL(-I/(2*lambda)*(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]));
        al21_prev = POW(2,W)*al21n_prev;

        UINT K = 0;
        UINT M = 0;
        for (i=1; i<L; i++) {
            Ei = E[0] + i*eps_E;
            lambda = E_to_lambda(Ei);
            // Catch that KdV scattering currently fails at lambda=0 because the AKNS change
            // of basis matrix becomes singular at this point.
            if (lambda == 0)
                continue;
            ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
                    &lambda, scatter_matrix, W_ptr, opts_ptr->discretization, 0/*derivative_flag*/);
            CHECK_RETCODE(ret_code, leave_fun);

            /* Main spectrum */
    
            REAL tmp = 0.5*CREAL(scatter_matrix[0] + scatter_matrix[3]);
            DEL = POW(2,W)*tmp;
            if (normalization_flag)
                DEL_LN = (LogNumber){.sign = get_sign(tmp), .log2A = LOG2(FABS(tmp)) + W};

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

            if (s == 2 || s == -2) { // two edge points between two grid points
                if (K+1>=*K_ptr) {
                    ret_code = E_OTHER("Found more than *K_ptr main spectrum points. Increase *K_ptr and try again.")
                    goto leave_fun;
                }

                // regula falsi for DEL-j with j=-1,1; the "-j" terms in the denominator cancel
                for (INT j=-1; j<=1; j+=2) {
                    main_spec[2*K] = (Ei_prev*(DEL-j) - Ei*(DEL_prev-j))/(DEL - DEL_prev);
                    main_spec[2*K+1] = j;
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
                K++;
            }

            /* Auxiliary spectrum */

            // The naive implementation would just be to use
            //   al21n = CREAL(-I/(2*lambda)*(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]));
            //   al21 = POW(2, W)*al21n;
            if (CIMAG(lambda) == 0) // Note: by construction, for some x>=0, either lambda=x or lambda=I*x
                al21n = CIMAG(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]);
            else
                al21n = -CREAL(scatter_matrix[0] + scatter_matrix[1] - scatter_matrix[2] - scatter_matrix[3]);
            al21 = POW(2,W)*al21n;
            if (CIMAG(lambda) == 0)
                al21 /= 2*CREAL(lambda);
            else
                al21 /= 2*CIMAG(lambda);
            if (al21n_prev*al21n<0) { // zero-crossing detected
                if (M>=*M_ptr) {
                    ret_code = E_OTHER("Found more than *M_ptr auxiliary spectrum points found. Increase *M_ptr and try again.")
                    goto leave_fun;
                }
                // regula falsi for al21
                aux_spec[M] = (Ei_prev*al21 - Ei*al21_prev)/(al21 - al21_prev);
                M++;
            }
 
            Ei_prev = Ei;
            DEL_prev = DEL;
            DEL_prev_LN = DEL_LN;
            W_prev = W;
            al21_prev = al21;
            al21n_prev = al21n;
        }
        *M_ptr = M;

        // We so far have the usual representation of the main spectrum (i.e., the Ei). If the user requested
        // another representation, we convert it now.
        switch (opts_ptr->mainspec_type) {
        case kdvp_mstype_EDGEPOINTS_AND_SIGNS:
            *K_ptr = K;
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
            REAL E_ref = 0; // reference level
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
