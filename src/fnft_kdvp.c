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

/**
 * Fast nonlinear Fourier transform for the nonlinear Schroedinger
 * equation with vanishing boundary conditions.
 */
INT fnft_kdvp(  const UINT D,
                COMPLEX const * const q,
                REAL const * const T,
                REAL const * const E,
                UINT * const K_ptr,
                REAL * const main_spec, 
                UINT * const L_ptr,
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
    if (aux_spec != NULL)
        return E_NOT_YET_IMPLEMENTED(aux_spec, Pass aux_spec="NULL");
    if (sheet_indices != NULL)
        return E_NOT_YET_IMPLEMENTED(sheet_indices, Pass sheet_indices="NULL");
    if (opts_ptr == NULL)
        opts_ptr = &default_opts;
    if (opts_ptr->mainspec_type != kdvp_mstype_FLOQUET && opts_ptr->grid_spacing <= 0)
        return E_INVALID_ARGUMENT(opts_ptr->grid_spacing);
    if (opts_ptr->discretization != kdv_discretization_BO)
        return E_NOT_YET_IMPLEMENTED(opts_ptr->discretization);

    (void) L_ptr; // to avoid unused parameter warning

    COMPLEX scatter_matrix[4] = {0};
    COMPLEX * r = NULL;
    COMPLEX lambda = 0;
    REAL Ei = 0, Ei_prev = 0;
    REAL DEL = 0, DEL_prev = 0;
    INT W = 0;
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
    
    switch (opts_ptr->mainspec_type) {
    case kdvp_mstype_AMPLITUDES_AND_MODULI:
        // fallthrough
    case kdvp_mstype_OPENBANDS:
        // fallthrough
    case kdvp_mstype_EDGEPOINTS_AND_SIGNS:
        Ei = E[0];
        if (Ei < 0)
            lambda = I*SQRT(-Ei);
        else
            lambda = SQRT(Ei);
        ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
            &lambda, scatter_matrix, &W, opts_ptr->discretization, 0/*derivative_flag*/);
        CHECK_RETCODE(ret_code, leave_fun);
        Ei_prev = E[0];
        DEL_prev = 0.5*POW(2, W)*CREAL(scatter_matrix[0] + scatter_matrix[3]);

        UINT K = 0;
        for (i=1; i<L; i++) {
            Ei = E[0] + i*eps_E;
            if (Ei < 0)
                lambda = I*SQRT(-Ei);
            else
                lambda = SQRT(Ei);
            ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*kappa*/, 1/*K*/,
                    &lambda, scatter_matrix, &W, opts_ptr->discretization, 0/*derivative_flag*/);
            CHECK_RETCODE(ret_code, leave_fun);

            DEL = 0.5*POW(2, W)*CREAL(scatter_matrix[0] + scatter_matrix[3]);
            REAL s = 0;
            if (DEL_prev<=-1 && DEL>=1)
                s = 2; // two edge points with delta increasing
            else if (DEL_prev>=1 && DEL<=-1)
                s = -2; // two edge points with delta decreasing
            else if ((DEL_prev<=1 && DEL>=1) || (DEL_prev>=1 && DEL<=1))
                s = 1; // one edge point s.t. delta=1
            else if ((DEL_prev>=-1 && DEL<=-1) || (DEL_prev<=-1 && DEL>=-1))
                s = -1; // one edge point s.t. delta=-1

            if (s == 2 || s == -2) { // two edge points between two grid points
                if (K+1>=*K_ptr) {
                    ret_code = E_OTHER("More than *K_ptr initial guesses for bound states found. Increase *K_ptr and try again.")
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
            Ei_prev = Ei;
            DEL_prev = DEL;
        }

        if (opts_ptr->mainspec_type == kdvp_mstype_EDGEPOINTS_AND_SIGNS) {
            *K_ptr = K;
            break;
        }

        if (opts_ptr->mainspec_type == kdvp_mstype_OPENBANDS) {
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
        }

        // mainspec_type == AMPLITUDES_AND_MODULI
 
        if (K < 3) { // nothing to do
            *K_ptr = 0;
            break;
        }

        // replace main spectrum points E_i with amplitudes and signs with moduli
        UINT cnt = 0;
        UINT i_ref = 0;
        REAL E_ref = FNFT_INF; // reference level
        for (i=1; i<K-1; i+=2) {
            const REAL E_2i = main_spec[2*i];
            const REAL E_2ip1 = main_spec[2*(i+1)];
            const REAL E_2im1 = main_spec[2*(i-1)];
            const REAL modulus = (E_2ip1 - E_2i)/(E_2ip1 - E_2im1);
            main_spec[2*cnt+1] = modulus;
            if (modulus >= 0.99) { // soliton
                i_ref = cnt;
                E_ref = E_2ip1;
                main_spec[2*cnt] = E_2i; // reference level is added later
            } else { // radiation
                main_spec[2*cnt] = 0.5*(E_2ip1 - E_2i);
            }
            cnt++;
        }

        if (E_ref < FNFT_INF) { // soliton reference level still needs to be accounted for
            for (i=0; i<=i_ref; i++)
                main_spec[2*i] = 2*(E_ref - main_spec[2*i]);
        }

        *K_ptr = cnt;
        break;

    case kdvp_mstype_FLOQUET:
        for (i=0; i<L; i++) {
            Ei = E[0] + i*eps_E;
            if (Ei < 0)
                lambda = I*SQRT(-Ei);
            else
                lambda = SQRT(Ei);
            ret_code = kdv_scatter_matrix(D, q, r, eps_t, 0/*E*/, 1/*K*/,
                    &lambda, scatter_matrix, &W, opts_ptr->discretization, 0/*derivative_flag*/);
            CHECK_RETCODE(ret_code, leave_fun);

            DEL = 0.5*POW(2, W)*CREAL(scatter_matrix[0] + scatter_matrix[3]);
            main_spec[i] = DEL;
        }
        break;
    }

leave_fun:
    free(r);
    return ret_code;
}
