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
 * Contibutor:
 * Lianne de Vries (TU Delft) 2021.
 * 
 * Contributors original file:
 * Sander Wahls (TU Delft) 2017-2018.
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J Prins (TU Delft) 2020.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakov_scatter.h"

/**
 * Auxiliary routines, used by the main routines below
 */
static inline void manakov_scatter_U_BO(COMPLEX const q1,       // calculates the matrix exp(P*eps_t) for 1 timestep
                                     COMPLEX const q2,
                                     COMPLEX const lam,
                                     COMPLEX const eps_t,
                                     INT const kappa,
                                     COMPLEX * const U)
{
 //   U[0] = (lam*cexp(-eps_t*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2))*I) - lam*cexp(eps_t*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2))*I) + cexp(-eps_t*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2))*I)*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2)) + cexp(eps_t*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2))*I)*csqrt(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2)))/csqrt(2*(lam*lam + kappa*q1*conj(q1) + kappa*q2*conj(q2)));
    COMPLEX x1 = CSQRT(lam*lam + kappa*q1*CONJ(q1) + kappa*q2*CONJ(q2));
    COMPLEX x2  = eps_t*x1*I;
    COMPLEX x3 = q1*CONJ(q1)+q2*CONJ(q2);

        U[0] = -(lam/x1)*CSINH(x2)+CCOSH(x2);
        U[1] = -(I*q1/x1)*CSINH(x2);
        U[2] = -(I*q2/x1)*CSINH(x2);

        U[3] = (I*kappa*CONJ(q1)/x1)*CSINH(x2);
        U[4] = (q1*CONJ(q1)/x3)*(CCOSH(x2) + (lam/x1)*CSINH(x2)) + q2*CONJ(q2)*CEXP(eps_t*lam*I)/x3;
        U[5] = (q2*CONJ(q1)/x3)*(CCOSH(x2) - CEXP(eps_t*lam*I) + (lam/x1)*CSINH(x2));

        U[6] = (I*kappa*CONJ(q2)/x1)*CSINH(x2);
        U[7] = (q1*CONJ(q2)/x3)*(CCOSH(x2) - CEXP(eps_t*lam*I) + (lam/x1)*CSINH(x2));
        U[8] = (q2*CONJ(q2)/x3)*(CCOSH(x2) + (lam/x1)*CSINH(x2)) + q1*CONJ(q1)*CEXP(eps_t*lam*I)/x3;

}


/**
 * If derivative_flag=0 returns [S11 S12 S21 S22] in result where
 * S = [S11, S12; S21, S22] is the scattering matrix computed using the
 * chosen scheme.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in
 * result where S11' is the derivative of S11 w.r.t to lambda.
 * Result should be preallocated with size 4*K or 8*K accordingly.
 */
// TODO: desription
INT manakov_scatter_matrix(UINT const D,
                        COMPLEX const * const q1,
                        COMPLEX const * const q2,
                        REAL const eps_t,
                        UINT const K,
                        COMPLEX const * const lambda,
                        INT const kappa,
                        COMPLEX * const result,
                        manakov_discretization_t const discretization)
{
    INT ret_code = SUCCESS;

    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q1 == NULL)
        return E_INVALID_ARGUMENT(q);
    if (q2 == NULL)
        return E_INVALID_ARGUMENT(r);
    if (!(eps_t > 0))
        return E_INVALID_ARGUMENT(eps_t);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (lambda == NULL)
        return E_INVALID_ARGUMENT(lambda);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    UINT const upsampling_factor = manakov_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL; //*eps_t_scaled = NULL;

    // Define stepsize constants that are often needed
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.
        COMPLEX *q1_preprocessed = NULL;
        COMPLEX *q2_preprocessed = NULL;
    switch (discretization) {
        case manakov_discretization_CF4_2:;         // commutator-free fourth-order. Semicolon added to resolve "a label can only be part of a statement and a declaration is not a statement" error
            // Do the resampling here
            UINT first_last_index[2] = {0};     // not using this somewhere after this, maybe leave out? TODO
            manakov_discretization_preprocess_signal(D, q1, q2, eps_t, kappa,
                    &D, &q1_preprocessed, &q2_preprocessed, first_last_index, discretization);
        break;
        case manakov_discretization_BO:            // bofetta-osborne scheme
            q1_preprocessed = q1;
            q2_preprocessed = q2;
        break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }


//    for (UINT n=0; n<upsampling_factor; n++ )       // TODO: check if this is ok.
//                eps_t_scaled[n] *= eps_t;
    UINT D_eff = upsampling_factor*D;
        // Calculate the scattering matrix without lambda-derivative
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][3][3] = { {{1,0,0}, {0,1,0}, {0,0,1}} }; // Initialize transfer matrix
            COMPLEX U[3][3];        // TM for current timestep
            UINT current = 0;

            switch (discretization) {
                case manakov_discretization_CF4_2:
                l_curr /=2;     // because we are doing BO with 2* the original number of samples,. Intentional fall through
                case manakov_discretization_BO:
                    for (UINT n = 0; n < D_eff; n++){
                        manakov_scatter_U_BO(q1_preprocessed[n],q2_preprocessed[n],l_curr,eps_t,kappa,*U);
                        printf("in loop %d\n",n);
                        misc_print_buf(9,U,"TM");
//                        misc_print_buf(9,H[!current][0][0],"H before mult");
                        misc_matrix_mult(3,3,3,&U[0][0],&H[current][0][0],&H[!current][0][0]);
//                        misc_print_buf(9,H[!current][0][0],"H after mult");
                        current = !current;
                    }
                    break;

                default: // Unknown discretization
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }

            // Copy result
            memcpy(&result[9*i],&H[current][0][0],9 * sizeof(COMPLEX));
        }
    
    
leave_fun:
//    free(eps_t_scaled);
// TODO: should probably do something (free variables?) here
leave_fun_no_eps_t_scaled:
    free(tmp1);
    return ret_code;
}

