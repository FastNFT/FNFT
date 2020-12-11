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
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J Prins (TU Delft) 2020.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__akns_scatter.h"

/**
 * Auxiliary routines, used by the main routines below
 */
static inline void akns_scatter_bound_states_U_BO(COMPLEX const qn,
                                                  COMPLEX const rn,
                                                  COMPLEX const ln,
                                                  COMPLEX const eps_t,
                                                  UINT const derivative_flag,
                                                  COMPLEX * const U)
{
    COMPLEX const ln2 = ln * ln;
    COMPLEX const ks = ((qn*rn)-ln2);
    COMPLEX const k = CSQRT(ks);
    COMPLEX const ch = CCOSH(k*eps_t);
    COMPLEX const sh = eps_t * misc_CSINC(I*k*eps_t);
    COMPLEX const u1 = ln*sh*I;

    U[0] = ch - u1;
    U[1] = qn*sh;
    if (derivative_flag) {
        U[4] = rn*sh;
        U[5] = ch + u1;
        memcpy(&U[10],&U[0],6 * sizeof(COMPLEX)); // lower right block

        COMPLEX const chi = ch/ks;
        COMPLEX const ud1 = eps_t*ln2*chi*I;
        COMPLEX const ud2 = misc_CSINC_derivative(I*eps_t*k)*I*ln*eps_t*eps_t/k;

        U[8]  = ud1 - ( ln*eps_t + I + (ln2*I)/ks )*sh;
        U[9]  = -qn*ud2;
        U[12] = -rn*ud2;
        U[13] = -ud1 - ( ln*eps_t - I - (ln2*I)/ks )*sh;
    } else {
        U[2] = rn*sh;
        U[3] = ch + u1;
    }
}

static inline void akns_scatter_bound_states_U_ES4(COMPLEX const a1,
                                                   COMPLEX const a2,
                                                   COMPLEX const a3,
                                                   UINT const derivative_flag,
                                                   COMPLEX * const U,
                                                   COMPLEX * const w,
                                                   COMPLEX * const s,
                                                   COMPLEX * const c)
{
    *w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
    *s = misc_CSINC(*w);
    *c = CCOS(*w);
    U[0]                   = *c + *s*a3;
    U[1]                   = *s*(a1 - I*a2);
    U[derivative_flag?4:2] = *s*(a1 + I*a2);
    U[derivative_flag?5:3] = *c - *s*a3;
}

/**
 * If derivative_flag=0 returns [S11 S12 S21 S22] in result where
 * S = [S11, S12; S21, S22] is the scattering matrix computed using the
 * chosen scheme.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in
 * result where S11' is the derivative of S11 w.r.t to lambda.
 * Result should be preallocated with size 4*K or 8*K accordingly.
 */
INT akns_scatter_matrix(UINT const D,
                        COMPLEX const * const q,
                        COMPLEX const * const r,
                        REAL const eps_t,
                        UINT const K,
                        COMPLEX const * const lambda,
                        COMPLEX * const result,
                        akns_discretization_t const discretization,
                        UINT const derivative_flag)
{
    INT ret_code = SUCCESS;

    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (!(eps_t > 0))
        return E_INVALID_ARGUMENT(eps_t);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (lambda == NULL)
        return E_INVALID_ARGUMENT(lambda);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (derivative_flag != 0 && derivative_flag != 1)
        return E_INVALID_ARGUMENT(derivative_flag);
    UINT const upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *l = NULL;

    // Define stepsize constants that are often needed
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.

    COMPLEX l_weights[4] = {0};
    UINT N = 0;
    switch (discretization) {
        case akns_discretization_ES4:
            tmp1 = derivative_flag ? malloc(2*D*sizeof(COMPLEX)) : malloc(D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            for (UINT n = 0; n < D; n+=3){
                tmp1[n] = eps_t_3*(q[n+2]+r[n+2])/48.0 + (eps_t*(q[n]+r[n]))*0.5;
                tmp1[n+1] = (eps_t*(q[n]-r[n])*I)*0.5 + (eps_t_3*(q[n+2]-r[n+2])*I)/48.0;
                tmp1[n+2] = -eps_t_3*(q[n]*r[n+1]- q[n+1]*r[n])/12.0;
            }
            if (derivative_flag == 1){
                tmp2 = &tmp1[D];
                for (UINT n = 0; n < D; n+=3){
                    tmp2[n] = I*eps_t_3*(q[n+1]-r[n+1])/12.0;
                    tmp2[n+1] = -eps_t_3*(q[n+1]+r[n+1])/12.0;
                    tmp2[n+2] = -I*eps_t;
                }
            }
            break;

        case akns_discretization_TES4:
            tmp1 = malloc(2*D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            tmp2 = &tmp1[D];
            for (UINT n = 0; n < D; n+=3){
                tmp1[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp1[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                tmp2[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp2[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
            }
            break;

        case akns_discretization_CF4_3:         // commutator-free fourth-order
        case akns_discretization_CF5_3:         // commutator-free fifth-order
        case akns_discretization_CF6_4:         // commutator-free sixth-order
            N++;                                // The previous three discretizations require N=3
            // fall through
        case akns_discretization_CF4_2:         // commutator-free fourth-order
            N++;                                // The previous discretization requires N=2
            // fall through
        case akns_discretization_BO:            // bofetta-osborne scheme
            N++;                                // The previous discretization requires N=1
            COMPLEX *weights = NULL;
            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, leave_fun); // if ret_code != SUCCESS, akns_discretization_method_weights frees weights if needed

            for (UINT i = 0; i < upsampling_factor; i++){
                for (UINT j = 0; j < N; j++)
                    l_weights[i] += weights[i*N+j];
            }
            free(weights);
#ifdef DEBUG
            misc_print_buf(upsampling_factor,l_weights,"l_weights");
#endif

            // Allocating memory for the array l, which will be used to store xi-samples
            // with the appropriate weight for the discretization.
            l = malloc(D*sizeof(COMPLEX));
            if (l == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            break;
            
        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

    if (derivative_flag){
        // Colculate the scattering matrix with lamda-derivative as in G. Boffetta an A.R. Osborne, 'Computation of the direct scattering transform for the nonlinear Schroedinger equation', www.doi.org/10.1016/0021-9991(92)90370-e .
        COMPLEX U[4][4] = {{ 0 }};
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][4][4] = { { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1} } }; // Init only first sixteen values
            UINT current = 0;

            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    for (UINT n = 0; n < D; n +=upsampling_factor ) {
                        for (UINT j = 0; j < upsampling_factor; j++ )
                            l[n+j] = l_curr*l_weights[j];
                    }
                    for (UINT n = 0; n < D; n++){
                        akns_scatter_bound_states_U_BO(q[n],r[n],l[n],eps_t,1,*U);
                        U[2][0] *= l_weights[n%upsampling_factor];
                        U[2][1] *= l_weights[n%upsampling_factor];
                        U[3][0] *= l_weights[n%upsampling_factor];
                        U[3][1] *= l_weights[n%upsampling_factor];
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                    }
                    break;

                case akns_discretization_ES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX w, s, c;
                        COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                        akns_scatter_bound_states_U_ES4(a1,a2,a3,1,*U,&w,&s,&c);
                        memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                        COMPLEX w_d = -(a1*tmp2[n]+a2*tmp2[n+1]+a3*tmp2[n+2]);
                        COMPLEX c_d = -misc_CSINC(w)*w_d;
                        w_d /= w;
                        COMPLEX s_d = w_d * misc_CSINC_derivative(w);
                        U[2][0] = c_d+s_d*a3+s*tmp2[n+2];
                        U[2][1] = s_d*a1+s*tmp2[n]-I*s_d*a2-I*s*tmp2[n+1];
                        U[3][0] = s_d*a1+s*tmp2[n]+I*s_d*a2+I*s*tmp2[n+1];
                        U[3][1] = c_d-s_d*a3-s*tmp2[n+2];

                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                    }
                    break;
                case akns_discretization_TES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX M[2][2], w, s, c;

                        akns_scatter_bound_states_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*M,&w,&s,&c);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][0][0],&H[!current][0][0]);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][2][0],&H[!current][2][0]);
                        current = !current;

                        akns_scatter_bound_states_U_ES4(0.5*eps_t*(q[n]+r[n]),
                                                        0.5*eps_t*(q[n]-r[n])*I,
                                                        -eps_t*l_curr*I, 1,*U,&w,&s,&c);
                        memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                        COMPLEX s_d = CSIN(w*eps_t)/w;
                        COMPLEX c_d = -eps_t*l_curr*s_d;
                        COMPLEX w_d = l_curr*(eps_t*w*CCOS(w*eps_t)-CSIN(w*eps_t))/(w*w*w);
                        U[2][0] = c_d-I*s_d;
                        U[2][1] = w_d*q[n];
                        U[3][0] = w_d*r[n];
                        U[3][1] = c_d+I*s_d;
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][2][0],&H[!current][2][0]);
                        current = !current;

                        akns_scatter_bound_states_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*M,&w,&s,&c);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][0][0],&H[!current][0][0]);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][2][0],&H[!current][2][0]);
                        current = !current;
                    }
                    break;

                default: // Unknown discretization
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
            // TODO: Change of basis

            // Copy result
            result[i*8] = H[current][0][0];
            result[i*8 + 1] = H[current][0][1];
            result[i*8 + 2] = H[current][1][0];
            result[i*8 + 3] = H[current][1][1];
            result[i*8 + 4] = H[current][2][0];
            result[i*8 + 5] = H[current][2][1];
            result[i*8 + 6] = H[current][3][0];
            result[i*8 + 7] = H[current][3][1];
        }
    } else {
        // Colculate the scattering matrix without lambda-derivative
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][2][2] = { { {1,0}, {0,1} } }; // Init only first four values
            UINT current = 0;

            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    for (UINT n = 0; n < D; n +=upsampling_factor ) {
                        for (UINT j = 0; j < upsampling_factor; j++ )
                            l[n+j] = l_curr*l_weights[j];
                    }
                    for (UINT n = 0; n < D; n++){
                        COMPLEX U[2][2];
                        akns_scatter_bound_states_U_BO(q[n],r[n],l[n],eps_t,0,*U);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                    }
                    break;

                case akns_discretization_ES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX U[2][2], w, s, c;

                        COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];

                        akns_scatter_bound_states_U_ES4(a1,a2,a3,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                    }
                    break;

                case akns_discretization_TES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX U[2][2], w, s, c;

                        akns_scatter_bound_states_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        akns_scatter_bound_states_U_ES4(0.5*eps_t*(q[n]+r[n]),
                                                        0.5*eps_t*(q[n]-r[n])*I,
                                                        -eps_t*l_curr*I, 0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        akns_scatter_bound_states_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                    }
                    break;

                default: // Unknown discretization
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
            // TODO: Change of basis

            // Copy result
            result[i*4] = H[current][0][0];
            result[i*4 + 1] = H[current][0][1];
            result[i*4 + 2] = H[current][1][0];
            result[i*4 + 3] = H[current][1][1];
        }
    }
    
leave_fun:
    free(l);
    free(tmp1);
    return ret_code;
}

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 */
INT akns_scatter_bound_states(UINT const D,
                              COMPLEX const * const q,
                              COMPLEX const * const r,
                              REAL const *const T,
                              UINT const K,
                              COMPLEX const * const bound_states,
                              COMPLEX * const a_vals,
                              COMPLEX * const aprime_vals,
                              COMPLEX * const b,
                              akns_discretization_t const discretization,
                              akns_pde_t const PDE,
                              UINT const vanilla_flag,
                              UINT const skip_b_flag)
{
    INT ret_code = SUCCESS;

    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (r == NULL)
        return E_INVALID_ARGUMENT(r);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (K <= 0.0)
        return E_INVALID_ARGUMENT(K);
    if (bound_states == NULL)
        return E_INVALID_ARGUMENT(bound_states);
    if (a_vals == NULL)
        return E_INVALID_ARGUMENT(a);
    if (aprime_vals == NULL)
        return E_INVALID_ARGUMENT(a_prime);
    if (b == NULL)
        return E_INVALID_ARGUMENT(b);
    UINT const upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    REAL const boundary_coeff = akns_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(>discretization);
    if (D%upsampling_factor != 0)
        return E_ASSERTION_FAILED;
    UINT const D_given = D/upsampling_factor;
    if (PDE!=akns_pde_KdV && PDE!=akns_pde_NSE)
        return E_INVALID_ARGUMENT(PDE);

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp4 = NULL, *l = NULL;

    // Allocating memory for storing PHI and PSI at all D_given points as
    // there are required to find the right value of b.
    // First, we will store the values of PHI and its xi-derivative as follows:
    // PSIPHI = [*,*,PHI1[0],PHI2[0],PHI1_D[0],PHI2_D[0],PHI1[1],PHI2[1],PHI1_D[1],PHI2_D[1], ... ,,PHI1[D_given-1],PHI2[D_given-1],PHI1_D[D_given-1],PHI2_D[D_given-1]]
    // Next, we will overwrite the derivatives that we don't need anymore
    // (all except for those at D_given) to store PSI:
    // PSIPHI = [PSI1[0],PSI2[0],PHI1[0],PHI2[0],PSI1[1],PSI2[1],PHI1[1],PHI2[1], ... PSI1[D_given-1],PSI2[D_given-1],PHI1[D_given-1],PHI2[D_given-1],PHI1_D[D_given-1],PHI2_D[D_given-1]]
    // This keeps all vectors in adjacent memory locations, such that we can
    // use matrix-vector multiplication.
    COMPLEX * const PSIPHI = malloc((4*(D_given+1)+2) * sizeof(COMPLEX));
    if (PSIPHI == NULL) {
        ret_code = E_NOMEM;
        CHECK_RETCODE(ret_code, leave_fun);
    }
    COMPLEX * const PSI = &PSIPHI[0];
    COMPLEX * const PHI = &PSIPHI[2];

    // Define stepsize constants that are often needed
    REAL const eps_t = (T[1] - T[0])/(D_given - 1);
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.

    REAL scl_factor = 1.0;
    COMPLEX l_weights[4] = {0};
    UINT N = 0;
    switch (discretization) {
        //  Fourth-order exponential method which requires
        // one matrix exponential. The matrix exponential is
        // implmented by using the expansion of the 2x2 matrix
        // in terms of Pauli matrices.
        case akns_discretization_ES4:
            tmp1 = malloc(2*D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            tmp2 = &tmp1[D];
            for (UINT n=0; n<D; n+=3){
                tmp1[n] = eps_t_3*(q[n+2]+r[n+2])/48.0 + (eps_t*(q[n]+r[n]))*0.5;
                tmp1[n+1] = (eps_t*(q[n]-r[n])*I)*0.5 + (eps_t_3*(q[n+2]-r[n+2])*I)/48.0;
                tmp1[n+2] = -eps_t_3*(q[n]*r[n+1]- q[n+1]*r[n])/12.0;

                tmp2[n] = I*eps_t_3*(q[n+1]-r[n+1])/12.0;
                tmp2[n+1] = -eps_t_3*(q[n+1]+r[n+1])/12.0;
                tmp2[n+2] = -I*eps_t;
            }
            break;

            //  Fourth-order exponential method which requires
            // three matrix exponentials. The matrix exponential is
            // implmented by using the expansion of the 2x2 matrix
            // in terms of Pauli matrices.
        case akns_discretization_TES4:
            tmp1 = skip_b_flag ? malloc(2*D*sizeof(COMPLEX)) : malloc(4*D*sizeof(COMPLEX));
            if (tmp1 == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            tmp2 = &tmp1[D];
            for (UINT n=0; n<D; n+=3){
                tmp1[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp1[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                tmp2[n] = (eps_t_3*(q[n+2]+r[n+2]))/96.0 + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                tmp2[n+1] = (eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
            }
            if (!skip_b_flag){
                tmp3 = &tmp1[2*D];
                tmp4 = &tmp1[3*D];
                if (tmp3 == NULL || tmp4 == NULL) {
                    ret_code = E_NOMEM;
                    CHECK_RETCODE(ret_code, leave_fun);
                }
                for (UINT n = 0; n < D; n+=3){
                    tmp3[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0 - (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp3[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0 + (eps_t_2*(r[n+1]-q[n+1])*I)/24.0;
                    tmp4[n] = (-eps_t_3*(q[n+2]+r[n+2]))/96.0  + (eps_t_2*(q[n+1]+r[n+1]))/24.0;
                    tmp4[n+1] = (-eps_t_3*(q[n+2]-r[n+2])*I)/96.0  + (eps_t_2*(q[n+1]-r[n+1])*I)/24.0;
                }
            }
            break;

        case akns_discretization_CF4_3:         // commutator-free fourth-order
        case akns_discretization_CF5_3:         // commutator-free fifth-order
        case akns_discretization_CF6_4:         // commutator-free sixth-order
            N++;                                // The previous three discretizations require N=3
            // fall through
        case akns_discretization_CF4_2:         // commutator-free fourth-order
            N++;                                // The previous discretization requires N=2
            // fall through
        case akns_discretization_BO:            // bofetta-osborne scheme
            N++;                                // The previous discretization requires N=1
            scl_factor /= upsampling_factor;
            COMPLEX *weights = NULL;
            ret_code = akns_discretization_method_weights(&weights,discretization);
            CHECK_RETCODE(ret_code, leave_fun); // if ret_code != SUCCESS, akns_discretization_method_weights frees weights if needed

            for (UINT i = 0; i < upsampling_factor; i++){
                for (UINT j = 0; j < N; j++)
                    l_weights[i] += weights[i*N+j];
            }
            free(weights);

            // Allocating memory for the array l, which will be used to store xi-samples
            // with the appropriate weight for the discretization.
            l = malloc(D*sizeof(COMPLEX));
            if (l == NULL) {
                ret_code = E_NOMEM;
                CHECK_RETCODE(ret_code, leave_fun);
            }
            break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

    for (UINT neig=0; neig<K; neig++) { // iterate over bound states
        COMPLEX l_curr = bound_states[neig];

        // Scattering PHI and PHI_D from T[0]-eps_t/2 to T[1]+eps_t/2
        // PHI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        // Set initial condition for PHI in S basis:
        COMPLEX f_S[4];
        f_S[0] = 1.0*CEXP(-I*l_curr*(T[0]-eps_t*boundary_coeff));
        f_S[1] = 0.0;
        f_S[2] = f_S[0]*(-I*(T[0]-eps_t*boundary_coeff));
        f_S[3] = 0.0;

        // Fetch the change of basis matrix from S to the basis of the discretization
        COMPLEX Tmx[4][4];
        UINT derivative_flag = 1;
        ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
        CHECK_RETCODE(ret_code, leave_fun);

        // Calculate the initial condition for PHI in the basis of the discretization
        misc_matrix_mult(4,4,1,&Tmx[0][0],&f_S[0],&PHI[4*0 + 0]);

        // Declaring chonge of state matrix here, to avoid letting them be
        // overwritten with zeros in every loop iteration.
        COMPLEX U[4][4] = {{0}};
        switch (discretization) {

            case akns_discretization_BO:
            case akns_discretization_CF4_2:
            case akns_discretization_CF4_3:
            case akns_discretization_CF5_3:
            case akns_discretization_CF6_4:
                for (UINT n=0; n<D; n+=upsampling_factor){
                    for (UINT i=0; i<upsampling_factor; i++)
                        l[n+i] = l_curr*l_weights[i];
                }
                COMPLEX phi_temp[2][4];
                UINT current = 0;
                memcpy(&phi_temp[current][0], PHI, 4 * sizeof(COMPLEX));
                for (UINT n_given=0; n_given<D_given; n_given++) {
                    for (UINT count=0; count<upsampling_factor; count++) {
                        UINT n = n_given * upsampling_factor + count;
                        akns_scatter_bound_states_U_BO(q[n],r[n],l[n],eps_t,1,*U);
                        U[2][0] *= l_weights[count];
                        U[2][1] *= l_weights[count];
                        U[3][0] *= l_weights[count];
                        U[3][1] *= l_weights[count];
                        misc_matrix_mult(4,4,1,&U[0][0],&phi_temp[current][0],&phi_temp[!current][0]);
                        current = !current;
                    }
                    memcpy(&PHI[4*(n_given+1)], &phi_temp[current][0], 4 * sizeof(COMPLEX));
                }
                break;
                //  Fourth-order exponential method which requires
                // one matrix exponential. The matrix exponential is
                // implmented by using the expansion of the 2x2 matrix
                // in terms of Pauli matrices.
            case akns_discretization_ES4:
                for (UINT n = 0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX w, s, c;
                    COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                    COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                    COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                    akns_scatter_bound_states_U_ES4(a1,a2,a3,1,*U,&w,&s,&c);
                    memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                    COMPLEX w_d = -(a1*tmp2[n]+a2*tmp2[n+1]+a3*tmp2[n+2]);
                    COMPLEX c_d = -misc_CSINC(w)*w_d;
                    w_d /= w;
                    COMPLEX s_d = w_d * misc_CSINC_derivative(w);
                    U[2][0] = c_d+s_d*a3+s*tmp2[n+2];
                    U[2][1] = s_d*a1+s*tmp2[n]-I*s_d*a2-I*s*tmp2[n+1];
                    U[3][0] = s_d*a1+s*tmp2[n]+I*s_d*a2+I*s*tmp2[n+1];
                    U[3][1] = c_d-s_d*a3-s*tmp2[n+2];
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                }
                break;
                // Fourth-order exponential method which requires
                // three matrix exponentials. The transfer metrix
                // needs to be built differently compared to the CF schemes.
            case akns_discretization_TES4:
                for (UINT n=0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX phi_temp[4], M[2][2], w, s, c;

                    akns_scatter_bound_states_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*M,&w,&s,&c);
                    misc_matrix_mult(2,2,1,*M,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    misc_matrix_mult(2,2,1,*M,&PHI[4*n_given+2],&PHI[4*(n_given+1)+2]);

                    akns_scatter_bound_states_U_ES4(0.5*eps_t*(q[n]+r[n]),
                                                    0.5*eps_t*(q[n]-r[n])*I,
                                                    -eps_t*l_curr*I, 1,*U,&w,&s,&c);
                    memcpy(&U[2][2],&U[0][0],6 * sizeof(COMPLEX)); // lower right block
                    COMPLEX s_d = eps_t * misc_CSINC(w);
                    COMPLEX c_d = -eps_t*l_curr*s_d;
                    COMPLEX w_d = l_curr*eps_t_2/w * misc_CSINC_derivative(eps_t*w);
                    U[2][0] = c_d-I*s_d;
                    U[2][1] = w_d*q[n];
                    U[3][0] = w_d*r[n];
                    U[3][1] = c_d+I*s_d;
                    misc_matrix_mult(4,4,1,*U,&PHI[4*(n_given+1)],phi_temp);

                    akns_scatter_bound_states_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*M,&w,&s,&c);
                    misc_matrix_mult(2,2,1,*M,&phi_temp[0],&PHI[4*(n_given+1)]);
                    misc_matrix_mult(2,2,1,*M,&phi_temp[2],&PHI[4*(n_given+1)+2]);
                }
                break;

            default: // Unknown discretization
                ret_code = E_INVALID_ARGUMENT(discretization);
                CHECK_RETCODE(ret_code, leave_fun);
        }

        // If b-coefficient is requested skip_b_flag will not be set.
        // Scattering PSI from T[1]+eps_t/2 to T[0]-eps_t/2.
        // PSI is stored at intermediate values as they are needed for the
        // accurate computation of b-coefficient.
        if (skip_b_flag == 0){

            // Set final condition for PSI in S basis:
            f_S[0] = 0.0;
            f_S[1] = 1.0*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));

            // Fetch the change of basis matrix from S to the basis of the discretization
            derivative_flag = 0;
            ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Calculate the final condition for PSI in the basis of the discretization
            misc_matrix_mult(2,2,1,&Tmx[0][0],&f_S[0],&PSI[4*D_given+0]);

            // Inverse transfer matrix at each step is built by taking
            // negative step -eps_t.
            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    {} // Stop the compiler from complaining about starting a case with a declaration rather than a statement.
                    COMPLEX psi_temp[2][2];
                    UINT current = 0;
                    memcpy(&psi_temp[current][0], &PSI[4*D_given], 2 * sizeof(COMPLEX));
                    for (UINT n_given=D_given; n_given-->0; ) {
                        for (UINT count=upsampling_factor; count-->0; ) {
                            COMPLEX U[2][2];
                            UINT n = n_given * upsampling_factor + count;
                            akns_scatter_bound_states_U_BO(q[n],r[n],l[n],-eps_t,0,*U);
                            misc_matrix_mult(2,2,1,&U[0][0],&psi_temp[current][0],&psi_temp[!current][0]);
                            current = !current;
                        }
                        memcpy(&PSI[4*n_given], psi_temp, 2 * sizeof(COMPLEX));
                    }
                    break;

                    //  Fourth-order exponential method which requires
                    // one matrix exponential. The matrix exponential is
                    // implmented by using the expansion of the 2x2 matrix
                    // in terms of Pauli matrices.
                case akns_discretization_ES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2], w, s, c;
                        COMPLEX a1 = -tmp1[n]- eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = -tmp1[n+1] + eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 =  eps_t*I*l_curr -tmp1[n+2];
                        akns_scatter_bound_states_U_ES4(a1,a2,a3,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                    }
                    break;
                    // Fourth-order exponential method which requires
                    // three matrix exponentials. The transfer metrix cannot
                    // needs to be built differently compared to the CF schemes.
                case akns_discretization_TES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2], psi_temp[2], w, s, c;

                        akns_scatter_bound_states_U_ES4(tmp3[n],tmp3[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);

                        akns_scatter_bound_states_U_ES4(-0.5*eps_t*(q[n]+r[n]),
                                                        -0.5*eps_t*(q[n]-r[n])*I,
                                                        eps_t*l_curr*I,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*n_given],psi_temp);

                        akns_scatter_bound_states_U_ES4(tmp4[n],tmp4[n+1],0.0,0,*U,&w,&s,&c);
                        misc_matrix_mult(2,2,1,*U,psi_temp,&PSI[4*n_given]);
                    }
                    break;

                default: // Unknown discretization

                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
        }

        // Fetch the change of basis matrix from the basis of the discretization to S
        derivative_flag = 1;
        ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
        CHECK_RETCODE(ret_code, leave_fun);

        // Calculate the final state for PHI in S basis
        misc_matrix_mult(4,4,1,&Tmx[0][0],&PHI[4*D_given + 0],&f_S[0]);

        // Calculate the final state for PHI in E basis
        a_vals[neig] = f_S[0] * CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff));
        if (PDE==akns_pde_KdV)
            a_vals[neig] = CREAL(a_vals[neig]);
        aprime_vals[neig] = f_S[2]*CEXP(I*l_curr*(T[1]+eps_t*boundary_coeff))+(I*(T[1]+eps_t*boundary_coeff))* a_vals[neig];
        if (PDE==akns_pde_KdV)
            aprime_vals[neig] = I * CIMAG(aprime_vals[neig]);

        if (skip_b_flag == 0){
            // Calculation of b assuming a=0
            // Uses the metric from DOI: 10.1109/ACCESS.2019.2932256 for choosing the
            // computation point

            // Fetch the change of basis matrix from the basis of the discretization to S
            derivative_flag = 0;
            ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            REAL error_metric = INFINITY, tmp = INFINITY;
            for (UINT n = 0; n <= D_given; n++){
                COMPLEX * const psi_S = &f_S[0], * const phi_S = &f_S[2], b_temp[2];
                // Calculate phi and psi for this sample in S basis
                misc_matrix_mult(2,2,1,&Tmx[0][0],&PHI[4*n + 0], phi_S);
                misc_matrix_mult(2,2,1,&Tmx[0][0],&PSI[4*n + 0], psi_S);
                if (PDE==akns_pde_KdV) {
                    for (UINT i=0; i<4; i++)
                        f_S[i] = CREAL(f_S[i]);
                }
                for (UINT i=0; i<2; i++)
                    b_temp[i] = phi_S[i]/psi_S[i];
                if (PDE!=akns_pde_KdV || CREAL(b_temp[0]*b_temp[1])>0) {
                    tmp = FABS( 0.5* LOG( (REAL)CABS( b_temp[1]/b_temp[0] ) ) );
                    if (tmp < error_metric){
                        b[neig] = b_temp[0];
                        error_metric = tmp;
                    }
                }
            }
        }

    }
leave_fun:
    free(tmp1);
    free(l);
    free(PSIPHI);
    return ret_code;
}
