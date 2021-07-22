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
#include "fnft__misc.h"
#include "fnft__manakov_discretization.h"


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
    if(CABS(q1)<EPSILON && CABS(q2)<EPSILON){
		U[0] = CEXP(-I*lam*eps_t);
		U[1] = 0;
		U[2] = 0;
		U[3] = 0;
		U[4] = CEXP(I*lam*eps_t);
		U[5] = 0;
		U[6] = 0;
		U[7] = 0;
		U[8] = CEXP(I*lam*eps_t);
	}else{
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

}


/**
 * Returns [S11 S12 S13 S21 S22 S23 S31 S32 S33] in result where
 * S = [S11, S12, S13; S21, S22, S23, S31, S32, S33] is the scattering
 * matrix computed using the chosen scheme.
 * result should be preallocated with size 9*K.
 */
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

        COMPLEX *q1_preprocessed = NULL;
        COMPLEX *q2_preprocessed = NULL;
UINT D_eff = D;
    switch (discretization) {
        case manakov_discretization_CF4_2:;         // commutator-free fourth-order. Fall through
        case manakov_discretization_BO:            // bofetta-osborne scheme
            q1_preprocessed = q1;
            q2_preprocessed = q2;
        break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

        // Calculate the scattering matrix
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][3][3] = { {{1,0,0}, {0,1,0}, {0,0,1}} }; // Initialize transfer matrix
            COMPLEX U[3][3];        // transfer matrix for current timestep
            UINT current = 0;

            switch (discretization) {
                case manakov_discretization_CF4_2:
                l_curr /=2;                     // Intional fall through
                case manakov_discretization_BO:
                    for (UINT n = 0; n < D_eff; n++){
                        manakov_scatter_U_BO(q1_preprocessed[n],q2_preprocessed[n],l_curr,eps_t,kappa,*U);
                        misc_matrix_mult(3,3,3,&U[0][0],&H[current][0][0],&H[!current][0][0]);
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
    return ret_code;
}

