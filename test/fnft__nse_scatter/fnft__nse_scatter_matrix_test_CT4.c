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
 * Igor Chekhovskoy (NSU, FRC ICT) 2026.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static void matrix2_mult(COMPLEX const * const A,
                         COMPLEX const * const B,
                         COMPLEX * const C)
{
    C[0] = A[0]*B[0] + A[1]*B[2];
    C[1] = A[0]*B[1] + A[1]*B[3];
    C[2] = A[2]*B[0] + A[3]*B[2];
    C[3] = A[2]*B[1] + A[3]*B[3];
}

static void matrix2_exp_q(COMPLEX const q, COMPLEX const r,
                          COMPLEX const lambda, REAL const h,
                          COMPLEX * const E)
{
    const COMPLEX k = CSQRT(q*r-lambda*lambda);
    const COMPLEX c = CCOSH(k*h);
    const COMPLEX s = h*misc_CSINC(I*k*h);
    E[0] = c-I*lambda*s;
    E[1] = q*s;
    E[2] = r*s;
    E[3] = c+I*lambda*s;
}

static INT matrix2_inverse(COMPLEX const * const A, COMPLEX * const Ainv)
{
    const COMPLEX det = A[0]*A[3]-A[1]*A[2];
    if (det == 0.0)
        return E_DIV_BY_ZERO;
    Ainv[0] = A[3]/det;
    Ainv[1] = -A[1]/det;
    Ainv[2] = -A[2]/det;
    Ainv[3] = A[0]/det;
    return SUCCESS;
}

/* Direct transcription of Eq. 17, independent of the library CT4 helper. */
static INT ct4_oracle(COMPLEX const * const q, COMPLEX const * const r,
                      COMPLEX const lambda, REAL const eps_t,
                      COMPLEX * const result)
{
    COMPLEX E[4], Em[4], H[4], Dp[4] = {0,q[1]-q[0],r[1]-r[0],0};
    COMPLEX Dm[4] = {0,q[2]-q[0],r[2]-r[0],0};
    COMPLEX tmp[4], Mp[4], Mm[4], X[4], Y[4], Xinv[4], C[4];
    INT ret_code;
    matrix2_exp_q(q[0],r[0],lambda,eps_t,E);
    matrix2_exp_q(q[0],r[0],lambda,-eps_t,Em);
    matrix2_exp_q(q[0],r[0],lambda,0.5*eps_t,H);
    matrix2_mult(Em,Dp,tmp);
    matrix2_mult(tmp,E,Mp);
    matrix2_mult(E,Dm,tmp);
    matrix2_mult(tmp,Em,Mm);
    for (UINT i=0; i<4; i++) {
        const COMPLEX A = eps_t*(Mp[i]+Mm[i])/48.0;
        X[i] = -A;
        Y[i] = A;
    }
    X[0] += 1.0; X[3] += 1.0;
    Y[0] += 1.0; Y[3] += 1.0;
    ret_code = matrix2_inverse(X,Xinv);
    CHECK_RETCODE(ret_code, leave_fun);
    matrix2_mult(Xinv,Y,C);
    matrix2_mult(H,C,tmp);
    matrix2_mult(tmp,H,result);
leave_fun:
    return ret_code;
}

static INT test_one_step(const INT kappa)
{
    const REAL eps_t = 0.17;
    const COMPLEX lambda = 0.37+0.19*I;
    COMPLEX q[3] = {0.51-0.13*I, 0.48+0.09*I, 0.55-0.21*I};
    COMPLEX r[3], exact[4], result[8], plus[4], minus[4];
    INT W = 0, ret_code;
    const REAL delta = 2e-6;
    for (UINT i=0; i<3; i++)
        r[i] = -kappa*CONJ(q[i]);
    ret_code = ct4_oracle(q,r,lambda,eps_t,exact);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nse_scatter_matrix(3,q,r,eps_t,kappa,1,&lambda,result,NULL,
            nse_discretization_CT4,1);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(4,result,exact) > 100*EPSILON)
        return E_TEST_FAILED;

    ret_code = ct4_oracle(q,r,lambda+delta,eps_t,plus);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = ct4_oracle(q,r,lambda-delta,eps_t,minus);
    CHECK_RETCODE(ret_code, leave_fun);
    for (UINT i=0; i<4; i++)
        exact[i] = (plus[i]-minus[i])/(2*delta);
    if (misc_rel_err(4,&result[4],exact) > 2e-9)
        return E_TEST_FAILED;

    ret_code = nse_scatter_matrix(3,q,r,eps_t,kappa,1,&lambda,plus,&W,
            nse_discretization_CT4,0);
    CHECK_RETCODE(ret_code, leave_fun);
    for (UINT i=0; i<4; i++)
        plus[i] *= POW(2,W);
    if (misc_rel_err(4,plus,result) > 100*EPSILON)
        return E_TEST_FAILED;
leave_fun:
    return ret_code;
}

static INT test_constant_and_invariant(const INT kappa)
{
    const REAL eps_t = 0.11;
    const COMPLEX lambda = 0.43;
    COMPLEX q[3] = {0.35+0.12*I,0.35+0.12*I,0.35+0.12*I};
    COMPLEX r[3], ct4[4], bo[4];
    INT ret_code;
    for (UINT i=0; i<3; i++)
        r[i] = -kappa*CONJ(q[i]);
    ret_code = nse_scatter_matrix(3,q,r,eps_t,kappa,1,&lambda,ct4,NULL,
            nse_discretization_CT4,0);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nse_scatter_matrix(1,q,r,eps_t,kappa,1,&lambda,bo,NULL,
            nse_discretization_BO,0);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(4,ct4,bo) > 100*EPSILON)
        return E_TEST_FAILED;
    q[1] = 0.31-0.09*I;
    q[2] = 0.42+0.17*I;
    for (UINT i=0; i<3; i++)
        r[i] = -kappa*CONJ(q[i]);
    ret_code = nse_scatter_matrix(3,q,r,eps_t,kappa,1,&lambda,ct4,NULL,
            nse_discretization_CT4,0);
    CHECK_RETCODE(ret_code, leave_fun);
    if (FABS(CABS(ct4[0])*CABS(ct4[0])
            + kappa*CABS(ct4[1])*CABS(ct4[1])-1.0) > 100*EPSILON)
        return E_TEST_FAILED;
leave_fun:
    return ret_code;
}

INT main()
{
    const COMPLEX q_bad[2] = {0.1,0.2};
    const COMPLEX r_bad[2] = {-0.1,-0.2};
    const COMPLEX lambda = 0.3;
    COMPLEX result[4];
    if (nse_scatter_matrix(2,q_bad,r_bad,0.1,+1,1,&lambda,result,NULL,
            nse_discretization_CT4,0) == SUCCESS)
        return EXIT_FAILURE;
    INT ret_code = test_one_step(+1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_one_step(-1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_constant_and_invariant(+1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_constant_and_invariant(-1);
    CHECK_RETCODE(ret_code, failure);
    return EXIT_SUCCESS;
failure:
    return EXIT_FAILURE;
}
