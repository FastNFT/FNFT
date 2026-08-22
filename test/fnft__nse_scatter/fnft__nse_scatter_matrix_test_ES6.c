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
 * Igor Chekhovskoy 2026.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__nse_scatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static void es6_pack(COMPLEX const samples[5], REAL const eps_t,
                     COMPLEX packed[5])
{
    packed[0] = eps_t*samples[2];
    packed[1] = eps_t*(-samples[4]+8.0*samples[3]
            -8.0*samples[1]+samples[0])/12.0;
    packed[2] = eps_t*(-samples[4]+16.0*samples[3]
            -30.0*samples[2]+16.0*samples[1]-samples[0])/12.0;
    packed[3] = eps_t*(samples[4]-2.0*samples[3]
            +2.0*samples[1]-samples[0])/2.0;
    packed[4] = eps_t*(samples[4]-4.0*samples[3]
            +6.0*samples[2]-4.0*samples[1]+samples[0]);
}

/* Direct transcription of Eqs. 65--67, independent of the library ES6
 * implementation and preprocessing helper. */
static void es6_oracle(COMPLEX const q[5], COMPLEX const r[5],
                       COMPLEX const lambda, REAL const eps_t,
                       COMPLEX result[4])
{
    COMPLEX Z[4], delta, c, s;
    const COMPLEX z = eps_t*lambda;
    const COMPLEX z2 = z*z;
    const COMPLEX z3 = z2*z;

    Z[0] = (r[0]*q[1]-q[0]*r[1])*z2/180.0
            -I*(1.0-(r[0]*q[2]+q[0]*r[2])/360.0
            +q[1]*r[1]/60.0)*z
            +(15.0-q[0]*r[0])*(r[0]*q[1]-q[0]*r[1])/180.0
            +(r[0]*q[3]-q[0]*r[3]+q[1]*r[2]-r[1]*q[2])/480.0;
    Z[1] = I*q[1]*z3/90.0-q[2]*z2/180.0
            +I*(q[1]/6.0+q[3]/240.0-q[0]*r[0]*q[1]/90.0)*z
            +q[0]+q[0]*(r[0]*q[2]-q[0]*r[2])/360.0
            +(q[0]*r[1]-r[0]*q[1])*q[1]/120.0
            +q[2]/24.0+q[4]/1920.0;
    Z[2] = -I*r[1]*z3/90.0-r[2]*z2/180.0
            -I*(r[1]/6.0+r[3]/240.0-q[0]*r[0]*r[1]/90.0)*z
            +r[0]+r[0]*(q[0]*r[2]-r[0]*q[2])/360.0
            +(r[0]*q[1]-q[0]*r[1])*r[1]/120.0
            +r[2]/24.0+r[4]/1920.0;
    Z[3] = -Z[0];
    delta = Z[0]*Z[0]+Z[1]*Z[2];
    if (CABS(delta) < 1e-12) {
        const COMPLEX delta2 = delta*delta;
        const COMPLEX delta3 = delta2*delta;
        c = 1.0+delta/2.0+delta2/24.0+delta3/720.0;
        s = 1.0+delta/6.0+delta2/120.0+delta3/5040.0;
    } else {
        const COMPLEX root = CSQRT(delta);
        c = CCOSH(root);
        s = CSINH(root)/root;
    }
    result[0] = c+s*Z[0];
    result[1] = s*Z[1];
    result[2] = s*Z[2];
    result[3] = c+s*Z[3];
}

static INT test_one_step(const INT kappa)
{
    const REAL eps_t = 0.17;
    const REAL delta = 2e-6;
    const COMPLEX lambda = 0.37+0.19*I;
    const COMPLEX samples[5] = {
        0.42-0.07*I, 0.48+0.09*I, 0.51-0.13*I,
        0.55-0.21*I, 0.46+0.16*I
    };
    COMPLEX rsamples[5], q[5], r[5], exact[4], result[8];
    COMPLEX plus[4], minus[4], normalized[4];
    INT W = 0, ret_code;
    UINT i;

    for (i=0; i<5; i++)
        rsamples[i] = -kappa*CONJ(samples[i]);
    es6_pack(samples,eps_t,q);
    es6_pack(rsamples,eps_t,r);
    es6_oracle(q,r,lambda,eps_t,exact);
    ret_code = nse_scatter_matrix(5,q,r,eps_t,kappa,1,&lambda,result,NULL,
            nse_discretization_ES6,1);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(4,result,exact) > 100*EPSILON)
        return E_TEST_FAILED;

    es6_oracle(q,r,lambda+delta,eps_t,plus);
    es6_oracle(q,r,lambda-delta,eps_t,minus);
    for (i=0; i<4; i++)
        exact[i] = (plus[i]-minus[i])/(2.0*delta);
    if (misc_rel_err(4,&result[4],exact) > 2e-9)
        return E_TEST_FAILED;

    ret_code = nse_scatter_matrix(5,q,r,eps_t,kappa,1,&lambda,normalized,&W,
            nse_discretization_ES6,0);
    CHECK_RETCODE(ret_code, leave_fun);
    for (i=0; i<4; i++)
        normalized[i] *= POW(2,W);
    if (misc_rel_err(4,normalized,result) > 100*EPSILON)
        return E_TEST_FAILED;
leave_fun:
    return ret_code;
}

static INT test_constant_zero_and_invariant(const INT kappa)
{
    const REAL eps_t = 0.11;
    const COMPLEX lambda = 0.43;
    COMPLEX samples[5], rsamples[5], q[5], r[5], es6[4], bo[4];
    INT ret_code;
    UINT i;

    for (i=0; i<5; i++) {
        samples[i] = 0.35+0.12*I;
        rsamples[i] = -kappa*CONJ(samples[i]);
    }
    es6_pack(samples,eps_t,q);
    es6_pack(rsamples,eps_t,r);
    ret_code = nse_scatter_matrix(5,q,r,eps_t,kappa,1,&lambda,es6,NULL,
            nse_discretization_ES6,0);
    CHECK_RETCODE(ret_code, leave_fun);
    ret_code = nse_scatter_matrix(1,samples,rsamples,eps_t,kappa,1,&lambda,
            bo,NULL,nse_discretization_BO,0);
    CHECK_RETCODE(ret_code, leave_fun);
    if (misc_rel_err(4,es6,bo) > 100*EPSILON)
        return E_TEST_FAILED;

    samples[0] = 0.42-0.07*I;
    samples[1] = 0.48+0.09*I;
    samples[2] = 0.51-0.13*I;
    samples[3] = 0.55-0.21*I;
    samples[4] = 0.46+0.16*I;
    for (i=0; i<5; i++)
        rsamples[i] = -kappa*CONJ(samples[i]);
    es6_pack(samples,eps_t,q);
    es6_pack(rsamples,eps_t,r);
    ret_code = nse_scatter_matrix(5,q,r,eps_t,kappa,1,&lambda,es6,NULL,
            nse_discretization_ES6,0);
    CHECK_RETCODE(ret_code, leave_fun);
    if (FABS(CABS(es6[0])*CABS(es6[0])
            +kappa*CABS(es6[1])*CABS(es6[1])-1.0) > 200*EPSILON)
        return E_TEST_FAILED;

    for (i=0; i<5; i++) {
        q[i] = 0.0;
        r[i] = 0.0;
    }
    ret_code = nse_scatter_matrix(5,q,r,eps_t,kappa,1,&lambda,es6,NULL,
            nse_discretization_ES6,0);
    CHECK_RETCODE(ret_code, leave_fun);
    bo[0] = CEXP(-I*lambda*eps_t);
    bo[1] = bo[2] = 0.0;
    bo[3] = CEXP(I*lambda*eps_t);
    if (misc_rel_err(4,es6,bo) > 100*EPSILON)
        return E_TEST_FAILED;
leave_fun:
    return ret_code;
}

INT main()
{
    const COMPLEX q_bad[4] = {0};
    const COMPLEX r_bad[4] = {0};
    const COMPLEX lambda = 0.3;
    COMPLEX result[4];
    INT ret_code;

    if (nse_scatter_matrix(4,q_bad,r_bad,0.1,+1,1,&lambda,result,NULL,
            nse_discretization_ES6,0) == SUCCESS)
        return EXIT_FAILURE;
    ret_code = test_one_step(+1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_one_step(-1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_constant_zero_and_invariant(+1);
    CHECK_RETCODE(ret_code, failure);
    ret_code = test_constant_zero_and_invariant(-1);
    CHECK_RETCODE(ret_code, failure);
    return EXIT_SUCCESS;
failure:
    return EXIT_FAILURE;
}
