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
 * Sander Wahls (TU Delft) 2017-2018, 2023.
 * Shrinivas Chimmalgi (TU Delft) 2017-2020.
 * Peter J Prins (TU Delft) 2020.
 * Sander Wahls (KIT) 2023.
 * Igor Chekhovskoy 2026.
 * Irina Vaseva 2026.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__akns_scatter.h"

/**
 * Auxiliary routines, used by the main routines below
 */
static inline void akns_scatter_U_BO(COMPLEX const qn,
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

static inline void akns_scatter_U_ES4(COMPLEX const a1,
                                      COMPLEX const a2,
                                      COMPLEX const a3,
                                      UINT const derivative_flag,
                                      COMPLEX * const U,
                                      COMPLEX const * const tmp2)
{
    COMPLEX const w = CSQRT(-(a1*a1)-(a2*a2)-(a3*a3));
    COMPLEX const s = misc_CSINC(w);
    COMPLEX const c = CCOS(w);
    U[0] = c + s*a3;
    U[1] = s*(a1 - I*a2);
    if (derivative_flag) {
        U[4] = s*(a1 + I*a2);
        U[5] = c - s*a3;
        memcpy(&U[10],&U[0],6 * sizeof(COMPLEX)); // lower right block

        COMPLEX w_d = -(a1*tmp2[0]+a2*tmp2[1]+a3*tmp2[2]);
        COMPLEX const c_d = -misc_CSINC(w)*w_d;
        w_d /= w;
        COMPLEX const s_d = w_d * misc_CSINC_derivative(w);
        U[8] = c_d+s_d*a3+s*tmp2[2];
        U[9] = s_d*a1+s*tmp2[0]-I*s_d*a2-I*s*tmp2[1];
        U[12] = s_d*a1+s*tmp2[0]+I*s_d*a2+I*s*tmp2[1];
        U[13] = c_d-s_d*a3-s*tmp2[2];
    } else {
        U[2] = s*(a1 + I*a2);
        U[3] = c - s*a3;
    }
}

static inline void akns_scatter_matrix2_mult(COMPLEX const * const A,
                                              COMPLEX const * const B,
                                              COMPLEX * const C)
{
    C[0] = A[0]*B[0] + A[1]*B[2];
    C[1] = A[0]*B[1] + A[1]*B[3];
    C[2] = A[2]*B[0] + A[3]*B[2];
    C[3] = A[2]*B[1] + A[3]*B[3];
}

static inline INT akns_scatter_matrix2_inverse(COMPLEX const * const A,
                                                COMPLEX * const Ainv)
{
    const COMPLEX det = A[0]*A[3] - A[1]*A[2];
    if (det == 0.0)
        return E_DIV_BY_ZERO;
    Ainv[0] = A[3]/det;
    Ainv[1] = -A[1]/det;
    Ainv[2] = -A[2]/det;
    Ainv[3] = A[0]/det;
    return SUCCESS;
}

/**
 * Fourth-order conservative transition matrix, Eq. 17 in
 * https://doi.org/10.1364/OL.44.002264 (also Eq. 50 in
 * https://doi.org/10.1364/OE.377140). The derivative is formed analytically
 * by applying the product and inverse rules to every matrix factor, as in
 * Eq. 58 of the latter reference.
 */
static inline INT akns_scatter_U_CT4(COMPLEX const q,
                                     COMPLEX const r,
                                     COMPLEX const q_plus,
                                     COMPLEX const r_plus,
                                     COMPLEX const q_minus,
                                     COMPLEX const r_minus,
                                     COMPLEX const lambda,
                                     REAL const eps_t,
                                     UINT const derivative_flag,
                                     UINT const inverse_flag,
                                     COMPLEX * const U)
{
    COMPLEX E4[16] = {0}, H4[16] = {0};
    COMPLEX E[4], Ed[4], Em[4] = {0}, Emd[4], H[4], Hd[4];
    COMPLEX Dp[4] = {0.0, q_plus-q, r_plus-r, 0.0};
    COMPLEX Dm[4] = {0.0, q_minus-q, r_minus-r, 0.0};
    COMPLEX tmp[4], tmp2[4], Mp[4], Mm[4], Mpd[4], Mmd[4];
    COMPLEX A[4], Ad[4], X[4], Xinv[4] = {0}, Y[4], C[4], Cd[4];
    COMPLEX T[4], Td[4], Tinv[4], Tdinv[4];
    INT ret_code = SUCCESS;

    if (derivative_flag) {
        akns_scatter_U_BO(q,r,lambda,eps_t,1,E4);
        akns_scatter_U_BO(q,r,lambda,0.5*eps_t,1,H4);
        for (UINT i=0; i<2; i++) {
            for (UINT j=0; j<2; j++) {
                const UINT k = 2*i+j;
                E[k] = E4[4*i+j];
                Ed[k] = E4[8+4*i+j];
                H[k] = H4[4*i+j];
                Hd[k] = H4[8+4*i+j];
            }
        }
    } else {
        akns_scatter_U_BO(q,r,lambda,eps_t,0,E);
        akns_scatter_U_BO(q,r,lambda,0.5*eps_t,0,H);
    }
    ret_code = akns_scatter_matrix2_inverse(E,Em);
    CHECK_RETCODE(ret_code, leave_fun);
    if (derivative_flag) {
        /* (E^-1)' = -E^-1 E' E^-1. */
        akns_scatter_matrix2_mult(Em,Ed,tmp);
        akns_scatter_matrix2_mult(tmp,Em,Emd);
        for (UINT i=0; i<4; i++)
            Emd[i] = -Emd[i];
    }

    akns_scatter_matrix2_mult(Em,Dp,tmp);
    akns_scatter_matrix2_mult(tmp,E,Mp);
    akns_scatter_matrix2_mult(E,Dm,tmp);
    akns_scatter_matrix2_mult(tmp,Em,Mm);

    if (derivative_flag) {
        akns_scatter_matrix2_mult(Emd,Dp,tmp);
        akns_scatter_matrix2_mult(tmp,E,Mpd);
        akns_scatter_matrix2_mult(Em,Dp,tmp);
        akns_scatter_matrix2_mult(tmp,Ed,tmp2);
        for (UINT i=0; i<4; i++)
            Mpd[i] += tmp2[i];

        akns_scatter_matrix2_mult(Ed,Dm,tmp);
        akns_scatter_matrix2_mult(tmp,Em,Mmd);
        akns_scatter_matrix2_mult(E,Dm,tmp);
        akns_scatter_matrix2_mult(tmp,Emd,tmp2);
        for (UINT i=0; i<4; i++)
            Mmd[i] += tmp2[i];
    }

    for (UINT i=0; i<4; i++) {
        A[i] = eps_t*(Mp[i]+Mm[i])/48.0;
        if (derivative_flag)
            Ad[i] = eps_t*(Mpd[i]+Mmd[i])/48.0;
        X[i] = -A[i];
        Y[i] = A[i];
    }
    X[0] += 1.0;
    X[3] += 1.0;
    Y[0] += 1.0;
    Y[3] += 1.0;
    ret_code = akns_scatter_matrix2_inverse(X,Xinv);
    CHECK_RETCODE(ret_code, leave_fun);
    akns_scatter_matrix2_mult(Xinv,Y,C);

    akns_scatter_matrix2_mult(H,C,tmp);
    akns_scatter_matrix2_mult(tmp,H,T);
    if (derivative_flag) {
        /* C' = X^-1 A' (C+I), where X=I-A. */
        tmp[0] = C[0]+1.0;
        tmp[1] = C[1];
        tmp[2] = C[2];
        tmp[3] = C[3]+1.0;
        akns_scatter_matrix2_mult(Ad,tmp,tmp2);
        akns_scatter_matrix2_mult(Xinv,tmp2,Cd);

        akns_scatter_matrix2_mult(Hd,C,tmp);
        akns_scatter_matrix2_mult(tmp,H,Td);
        akns_scatter_matrix2_mult(H,Cd,tmp);
        akns_scatter_matrix2_mult(tmp,H,tmp2);
        for (UINT i=0; i<4; i++)
            Td[i] += tmp2[i];
        akns_scatter_matrix2_mult(H,C,tmp);
        akns_scatter_matrix2_mult(tmp,Hd,tmp2);
        for (UINT i=0; i<4; i++)
            Td[i] += tmp2[i];
    }

    if (inverse_flag) {
        ret_code = akns_scatter_matrix2_inverse(T,Tinv);
        CHECK_RETCODE(ret_code, leave_fun);
        if (derivative_flag) {
            akns_scatter_matrix2_mult(Tinv,Td,tmp);
            akns_scatter_matrix2_mult(tmp,Tinv,Tdinv);
        }
        for (UINT i=0; i<4; i++) {
            T[i] = Tinv[i];
            if (derivative_flag)
                Td[i] = -Tdinv[i];
        }
    }

    if (derivative_flag) {
        memset(U,0,16*sizeof(COMPLEX));
        U[0] = U[10] = T[0];
        U[1] = U[11] = T[1];
        U[4] = U[14] = T[2];
        U[5] = U[15] = T[3];
        U[8] = Td[0];
        U[9] = Td[1];
        U[12] = Td[2];
        U[13] = Td[3];
    } else {
        memcpy(U,T,4*sizeof(COMPLEX));
    }

leave_fun:
    return ret_code;
}

/**
 * Sixth-order exponential transition matrix exp(Z), with Z given by
 * Eqs. 65--67 of https://doi.org/10.1016/j.jcp.2021.110764. The five
 * entries of q and r are the scaled value and first through fourth
 * differences produced during preprocessing.
 */
static inline void akns_scatter_U_ES6(COMPLEX const q[5],
                                     COMPLEX const r[5],
                                     COMPLEX const lambda,
                                     REAL const eps_t,
                                     UINT const derivative_flag,
                                     UINT const inverse_flag,
                                     COMPLEX * const U)
{
    const COMPLEX z = eps_t*lambda;
    const COMPLEX z2 = z*z;
    const COMPLEX z3 = z2*z;
    COMPLEX Z[4], Zd[4] = {0};
    COMPLEX delta, delta_d = 0.0, c, s, s_delta;
    COMPLEX c_d = 0.0, s_d = 0.0;
    const REAL sign = inverse_flag ? -1.0 : 1.0;

    Z[0] = (r[0]*q[1]-q[0]*r[1])*z2/180.0
            - I*(1.0-(r[0]*q[2]+q[0]*r[2])/360.0
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

    if (derivative_flag) {
        Zd[0] = eps_t*(2.0*(r[0]*q[1]-q[0]*r[1])*z/180.0
                -I*(1.0-(r[0]*q[2]+q[0]*r[2])/360.0
                +q[1]*r[1]/60.0));
        Zd[1] = eps_t*(I*q[1]*z2/30.0-q[2]*z/90.0
                +I*(q[1]/6.0+q[3]/240.0
                -q[0]*r[0]*q[1]/90.0));
        Zd[2] = eps_t*(-I*r[1]*z2/30.0-r[2]*z/90.0
                -I*(r[1]/6.0+r[3]/240.0
                -q[0]*r[0]*r[1]/90.0));
        Zd[3] = -Zd[0];
    }

    delta = Z[0]*Z[0]+Z[1]*Z[2];
    if (CABS(delta) < 1e-8) {
        const COMPLEX delta2 = delta*delta;
        const COMPLEX delta3 = delta2*delta;
        const COMPLEX delta4 = delta3*delta;
        c = 1.0+delta/2.0+delta2/24.0+delta3/720.0
                +delta4/40320.0;
        s = 1.0+delta/6.0+delta2/120.0+delta3/5040.0
                +delta4/362880.0;
        s_delta = 1.0/6.0+delta/60.0+delta2/1680.0
                +delta3/90720.0;
    } else {
        const COMPLEX root = CSQRT(delta);
        c = CCOSH(root);
        s = misc_CSINC(I*root);
        s_delta = (c-s)/(2.0*delta);
    }
    if (derivative_flag) {
        delta_d = 2.0*Z[0]*Zd[0]+Zd[1]*Z[2]+Z[1]*Zd[2];
        c_d = 0.5*s*delta_d;
        s_d = s_delta*delta_d;
    }

    if (derivative_flag) {
        memset(U,0,16*sizeof(COMPLEX));
        U[0] = U[10] = c+sign*s*Z[0];
        U[1] = U[11] = sign*s*Z[1];
        U[4] = U[14] = sign*s*Z[2];
        U[5] = U[15] = c+sign*s*Z[3];
        U[8] = c_d+sign*(s_d*Z[0]+s*Zd[0]);
        U[9] = sign*(s_d*Z[1]+s*Zd[1]);
        U[12] = sign*(s_d*Z[2]+s*Zd[2]);
        U[13] = c_d+sign*(s_d*Z[3]+s*Zd[3]);
    } else {
        U[0] = c+sign*s*Z[0];
        U[1] = sign*s*Z[1];
        U[2] = sign*s*Z[2];
        U[3] = c+sign*s*Z[3];
    }
}

/**
 * Eighth-order exponential transition matrix exp(Z), with the degree-five
 * matrix polynomial Z from Eqs. 51--60 of arXiv:2608.11892v1. The coefficients
 * are precomputed once per time node before iterating over lambda.
 */
static inline void akns_scatter_U_ES8(COMPLEX const coeff[24],
                                     COMPLEX const lambda,
                                     REAL const eps_t,
                                     UINT const derivative_flag,
                                     UINT const inverse_flag,
                                     COMPLEX * const U)
{
    const COMPLEX z = eps_t*lambda;
    COMPLEX Z[4], Zd[4] = {0};
    COMPLEX delta, delta_d = 0.0, c, s, s_delta;
    COMPLEX c_d = 0.0, s_d = 0.0;
    const REAL sign = inverse_flag ? -1.0 : 1.0;

    for (UINT j=0; j<4; j++) {
        Z[j] = coeff[20+j];
        for (UINT k=5; k-->0; )
            Z[j] = Z[j]*z+coeff[4*k+j];
        if (derivative_flag) {
            Zd[j] = 5.0*coeff[20+j];
            for (UINT k=5; k-->1; )
                Zd[j] = Zd[j]*z+k*coeff[4*k+j];
            Zd[j] *= eps_t;
        }
    }

    delta = Z[0]*Z[0]+Z[1]*Z[2];
    if (CABS(delta) < 1e-8) {
        const COMPLEX delta2 = delta*delta;
        const COMPLEX delta3 = delta2*delta;
        const COMPLEX delta4 = delta3*delta;
        c = 1.0+delta/2.0+delta2/24.0+delta3/720.0
                +delta4/40320.0;
        s = 1.0+delta/6.0+delta2/120.0+delta3/5040.0
                +delta4/362880.0;
        s_delta = 1.0/6.0+delta/60.0+delta2/1680.0
                +delta3/90720.0;
    } else {
        const COMPLEX root = CSQRT(delta);
        c = CCOSH(root);
        s = misc_CSINC(I*root);
        s_delta = (c-s)/(2.0*delta);
    }
    if (derivative_flag) {
        delta_d = 2.0*Z[0]*Zd[0]+Zd[1]*Z[2]+Z[1]*Zd[2];
        c_d = 0.5*s*delta_d;
        s_d = s_delta*delta_d;
        memset(U,0,16*sizeof(COMPLEX));
        U[0] = U[10] = c+sign*s*Z[0];
        U[1] = U[11] = sign*s*Z[1];
        U[4] = U[14] = sign*s*Z[2];
        U[5] = U[15] = c+sign*s*Z[3];
        U[8] = c_d+sign*(s_d*Z[0]+s*Zd[0]);
        U[9] = sign*(s_d*Z[1]+s*Zd[1]);
        U[12] = sign*(s_d*Z[2]+s*Zd[2]);
        U[13] = c_d+sign*(s_d*Z[3]+s*Zd[3]);
    } else {
        U[0] = c+sign*s*Z[0];
        U[1] = sign*s*Z[1];
        U[2] = sign*s*Z[2];
        U[3] = c+sign*s*Z[3];
    }
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
                        INT * const W,
                        akns_discretization_t const discretization,
                        akns_pde_t const PDE,
                        UINT const vanilla_flag,
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
    if (D%upsampling_factor != 0)
        return E_ASSERTION_FAILED;

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *eps_t_scaled = NULL;

    // Define stepsize constants that are often needed
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.

    UINT N = 0;
    switch (discretization) {
        case akns_discretization_CT4:
        case akns_discretization_ES6:
            break;
        case akns_discretization_ES8:
            tmp1 = malloc(24*(D/7)*sizeof(COMPLEX));
            CHECK_NOMEM(tmp1,ret_code,leave_fun);
            for (UINT n=0, node=0; n<D; n+=7, node++)
                fnft__akns_es8_z_coefficients(&q[n],&r[n],&tmp1[24*node]);
            break;
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
            COMPLEX *qr_weights = NULL;
            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_scaled,discretization);
            CHECK_RETCODE(ret_code, leave_fun_no_eps_t_scaled); // if ret_code != SUCCESS, akns_discretization_method_weights frees qr_weights and eps_t_scaled if needed
            free(qr_weights);
            for (UINT n=0; n<upsampling_factor; n++ )
                eps_t_scaled[n] *= eps_t;
            break;
            
        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(>discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

    if (derivative_flag){
        // Calculate the scattering matrix with lamda-derivative as in G. Boffetta an A.R. Osborne, 'Computation of the direct scattering transform for the nonlinear Schroedinger equation', www.doi.org/10.1016/0021-9991(92)90370-e .
        COMPLEX U[4][4] = {{ 0 }};
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][4][4] = { { {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1} } }; // Initiate only first sixteen values
            UINT current = 0;
            INT Wi = 0;

            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    for (UINT n = 0; n < D; n++){
                        COMPLEX h = eps_t_scaled[n%upsampling_factor];
                        akns_scatter_U_BO(q[n],r[n],l_curr,h,1,*U);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(16, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;

                case akns_discretization_ES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                        akns_scatter_U_ES4(a1,a2,a3,1,*U,&tmp2[n]);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(16, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_TES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX M[2][2];

                        // First substep, block diagonal matrix
                        akns_scatter_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*M,NULL);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][0][0],&H[!current][0][0]);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][2][0],&H[!current][2][0]);
                        current = !current;

                        // Second substep
                        akns_scatter_U_BO(q[n],r[n],l_curr,eps_t,1,*U);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        // Third substep, block diagonal matrix
                        akns_scatter_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*M,NULL);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][0][0],&H[!current][0][0]);
                        misc_matrix_mult(2,2,4,&M[0][0],&H[current][2][0],&H[!current][2][0]);
                        current = !current;

                        if (W != NULL)
                            Wi += misc_normalize_vector(16, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_CT4:
                    for (UINT n = 0; n < D; n+=3) {
                        ret_code = akns_scatter_U_CT4(q[n],r[n],q[n+1],r[n+1],
                                q[n+2],r[n+2],l_curr,eps_t,1,0,*U);
                        CHECK_RETCODE(ret_code, leave_fun);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(16, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_ES6:
                    for (UINT n = 0; n < D; n+=5) {
                        akns_scatter_U_ES6(&q[n],&r[n],l_curr,eps_t,1,0,*U);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(16, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_ES8:
                    for (UINT n = 0, node=0; n < D; n+=7, node++) {
                        akns_scatter_U_ES8(&tmp1[24*node],l_curr,eps_t,1,0,*U);
                        misc_matrix_mult(4,4,4,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(16,&H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;

                default: // Unknown discretization
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
            COMPLEX Tmx[4][4];

            // Fetch the change of basis matrix from the basis of the discretization to S.
            ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Left-multiply the change of state matrix by this state of basis matrix.
            misc_matrix_mult(4,4,4,&Tmx[0][0],&H[current][0][0],&H[!current][0][0]);
            current = !current;

            // Fetch the change of basis matrix from the basis of the discretization to S.
            ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Right-multiply the change of state matrix by this state of basis matrix to obtain the scattering matrix in S basis.
            misc_matrix_mult(4,4,4,&H[current][0][0],&Tmx[0][0],&H[!current][0][0]);
            current = !current;

            // Copy result
            result[8*i + 0] = 0.5*H[current][0][0] + 0.5*H[current][2][2];
            result[8*i + 1] = 0.5*H[current][0][1] + 0.5*H[current][2][3];
            result[8*i + 2] = 0.5*H[current][1][0] + 0.5*H[current][3][2];
            result[8*i + 3] = 0.5*H[current][1][1] + 0.5*H[current][3][3];
            result[8*i + 4] = H[current][2][0];
            result[8*i + 5] = H[current][2][1];
            result[8*i + 6] = H[current][3][0];
            result[8*i + 7] = H[current][3][1];
        }
    } else {
        // Calculate the scattering matrix without lambda-derivative
        for (UINT i = 0; i < K; i++) { // iterate over lambda
            // Initialize scattering matrix
            COMPLEX l_curr = lambda[i];
            COMPLEX H[2][2][2] = { { {1,0}, {0,1} } }; // Initiate only first four values
            UINT current = 0;
            INT Wi = 0;

            switch (discretization) {
                case akns_discretization_BO:
                case akns_discretization_CF4_2:
                case akns_discretization_CF4_3:
                case akns_discretization_CF5_3:
                case akns_discretization_CF6_4:
                    for (UINT n = 0; n < D; n++){
                        COMPLEX U[2][2];
                        COMPLEX h = eps_t_scaled[n%upsampling_factor];
                        akns_scatter_U_BO(q[n],r[n],l_curr,h,0,*U);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(4, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;

                case akns_discretization_ES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX U[2][2];

                        COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                        akns_scatter_U_ES4(a1,a2,a3,0,*U,NULL);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(4, &H[current][0][0]);                      
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;

                case akns_discretization_TES4:
                    for (UINT n = 0; n < D; n+=3){
                        COMPLEX U[2][2];

                        // First substep
                        akns_scatter_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*U,NULL);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        // Second substep
                        akns_scatter_U_BO(q[n],r[n],l_curr,eps_t,0,*U);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        // Third substep
                        akns_scatter_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*U,NULL);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;

                        if (W != NULL)
                            Wi += misc_normalize_vector(4, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_CT4:
                    for (UINT n = 0; n < D; n+=3) {
                        COMPLEX U[2][2] = {{0}};
                        ret_code = akns_scatter_U_CT4(q[n],r[n],q[n+1],r[n+1],
                                q[n+2],r[n+2],l_curr,eps_t,0,0,*U);
                        CHECK_RETCODE(ret_code, leave_fun);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(4, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_ES6:
                    for (UINT n = 0; n < D; n+=5) {
                        COMPLEX U[2][2] = {{0}};
                        akns_scatter_U_ES6(&q[n],&r[n],l_curr,eps_t,0,0,*U);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(4, &H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;
                case akns_discretization_ES8:
                    for (UINT n = 0, node=0; n < D; n+=7, node++) {
                        COMPLEX U[2][2] = {{0}};
                        akns_scatter_U_ES8(&tmp1[24*node],l_curr,eps_t,0,0,*U);
                        misc_matrix_mult(2,2,2,&U[0][0],&H[current][0][0],&H[!current][0][0]);
                        current = !current;
                        if (W != NULL)
                            Wi += misc_normalize_vector(4,&H[current][0][0]);
                    }
                    if (W != NULL)
                        W[i] = Wi;
                    break;

                default: // Unknown discretization
                    ret_code = E_INVALID_ARGUMENT(discretization);
                    CHECK_RETCODE(ret_code, leave_fun);
            }
            COMPLEX Tmx[2][2];

            // Fetch the change of basis matrix from the basis of the discretization to S.
            ret_code = akns_discretization_change_of_basis_matrix_to_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Left-multiply the change of state matrix by this state of basis matrix.
            misc_matrix_mult(2,2,2,&Tmx[0][0],&H[current][0][0],&H[!current][0][0]);
            current = !current;

            // Fetch the change of basis matrix from the basis of the discretization to S.
            ret_code = akns_discretization_change_of_basis_matrix_from_S(&Tmx[0][0],l_curr,derivative_flag,eps_t,discretization,vanilla_flag,PDE);
            CHECK_RETCODE(ret_code, leave_fun);

            // Right-multiply the change of state matrix by this state of basis matrix to obtain the scattering matrix in S basis.
            misc_matrix_mult(2,2,2,&H[current][0][0],&Tmx[0][0],&H[!current][0][0]);
            current = !current;

            // Copy result
            memcpy(&result[4*i],&H[current][0][0],4 * sizeof(COMPLEX));
        }
    }
    
leave_fun:
    free(eps_t_scaled);
leave_fun_no_eps_t_scaled:
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
                              COMPLEX * const b_vals,
                              INT * const Ws,
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
    if (!skip_b_flag && b_vals == NULL)
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
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp4 = NULL, *eps_t_scaled = NULL;

    INT * WPHI = NULL; // for storing intermediate scaling factors 
    INT * WPSI = NULL; // if ntormalization is enabled
                       
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
    CHECK_NOMEM(PSIPHI, ret_code, leave_fun);
    COMPLEX * const PSI = &PSIPHI[0];
    COMPLEX * const PHI = &PSIPHI[2];

    // We need to store many intermediate scaling factors for the
    // forward-backward computation of b if normalization is on
    const INT normalization_flag = Ws != NULL;
    if (normalization_flag && !skip_b_flag) {
        WPHI = calloc((D_given + 1), sizeof(COMPLEX)); // calloc initializes to zero
        CHECK_NOMEM(WPHI, ret_code, leave_fun);
        WPSI = calloc((D_given + 1), sizeof(COMPLEX)); // calloc initializes to zero
        CHECK_NOMEM(WPHI, ret_code, leave_fun);
    }

    // Define stepsize constants that are often needed
    REAL const eps_t = (T[1] - T[0])/(D_given - 1);
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // In the case of ES4 and TES4 computing values that are functions of
    // q and r but not l.

    UINT N = 0;
    switch (discretization) {
        case akns_discretization_CT4:
        case akns_discretization_ES6:
            break;
        case akns_discretization_ES8:
            tmp1 = malloc(24*(D/7)*sizeof(COMPLEX));
            CHECK_NOMEM(tmp1,ret_code,leave_fun);
            for (UINT n=0, node=0; n<D; n+=7, node++)
                fnft__akns_es8_z_coefficients(&q[n],&r[n],&tmp1[24*node]);
            break;
        // Fourth-order exponential method which requires
        // one matrix exponential. The matrix exponential is
        // implmented by using the expansion of the 2x2 matrix
        // in terms of Pauli matrices.
        case akns_discretization_ES4:
            tmp1 = malloc(2*D*sizeof(COMPLEX));
            CHECK_NOMEM(tmp1, ret_code, leave_fun);
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
            CHECK_NOMEM(tmp1, ret_code, leave_fun);
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
                CHECK_NOMEM(tmp3, ret_code, leave_fun);
                CHECK_NOMEM(tmp4, ret_code, leave_fun);
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
            COMPLEX *qr_weights = NULL;
            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_scaled,discretization);
            CHECK_RETCODE(ret_code, leave_fun_no_eps_t_scaled); // if ret_code != SUCCESS, akns_discretization_method_weights frees qr_weights and eps_t_scaled if needed
            free(qr_weights);
            for (UINT n=0; n<upsampling_factor; n++ )
                eps_t_scaled[n] *= eps_t;
            break;

        default: // Unknown discretization
            ret_code = E_INVALID_ARGUMENT(discretization);
            CHECK_RETCODE(ret_code, leave_fun);
    }

    for (UINT neig=0; neig<K; neig++) { // iterate over bound states
        COMPLEX l_curr = bound_states[neig];
        INT WPHI_acc = 0; // accumulated scaling factor for phi
        INT WPSI_acc = 0; // ... for psi

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

        // Declaring change of state matrix here, to avoid letting them be
        // overwritten with zeros in every loop iteration.
        COMPLEX U[4][4] = {{0}};
        switch (discretization) {

            case akns_discretization_BO:
            case akns_discretization_CF4_2:
            case akns_discretization_CF4_3:
            case akns_discretization_CF5_3:
            case akns_discretization_CF6_4:
            {
                COMPLEX phi_temp[2][4];
                UINT current = 0;
                memcpy(&phi_temp[current][0], PHI, 4 * sizeof(COMPLEX));
                for (UINT n_given=0; n_given<D_given; n_given++) {
                    for (UINT count=0; count<upsampling_factor; count++) {
                        UINT n = n_given * upsampling_factor + count;
                        akns_scatter_U_BO(q[n],r[n],l_curr,eps_t_scaled[count],1,*U);
                        misc_matrix_mult(4,4,1,&U[0][0],&phi_temp[current][0],&phi_temp[!current][0]);
                        current = !current;
                    }
                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4, &phi_temp[current][0]);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }
                    memcpy(&PHI[4*(n_given+1)], &phi_temp[current][0], 4 * sizeof(COMPLEX));
                }
            }
            break;

            // Fourth-order exponential method which requires
            // one matrix exponential. The matrix exponential is
            // implmented by using the expansion of the 2x2 matrix
            // in terms of Pauli matrices.
            case akns_discretization_ES4:
                for (UINT n = 0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX a1 = tmp1[n]+ eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                    COMPLEX a2 = tmp1[n+1] - eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                    COMPLEX a3 = - eps_t*I*l_curr +tmp1[n+2];
                    akns_scatter_U_ES4(a1,a2,a3,1,*U,&tmp2[n]);                   
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4, &PHI[4*(n_given+1)]);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }                  
                }
                break;

            // Fourth-order exponential method which requires
            // three matrix exponentials. The outer two transfer metrices
            // need to be built differently compared to the CF schemes.
            case akns_discretization_TES4:
                for (UINT n=0, n_given=0; n<D; n+=3, n_given++) {
                    COMPLEX phi_temp[4], M[2][2];

                    // First substep, block diagonal matrix
                    akns_scatter_U_ES4(tmp1[n],tmp1[n+1],0.0,0,*M,NULL);
                    misc_matrix_mult(2,2,1,*M,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    misc_matrix_mult(2,2,1,*M,&PHI[4*n_given+2],&PHI[4*(n_given+1)+2]);

                    // Second substep
                    akns_scatter_U_BO(q[n],r[n],l_curr,eps_t,1,*U);
                    misc_matrix_mult(4,4,1,*U,&PHI[4*(n_given+1)],phi_temp);

                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4, phi_temp);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }

                    // Third substep, block diagonal matrix
                    akns_scatter_U_ES4(tmp2[n],tmp2[n+1],0.0,0,*M,NULL);
                    misc_matrix_mult(2,2,1,*M,&phi_temp[0],&PHI[4*(n_given+1)]);
                    misc_matrix_mult(2,2,1,*M,&phi_temp[2],&PHI[4*(n_given+1)+2]);
                }
                break;

            case akns_discretization_CT4:
                for (UINT n=0, n_given=0; n<D; n+=3, n_given++) {
                    ret_code = akns_scatter_U_CT4(q[n],r[n],q[n+1],r[n+1],
                            q[n+2],r[n+2],l_curr,eps_t,1,0,*U);
                    CHECK_RETCODE(ret_code, leave_fun);
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4, &PHI[4*(n_given+1)]);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }
                }
                break;
            case akns_discretization_ES6:
                for (UINT n=0, n_given=0; n<D; n+=5, n_given++) {
                    akns_scatter_U_ES6(&q[n],&r[n],l_curr,eps_t,1,0,*U);
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4, &PHI[4*(n_given+1)]);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }
                }
                break;
            case akns_discretization_ES8:
                for (UINT n=0, n_given=0; n<D; n+=7, n_given++) {
                    akns_scatter_U_ES8(&tmp1[24*n_given],l_curr,eps_t,1,0,*U);
                    misc_matrix_mult(4,4,1,*U,&PHI[4*n_given],&PHI[4*(n_given+1)]);
                    if (normalization_flag) {
                        WPHI_acc += misc_normalize_vector(4,&PHI[4*(n_given+1)]);
                        if (WPHI != NULL)
                            WPHI[n_given+1] = WPHI_acc;
                    }
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
        if (!skip_b_flag) {

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
                {
                    COMPLEX psi_temp[2][2];
                    UINT current = 0;
                    memcpy(&psi_temp[current][0], &PSI[4*D_given], 2 * sizeof(COMPLEX));
                    for (UINT n_given=D_given; n_given-->0; ) {
                        for (UINT count=upsampling_factor; count-->0; ) {
                            COMPLEX U[2][2];
                            UINT n = n_given * upsampling_factor + count;
                            akns_scatter_U_BO(q[n],r[n],l_curr,-eps_t_scaled[count],0,*U);
                            misc_matrix_mult(2,2,1,&U[0][0],&psi_temp[current][0],&psi_temp[!current][0]);
                            current = !current;
                        }
                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2, &psi_temp[current][0]);
                            WPSI[n_given] = WPSI_acc;
                        }             
                        memcpy(&PSI[4*n_given], &psi_temp[current][0], 2 * sizeof(COMPLEX));
                    }
                    break;
                }

                // Fourth-order exponential method which requires
                // one matrix exponential. The matrix exponential is
                // implmented by using the expansion of the 2x2 matrix
                // in terms of Pauli matrices.
                case akns_discretization_ES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2];
                        COMPLEX a1 = -tmp1[n]- eps_t_3*(l_curr*I*(q[n+1]-r[n+1]))/12.0;
                        COMPLEX a2 = -tmp1[n+1] + eps_t_3*l_curr*(q[n+1]+r[n+1])/12.0;
                        COMPLEX a3 =  eps_t*I*l_curr -tmp1[n+2];
                        akns_scatter_U_ES4(a1,a2,a3,0,*U,NULL);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2, &PSI[4*n_given]);
                            WPSI[n_given] = WPSI_acc;
                        }
                    }
                    break;

                // Fourth-order exponential method which requires
                // three matrix exponentials. The transfer metrix cannot
                // needs to be built differently compared to the CF schemes.
                case akns_discretization_TES4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2], psi_temp[2];

                        // First substep
                        akns_scatter_U_ES4(tmp3[n],tmp3[n+1],0.0,0,*U,NULL);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);

                        // Second substep
                        akns_scatter_U_BO(q[n],r[n],l_curr,-eps_t,0,*U);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*n_given],psi_temp);

                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2, psi_temp);
                            WPSI[n_given] = WPSI_acc;
                        }

                        // Third substep
                        akns_scatter_U_ES4(tmp4[n],tmp4[n+1],0.0,0,*U,NULL);
                        misc_matrix_mult(2,2,1,*U,psi_temp,&PSI[4*n_given]);
                    }
                    break;

                case akns_discretization_CT4:
                    for (UINT n_given=D_given, n=D-3; n_given-->0; n-=3) {
                        COMPLEX U[2][2] = {{0}};
                        ret_code = akns_scatter_U_CT4(q[n],r[n],q[n+1],r[n+1],
                                q[n+2],r[n+2],l_curr,eps_t,0,1,*U);
                        CHECK_RETCODE(ret_code, leave_fun);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2, &PSI[4*n_given]);
                            WPSI[n_given] = WPSI_acc;
                        }
                    }
                    break;
                case akns_discretization_ES6:
                    for (UINT n_given=D_given, n=D-5; n_given-->0; n-=5) {
                        COMPLEX U[2][2] = {{0}};
                        akns_scatter_U_ES6(&q[n],&r[n],l_curr,eps_t,0,1,*U);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2, &PSI[4*n_given]);
                            WPSI[n_given] = WPSI_acc;
                        }
                    }
                    break;
                case akns_discretization_ES8:
                    for (UINT n_given=D_given, n=D-7; n_given-->0; n-=7) {
                        COMPLEX U[2][2] = {{0}};
                        akns_scatter_U_ES8(&tmp1[24*n_given],l_curr,eps_t,0,1,*U);
                        misc_matrix_mult(2,2,1,*U,&PSI[4*(n_given+1)],&PSI[4*n_given]);
                        if (normalization_flag) {
                            WPSI_acc += misc_normalize_vector(2,&PSI[4*n_given]);
                            WPSI[n_given] = WPSI_acc;
                        }
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
        const COMPLEX exponent = I*l_curr*(T[1]+eps_t*boundary_coeff);
        if (!normalization_flag) {
            a_vals[neig] = f_S[0]*CEXP(exponent);
            aprime_vals[neig] = f_S[2]*CEXP(exponent) + I*(T[1]+eps_t*boundary_coeff)*a_vals[neig];
        } else {
            const REAL l2 = LOG(2); // to mix powers of 2 and e
            // We'll have to multiply with CEXP(exponent), which can drastically
            // change the absolute values of the results. Therefore, we write the
            // amplitude of that term as 2^u*2^v, with -1/2<=u<=1/2 and v integer,
            // and we shift the 2^v part into the scaling factor 2^Ws.
            const INT v = ROUND(CREAL(exponent)/l2);
            const COMPLEX u = CREAL(exponent)/l2 - v; // between -1/2 and 1/2
            a_vals[neig] = f_S[0] * CEXP(I*CIMAG(exponent) + u*l2);
            Ws[neig] = v + WPHI_acc;
            aprime_vals[neig] = f_S[2] + I*(T[1]+eps_t*boundary_coeff)*f_S[0];
            aprime_vals[neig] *= CEXP(I*CIMAG(exponent) + u*l2);
        }
        if (PDE==akns_pde_KdV) {
            a_vals[neig] = CREAL(a_vals[neig]);
            aprime_vals[neig] = I * CIMAG(aprime_vals[neig]);
        }

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
                for (UINT i=0; i<2; i++) {
                    b_temp[i] = phi_S[i]/psi_S[i];
                    if (normalization_flag) {
                        b_temp[i] *= POW(2, WPHI[n] - WPSI[n]);
                    }
                }
                if (PDE!=akns_pde_KdV || CREAL(b_temp[0]*b_temp[1])>0) {
                    tmp = FABS( 0.5* LOG( (REAL)CABS( b_temp[1]/b_temp[0] ) ) );
                    if (tmp < error_metric){
                        b_vals[neig] = b_temp[0];
                        error_metric = tmp;
                    }
                }
            }
        }
    }

leave_fun:
    free(eps_t_scaled);
leave_fun_no_eps_t_scaled:
    free(tmp1);
    free(PSIPHI);
    free(WPHI);
    free(WPSI);
    return ret_code;
}
