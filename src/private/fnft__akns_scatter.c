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
 * Igor Chekhovskoy (NSU, FRC ICT) 2026.
 * Irina Vaseva (FRC ICT, NSU) 2026.
 */
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__akns_scatter.h"
#include "fnft__akns_fscatter_pade.h"

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
 * Exponential transition matrix exp(Z) for ES4, ES6 and ES8. The
 * order-dependent coefficients of Z(eps_t*lambda) are precomputed once
 * per time node. This routine also evaluates the lambda derivative and
 * the inverse transition matrix from the same polynomial.
 */
typedef struct {
    UINT degree;
    REAL const *F;
    REAL const *G;
    REAL const *denominator;
} akns_scatter_pade_coefficients_t;

/* The common slow-ES state already lies in the no-op interval of
 * misc_normalize_vector. Avoid cabs, log2 and pow in that case. */
static inline INT akns_scatter_normalize_ES_vector(UINT const len,
                                                    COMPLEX * const v)
{
    REAL max_abs_squared = 0.0;

    for (UINT i=0; i<len; i++) {
        const REAL real_part = CREAL(v[i]);
        const REAL imag_part = CIMAG(v[i]);
        const REAL abs_squared = real_part*real_part+imag_part*imag_part;

        if (!isfinite(abs_squared))
            return misc_normalize_vector(len,v);
        if (abs_squared > max_abs_squared)
            max_abs_squared = abs_squared;
    }
    if (max_abs_squared > 1.0+8.0*EPSILON
            && max_abs_squared < 4.0-32.0*EPSILON)
        return 0;
    return misc_normalize_vector(len,v);
}

/* Some C toolchains lower every complex product to a runtime call.
 * Keep finite products inline without enabling unsafe global math flags. */
static inline COMPLEX akns_scatter_ES_multiply(COMPLEX const a,
                                                COMPLEX const b)
{
#if defined(__MINGW32__) && defined(__GNUC__) && !defined(__clang__)
    const REAL ar = CREAL(a), ai = CIMAG(a);
    const REAL br = CREAL(b), bi = CIMAG(b);
    const REAL parts[2] = {ar*br-ai*bi,ar*bi+ai*br};
    COMPLEX value;

    if (!isfinite(parts[0]) || !isfinite(parts[1]))
        return a*b;
    memcpy(&value,parts,sizeof(value));
    return value;
#else
    return a*b;
#endif
}

static inline void akns_scatter_apply_U_ES(
        COMPLEX const * const U, UINT const U_stride,
        COMPLEX const * const v, UINT const v_stride,
        COMPLEX * const result, UINT const result_stride)
{
    result[0] = akns_scatter_ES_multiply(U[0],v[0])
            +akns_scatter_ES_multiply(U[1],v[v_stride]);
    result[result_stride] = akns_scatter_ES_multiply(U[U_stride],v[0])
            +akns_scatter_ES_multiply(U[U_stride+1],v[v_stride]);
}

static inline void akns_scatter_apply_U_ES_derivative(
        COMPLEX const * const U, COMPLEX const v[4], COMPLEX result[4])
{
    COMPLEX derivative_part[2];

    akns_scatter_apply_U_ES(U,4,v,1,result,1);
    akns_scatter_apply_U_ES(&U[8],4,v,1,derivative_part,1);
    akns_scatter_apply_U_ES(&U[10],4,&v[2],1,&result[2],1);
    result[2] += derivative_part[0];
    result[3] += derivative_part[1];
}

enum { AKNS_SCATTER_ES_BLOCK_SIZE = 32 };

/* Independent local transitions are built in small blocks. The real and
 * imaginary lanes let the compiler vectorize their polynomial arithmetic;
 * propagation and normalization still occur in the original node order. */
static inline void akns_scatter_ES_polynomial_block(
        COMPLEX const * const coeff, UINT const degree, COMPLEX const z,
        UINT const count, COMPLEX result[AKNS_SCATTER_ES_BLOCK_SIZE][4])
{
    REAL re[3][AKNS_SCATTER_ES_BLOCK_SIZE];
    REAL im[3][AKNS_SCATTER_ES_BLOCK_SIZE];
    const UINT stride = 4*(degree+1);
    const REAL zr = CREAL(z), zi = CIMAG(z);

    /* A linear polynomial needs no intermediate Horner arrays. */
    if (degree == 1 && zi == 0.0 && isfinite(zr)) {
        INT finite = 1;
        for (UINT n=0; n<count; n++) {
            for (UINT j=0; j<3; j++) {
                const REAL parts[2] = {
                    CREAL(coeff[stride*n+4+j])*zr+CREAL(coeff[stride*n+j]),
                    CIMAG(coeff[stride*n+4+j])*zr+CIMAG(coeff[stride*n+j])
                };
                memcpy(&result[n][j],parts,sizeof(COMPLEX));
                finite &= isfinite(parts[0]) & isfinite(parts[1]);
            }
        }
        if (!finite) {
            for (UINT n=0; n<count; n++) {
                for (UINT j=0; j<3; j++) {
                    if (!isfinite(CREAL(result[n][j]))
                            || !isfinite(CIMAG(result[n][j])))
                        result[n][j] = coeff[stride*n+4+j]*z+coeff[stride*n+j];
                }
            }
        }
        for (UINT n=0; n<count; n++)
            result[n][3] = -result[n][0];
        return;
    }
    for (UINT j=0; j<3; j++) {
        for (UINT n=0; n<count; n++) {
            re[j][n] = CREAL(coeff[stride*n+4*degree+j]);
            im[j][n] = CIMAG(coeff[stride*n+4*degree+j]);
        }
    }
    if (zi == 0.0 && isfinite(zr)) {
        for (UINT k=degree; k-->0; ) {
            for (UINT j=0; j<3; j++) {
                for (UINT n=0; n<count; n++) {
                    re[j][n] = re[j][n]*zr+CREAL(coeff[stride*n+4*k+j]);
                    im[j][n] = im[j][n]*zr+CIMAG(coeff[stride*n+4*k+j]);
                }
            }
        }
    } else {
        for (UINT k=degree; k-->0; ) {
            for (UINT j=0; j<3; j++) {
                for (UINT n=0; n<count; n++) {
                    const REAL r = re[j][n], i = im[j][n];
                    re[j][n] = r*zr-i*zi+CREAL(coeff[stride*n+4*k+j]);
                    im[j][n] = r*zi+i*zr+CIMAG(coeff[stride*n+4*k+j]);
                }
            }
        }
    }
    INT finite = 1;
    for (UINT j=0; j<3; j++) {
        for (UINT n=0; n<count; n++) {
            const REAL parts[2] = {re[j][n],im[j][n]};
            memcpy(&result[n][j],parts,sizeof(COMPLEX));
            finite &= isfinite(parts[0]) & isfinite(parts[1]);
        }
    }
    if (!finite) {
        for (UINT n=0; n<count; n++) {
            for (UINT j=0; j<3; j++) {
                if (!isfinite(re[j][n]) || !isfinite(im[j][n])) {
                    /* Retain C complex arithmetic at exceptional magnitudes. */
                    result[n][j] = coeff[stride*n+4*degree+j];
                    for (UINT k=degree; k-->0; )
                        result[n][j] = result[n][j]*z+coeff[stride*n+4*k+j];
                }
            }
        }
    }
    for (UINT n=0; n<count; n++)
        result[n][3] = -result[n][0];
}

static inline void akns_scatter_ES_pade_block(COMPLEX const * const delta,
        UINT const count, akns_scatter_pade_coefficients_t const * const pade,
        COMPLEX result[AKNS_SCATTER_ES_BLOCK_SIZE][3])
{
    REAL re[3][AKNS_SCATTER_ES_BLOCK_SIZE];
    REAL im[3][AKNS_SCATTER_ES_BLOCK_SIZE];
    REAL dr[AKNS_SCATTER_ES_BLOCK_SIZE], di[AKNS_SCATTER_ES_BLOCK_SIZE];
    const UINT degree = pade->degree;
    REAL const * const coefficients[3] = {pade->F,pade->G,pade->denominator};

    for (UINT n=0; n<count; n++) {
        dr[n] = CREAL(delta[n]);
        di[n] = CIMAG(delta[n]);
    }
    for (UINT j=0; j<3; j++) {
        const UINT d = j == 1 ? degree-1 : degree;
        for (UINT n=0; n<count; n++) {
            re[j][n] = coefficients[j][d];
            im[j][n] = 0.0;
        }
        for (UINT k=d; k-->0; ) {
            for (UINT n=0; n<count; n++) {
                const REAL r = re[j][n], i = im[j][n];
                re[j][n] = r*dr[n]-i*di[n]+coefficients[j][k];
                im[j][n] = r*di[n]+i*dr[n];
            }
        }
        INT finite = 1;
        for (UINT n=0; n<count; n++) {
            const REAL parts[2] = {re[j][n],im[j][n]};
            memcpy(&result[n][j],parts,sizeof(COMPLEX));
            finite &= isfinite(parts[0]) & isfinite(parts[1]);
        }
        if (!finite) {
            for (UINT n=0; n<count; n++) {
                if (!isfinite(re[j][n]) || !isfinite(im[j][n])) {
                    result[n][j] = coefficients[j][d];
                    for (UINT k=d; k-->0; )
                        result[n][j] = result[n][j]*delta[n]+coefficients[j][k];
                }
            }
        }
    }
}

static inline INT akns_scatter_U_ES(COMPLEX const * const coeff,
                                     UINT const degree,
                                     COMPLEX const lambda,
                                     REAL const eps_t,
                                     akns_scatter_pade_coefficients_t const *
                                        const pade,
                                     UINT const derivative_flag,
                                     UINT const inverse_flag,
                                     COMPLEX * const result, UINT const count)
{
    const COMPLEX z = eps_t*lambda;
    COMPLEX Z_values[AKNS_SCATTER_ES_BLOCK_SIZE][4];
    COMPLEX deltas[AKNS_SCATTER_ES_BLOCK_SIZE];
    COMPLEX pade_values[AKNS_SCATTER_ES_BLOCK_SIZE][3];
    const REAL sign = inverse_flag ? -1.0 : 1.0;
    INT ret_code = SUCCESS;

    akns_scatter_ES_polynomial_block(coeff,degree,z,count,Z_values);
    for (UINT n=0; n<count; n++) {
        COMPLEX const * const Z = Z_values[n];
        deltas[n] = akns_scatter_ES_multiply(Z[0],Z[0])
                +akns_scatter_ES_multiply(Z[1],Z[2]);
    }
    if (pade != NULL)
        akns_scatter_ES_pade_block(deltas,count,pade,pade_values);

    for (UINT n=0; n<count; n++) {
        COMPLEX const * const Z = Z_values[n];
        COMPLEX const * const node_coeff = coeff+4*(degree+1)*n;
        COMPLEX * const U = result+(derivative_flag ? 16 : 4)*n;
        const COMPLEX delta = deltas[n];
        COMPLEX Zd[4] = {0}, f, g, f_d = 0.0, g_d = 0.0;
        COMPLEX delta_d = 0.0;
        if (derivative_flag) {
            for (UINT j=0; j<3; j++) {
                Zd[j] = degree*node_coeff[4*degree+j];
                for (UINT k=degree; k-->1; )
                    Zd[j] = Zd[j]*z+k*node_coeff[4*k+j];
                Zd[j] *= eps_t;
            }
            Zd[3] = -Zd[0];
            delta_d = 2.0*Z[0]*Zd[0]+Zd[1]*Z[2]+Z[1]*Zd[2];
        }
        if (pade == NULL) {
            COMPLEX f_delta;

            if (CABS(delta) < 1e-8) {
                const COMPLEX delta2 = delta*delta;
                const COMPLEX delta3 = delta2*delta;
                const COMPLEX delta4 = delta3*delta;
                f = 1.0+delta/2.0+delta2/24.0+delta3/720.0
                        +delta4/40320.0;
                g = 1.0+delta/6.0+delta2/120.0+delta3/5040.0
                        +delta4/362880.0;
                f_delta = 1.0/6.0+delta/60.0+delta2/1680.0
                        +delta3/90720.0;
            } else {
                const COMPLEX root = CSQRT(delta);
                fnft__akns_scatter_exact_scalar_pair(root,&f,&g);
                if (derivative_flag)
                    f_delta = (f-g)/(2.0*delta);
            }
            if (derivative_flag) {
                f_d = 0.5*g*delta_d;
                g_d = f_delta*delta_d;
            }
        } else {
            const UINT pade_degree = pade->degree;
            REAL const * const F_coeff = pade->F;
            REAL const * const G_coeff = pade->G;
            REAL const * const denominator_coeff = pade->denominator;
            COMPLEX F, G, denominator, denominator_inverse;
            COMPLEX F_delta = 0.0, G_delta = 0.0;
            COMPLEX denominator_delta = 0.0;
            F = pade_values[n][0];
            G = pade_values[n][1];
            denominator = pade_values[n][2];
            if (denominator == 0.0) {
                ret_code = E_DIV_BY_ZERO;
                goto leave_fun;
            }
            {
                const REAL ar = CREAL(denominator), ai = CIMAG(denominator);
                const REAL square = ar*ar+ai*ai;
                if (square >= DBL_MIN && square <= 1.0/DBL_MIN) {
                    const REAL inv = 1.0/square;
                    const REAL parts[2] = {ar*inv,-ai*inv};
                    memcpy(&denominator_inverse,parts,sizeof(denominator_inverse));
                } else {
                    denominator_inverse = 1.0/denominator;
                }
            }
            if (derivative_flag) {
                F_delta = pade_degree*F_coeff[pade_degree];
                denominator_delta = pade_degree
                        *denominator_coeff[pade_degree];
                for (UINT k=pade_degree; k-->1; ) {
                    F_delta = F_delta*delta+k*F_coeff[k];
                    denominator_delta = denominator_delta*delta
                            +k*denominator_coeff[k];
                }
                if (pade_degree > 1) {
                    G_delta = (pade_degree-1)*G_coeff[pade_degree-1];
                    for (UINT k=pade_degree-1; k-->1; )
                        G_delta = G_delta*delta+k*G_coeff[k];
                }
                f_d = (F_delta*denominator-F*denominator_delta)*delta_d
                        *denominator_inverse*denominator_inverse;
                g_d = (G_delta*denominator-G*denominator_delta)*delta_d
                        *denominator_inverse*denominator_inverse;
            }
            f = akns_scatter_ES_multiply(F,denominator_inverse);
            g = akns_scatter_ES_multiply(G,denominator_inverse);
        }
        if (derivative_flag) {
            memset(U,0,16*sizeof(COMPLEX));
            U[0] = U[10] = f+akns_scatter_ES_multiply(sign*g,Z[0]);
            U[1] = U[11] = akns_scatter_ES_multiply(sign*g,Z[1]);
            U[4] = U[14] = akns_scatter_ES_multiply(sign*g,Z[2]);
            U[5] = U[15] = f+akns_scatter_ES_multiply(sign*g,Z[3]);
            U[8] = f_d+sign*(g_d*Z[0]+g*Zd[0]);
            U[9] = sign*(g_d*Z[1]+g*Zd[1]);
            U[12] = sign*(g_d*Z[2]+g*Zd[2]);
            U[13] = f_d+sign*(g_d*Z[3]+g*Zd[3]);
        } else {
            U[0] = f+akns_scatter_ES_multiply(sign*g,Z[0]);
            U[1] = akns_scatter_ES_multiply(sign*g,Z[1]);
            U[2] = akns_scatter_ES_multiply(sign*g,Z[2]);
            U[3] = f+akns_scatter_ES_multiply(sign*g,Z[3]);
        }
    }
    return SUCCESS;

leave_fun:
    return ret_code;
}

/**
 * If derivative_flag=0 returns [S11 S12 S21 S22] in result where
 * S = [S11, S12; S21, S22] is the scattering matrix computed using the
 * chosen scheme.
 * If derivative_flag=1 returns [S11 S12 S21 S22 S11' S12' S21' S22'] in
 * result where S11' is the derivative of S11 w.r.t to lambda.
 * Result should be preallocated with size 4*K or 8*K accordingly.
 */
static INT akns_scatter_matrix_impl(UINT const D,
                        COMPLEX const * const q,
                        COMPLEX const * const r,
                        REAL const eps_t,
                        UINT const K,
                        COMPLEX const * const lambda,
                        COMPLEX * const result,
                        COMPLEX * const result_second,
                        INT * const W,
                         akns_discretization_t const discretization,
                         akns_pde_t const PDE,
                         UINT const vanilla_flag,
                         UINT const pade_degree,
                         UINT const derivative_flag,
                         UINT const first_column_flag)
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
    if (first_column_flag && result_second == NULL)
        return E_INVALID_ARGUMENT(result_second);
    if (derivative_flag != 0 && derivative_flag != 1)
        return E_INVALID_ARGUMENT(derivative_flag);
    if (first_column_flag && (derivative_flag
            || PDE != akns_pde_NSE || vanilla_flag
            || (discretization != akns_discretization_ES4
                && discretization != akns_discretization_ES6
                && discretization != akns_discretization_ES8)))
        return E_INVALID_ARGUMENT(first_column_flag);
    if (first_column_flag) {
        for (UINT i=0; i<K; i++) {
            if (CIMAG(lambda[i]) != 0.0)
                return E_INVALID_ARGUMENT(lambda);
        }
    }
    if (pade_degree > 7 || (pade_degree != 0
            && discretization != akns_discretization_ES4
            && discretization != akns_discretization_ES6
            && discretization != akns_discretization_ES8))
        return E_INVALID_ARGUMENT(pade_degree);
    UINT const upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    if (D%upsampling_factor != 0)
        return E_ASSERTION_FAILED;
    const UINT es_degree = upsampling_factor >= 3 ? upsampling_factor-2 : 0;
    const UINT es_numel = 4*(es_degree+1);
    REAL pade_F[8], pade_G[7], pade_denominator[8];
    akns_scatter_pade_coefficients_t pade_coefficients;
    akns_scatter_pade_coefficients_t const *pade = NULL;

    if (pade_degree != 0) {
        ret_code = akns_pade_coefficients(pade_degree,pade_F,pade_G,
                pade_denominator);
        if (ret_code != SUCCESS)
            return ret_code;
        pade_coefficients.degree = pade_degree;
        pade_coefficients.F = pade_F;
        pade_coefficients.G = pade_G;
        pade_coefficients.denominator = pade_denominator;
        pade = &pade_coefficients;
    }

    // Declare pointers that may or may not be used, depending on the discretization.
    // We must do so before possibly jumping to leave_fun.
    COMPLEX *tmp1 = NULL, *tmp2 = NULL, *eps_t_scaled = NULL;

    // Define stepsize constants that are often needed
    REAL const eps_t_2 = eps_t * eps_t;
    REAL const eps_t_3 = eps_t_2 * eps_t;

    // Pre-computing weights required for higher-order CF methods that are
    // independent of q, r and l.
    // For ES4/ES6/ES8, precompute the coefficients of Z(eps_t*l);
    // for TES4, precompute the outer matrix exponentials.

    switch (discretization) {
        case akns_discretization_CT4:
            break;
        case akns_discretization_ES4:
        case akns_discretization_ES6:
        case akns_discretization_ES8:
            if (D/upsampling_factor > (UINT)-1/(es_numel*sizeof(COMPLEX))) {
                ret_code = E_INVALID_ARGUMENT(D);
                goto leave_fun;
            }
            tmp1 = malloc(es_numel*(D/upsampling_factor)*sizeof(COMPLEX));
            CHECK_NOMEM(tmp1,ret_code,leave_fun);
            for (UINT n=0, node=0; n<D; n+=upsampling_factor, node++)
                fnft__akns_es_z_coefficients(upsampling_factor+1,
                        &q[n],&r[n],&tmp1[es_numel*node]);
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
            // fall through
        case akns_discretization_CF4_2:         // commutator-free fourth-order
            // fall through
        case akns_discretization_BO:            // bofetta-osborne scheme
        {
            COMPLEX *qr_weights = NULL;
            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_scaled,discretization);
            CHECK_RETCODE(ret_code, leave_fun_no_eps_t_scaled); // if ret_code != SUCCESS, akns_discretization_method_weights frees qr_weights and eps_t_scaled if needed
            free(qr_weights);
            for (UINT n=0; n<upsampling_factor; n++ )
                eps_t_scaled[n] *= eps_t;
            break;
        }
            
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
                case akns_discretization_ES4:
                case akns_discretization_ES6:
                case akns_discretization_ES8:
                    for (UINT node=0; node<D/upsampling_factor; ) {
                        COMPLEX steps[AKNS_SCATTER_ES_BLOCK_SIZE][16];
                        UINT count = D/upsampling_factor-node;
                        if (count > AKNS_SCATTER_ES_BLOCK_SIZE)
                            count = AKNS_SCATTER_ES_BLOCK_SIZE;
                        ret_code = akns_scatter_U_ES(&tmp1[es_numel*node],
                                es_degree,l_curr,eps_t,pade,1,0,steps[0],count);
                        CHECK_RETCODE(ret_code, leave_fun);
                        for (UINT step=0; step<count; step++) {
                            misc_matrix_mult(4,4,4,steps[step],&H[current][0][0],&H[!current][0][0]);
                            current = !current;
                            if (W != NULL)
                                Wi += akns_scatter_normalize_ES_vector(16,
                                        &H[current][0][0]);
                        }
                        node += count;
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
    } else if (first_column_flag) {
        for (UINT i=0; i<K; i++) {
            const COMPLEX l_curr = lambda[i];
            COMPLEX H[2][2] = {{1.0,0.0}};
            UINT current = 0;
            INT Wi = 0;

            for (UINT node=0; node<D/upsampling_factor; ) {
                COMPLEX steps[AKNS_SCATTER_ES_BLOCK_SIZE][4];
                UINT count = D/upsampling_factor-node;
                if (count > AKNS_SCATTER_ES_BLOCK_SIZE)
                    count = AKNS_SCATTER_ES_BLOCK_SIZE;
                ret_code = akns_scatter_U_ES(&tmp1[es_numel*node],
                        es_degree,l_curr,eps_t,pade,0,0,steps[0],count);
                CHECK_RETCODE(ret_code, leave_fun);
                for (UINT step=0; step<count; step++) {
                    akns_scatter_apply_U_ES(steps[step],2,H[current],1,
                            H[!current],1);
                    current = !current;
                    if (W != NULL)
                        Wi += akns_scatter_normalize_ES_vector(2,H[current]);
                }
                node += count;
            }
            result[i] = H[current][0];
            result_second[i] = H[current][1];
            if (W != NULL)
                W[i] = Wi;
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
                case akns_discretization_ES4:
                case akns_discretization_ES6:
                case akns_discretization_ES8:
                    for (UINT node=0; node<D/upsampling_factor; ) {
                        COMPLEX steps[AKNS_SCATTER_ES_BLOCK_SIZE][4];
                        UINT count = D/upsampling_factor-node;
                        if (count > AKNS_SCATTER_ES_BLOCK_SIZE)
                            count = AKNS_SCATTER_ES_BLOCK_SIZE;
                        ret_code = akns_scatter_U_ES(&tmp1[es_numel*node],
                                es_degree,l_curr,eps_t,pade,0,0,steps[0],count);
                        CHECK_RETCODE(ret_code, leave_fun);
                        for (UINT step=0; step<count; step++) {
                            akns_scatter_apply_U_ES(steps[step],2,
                                    &H[current][0][0],2,&H[!current][0][0],2);
                            akns_scatter_apply_U_ES(steps[step],2,
                                    &H[current][0][1],2,&H[!current][0][1],2);
                            current = !current;
                            if (W != NULL)
                                Wi += akns_scatter_normalize_ES_vector(4,
                                        &H[current][0][0]);
                        }
                        node += count;
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
    return akns_scatter_matrix_impl(D,q,r,eps_t,K,lambda,result,NULL,W,
            discretization,PDE,vanilla_flag,0,derivative_flag,0);
}

INT akns_scatter_matrix_pade(UINT const D,
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
                        UINT const pade_degree,
                        UINT const derivative_flag)
{
    if (pade_degree == 0)
        return E_INVALID_ARGUMENT(pade_degree);
    return akns_scatter_matrix_impl(D,q,r,eps_t,K,lambda,result,NULL,W,
            discretization,PDE,vanilla_flag,pade_degree,
            derivative_flag,0);
}

INT akns_scatter_matrix_first_column(UINT const D,
                        COMPLEX const * const q,
                        COMPLEX const * const r,
                        REAL const eps_t,
                        UINT const K,
                        COMPLEX const * const lambda,
                        COMPLEX * const H11,
                        COMPLEX * const H21,
                        INT * const W,
                        akns_discretization_t const discretization,
                        UINT const pade_degree)
{
    return akns_scatter_matrix_impl(D,q,r,eps_t,K,lambda,H11,H21,W,
            discretization,akns_pde_NSE,0,pade_degree,0,1);
}

/**
 * Returns the a, a_prime and b computed using the chosen scheme.
 */
static INT akns_scatter_bound_states_impl(UINT const D,
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
                               UINT const pade_degree,
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
    if (pade_degree > 7 || (pade_degree != 0
            && discretization != akns_discretization_ES4
            && discretization != akns_discretization_ES6
            && discretization != akns_discretization_ES8))
        return E_INVALID_ARGUMENT(pade_degree);
    UINT const upsampling_factor = akns_discretization_upsampling_factor(discretization);
    if (upsampling_factor == 0)
        return E_INVALID_ARGUMENT(discretization);
    REAL const boundary_coeff = akns_discretization_boundary_coeff(discretization);
    if (boundary_coeff == NAN)
        return E_INVALID_ARGUMENT(>discretization);
    if (D%upsampling_factor != 0)
        return E_ASSERTION_FAILED;
    const UINT es_degree = upsampling_factor >= 3 ? upsampling_factor-2 : 0;
    const UINT es_numel = 4*(es_degree+1);
    UINT const D_given = D/upsampling_factor;
    if (PDE!=akns_pde_KdV && PDE!=akns_pde_NSE)
        return E_INVALID_ARGUMENT(PDE);
    REAL pade_F[8], pade_G[7], pade_denominator[8];
    akns_scatter_pade_coefficients_t pade_coefficients;
    akns_scatter_pade_coefficients_t const *pade = NULL;

    if (pade_degree != 0) {
        ret_code = akns_pade_coefficients(pade_degree,pade_F,pade_G,
                pade_denominator);
        if (ret_code != SUCCESS)
            return ret_code;
        pade_coefficients.degree = pade_degree;
        pade_coefficients.F = pade_F;
        pade_coefficients.G = pade_G;
        pade_coefficients.denominator = pade_denominator;
        pade = &pade_coefficients;
    }

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
    // For ES4/ES6/ES8, precompute the coefficients of Z(eps_t*l);
    // for TES4, precompute the outer matrix exponentials.

    switch (discretization) {
        case akns_discretization_CT4:
            break;
        case akns_discretization_ES4:
        case akns_discretization_ES6:
        case akns_discretization_ES8:
            if (D/upsampling_factor > (UINT)-1/(es_numel*sizeof(COMPLEX))) {
                ret_code = E_INVALID_ARGUMENT(D);
                goto leave_fun;
            }
            tmp1 = malloc(es_numel*(D/upsampling_factor)*sizeof(COMPLEX));
            CHECK_NOMEM(tmp1,ret_code,leave_fun);
            for (UINT n=0, node=0; n<D; n+=upsampling_factor, node++)
                fnft__akns_es_z_coefficients(upsampling_factor+1,
                        &q[n],&r[n],&tmp1[es_numel*node]);
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
            // fall through
        case akns_discretization_CF4_2:         // commutator-free fourth-order
            // fall through
        case akns_discretization_BO:            // bofetta-osborne scheme
        {
            COMPLEX *qr_weights = NULL;
            ret_code = akns_discretization_method_weights(&qr_weights,&eps_t_scaled,discretization);
            CHECK_RETCODE(ret_code, leave_fun_no_eps_t_scaled); // if ret_code != SUCCESS, akns_discretization_method_weights frees qr_weights and eps_t_scaled if needed
            free(qr_weights);
            for (UINT n=0; n<upsampling_factor; n++ )
                eps_t_scaled[n] *= eps_t;
            break;
        }

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
            case akns_discretization_ES4:
            case akns_discretization_ES6:
            case akns_discretization_ES8:
                for (UINT node=0; node<D_given; ) {
                    COMPLEX steps[AKNS_SCATTER_ES_BLOCK_SIZE][16];
                    UINT count = D_given-node;
                    if (count > AKNS_SCATTER_ES_BLOCK_SIZE)
                        count = AKNS_SCATTER_ES_BLOCK_SIZE;
                    ret_code = akns_scatter_U_ES(&tmp1[es_numel*node],
                            es_degree,l_curr,eps_t,pade,1,0,steps[0],count);
                    CHECK_RETCODE(ret_code, leave_fun);
                    for (UINT step=0; step<count; step++) {
                        const UINT n_given = node+step;
                        akns_scatter_apply_U_ES_derivative(steps[step],
                                &PHI[4*n_given],&PHI[4*(n_given+1)]);
                        if (normalization_flag) {
                            WPHI_acc += akns_scatter_normalize_ES_vector(4,
                                    &PHI[4*(n_given+1)]);
                            if (WPHI != NULL)
                                WPHI[n_given+1] = WPHI_acc;
                        }
                    }
                    node += count;
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
                case akns_discretization_ES4:
                case akns_discretization_ES6:
                case akns_discretization_ES8:
                    for (UINT end=D_given; end>0; ) {
                        COMPLEX steps[AKNS_SCATTER_ES_BLOCK_SIZE][4];
                        const UINT count = end < AKNS_SCATTER_ES_BLOCK_SIZE
                                ? end : AKNS_SCATTER_ES_BLOCK_SIZE;
                        const UINT node = end-count;
                        ret_code = akns_scatter_U_ES(&tmp1[es_numel*node],
                                es_degree,l_curr,eps_t,pade,0,1,steps[0],count);
                        CHECK_RETCODE(ret_code, leave_fun);
                        for (UINT step=count; step-->0; ) {
                            const UINT n_given = node+step;
                            akns_scatter_apply_U_ES(steps[step],2,
                                    &PSI[4*(n_given+1)],1,&PSI[4*n_given],1);
                            if (normalization_flag) {
                                WPSI_acc += akns_scatter_normalize_ES_vector(2,
                                        &PSI[4*n_given]);
                                WPSI[n_given] = WPSI_acc;
                            }
                        }
                        end = node;
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
    return akns_scatter_bound_states_impl(D,q,r,T,K,bound_states,a_vals,
            aprime_vals,b_vals,Ws,discretization,PDE,vanilla_flag,0,
            skip_b_flag);
}

INT akns_scatter_bound_states_pade(UINT const D,
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
                              UINT const pade_degree,
                              UINT const skip_b_flag)
{
    if (pade_degree == 0)
        return E_INVALID_ARGUMENT(pade_degree);
    return akns_scatter_bound_states_impl(D,q,r,T,K,bound_states,a_vals,
            aprime_vals,b_vals,Ws,discretization,PDE,vanilla_flag,pade_degree,
            skip_b_flag);
}
