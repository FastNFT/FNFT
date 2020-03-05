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
*/
#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include <string.h>

#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"

// These macros allow us to use the same indices as in the Fortran
// sources, where the first element on ar array is at index 1 and not 0
// and matrices are stored column-wise and not row-wise.
#define VEC_EL(A, n) (A[(n)-1]) // vector element
#define MAT_EL(A, M, m, n) (A[((M)*((n)-1)) + (m)-1]) // matrix element,
// M = number of rows, m = row, n = column
#define DO(i,from,to) for((i)=(from); (i)<=(to); (i)++)
#define RANDOM_NUMBER() ((rand() % (RAND_MAX-1))/RAND_MAX) // random number in
// the half-open interval [0,1) -- we use the modulus since rand()/(RAND_MAX+1)
// might result in integer overflow; the probability of zero is slightly
// increased by this, but let's hope that's ok...

// Interface to the EISCOR routies
extern int z_poly_roots_modified_(INT *N, double complex const * const coeffs,
    double complex * const roots, double *threshold, INT *info);
extern int z_compmat_compress_(int *N_ptr, int *P, double complex *coeffs,
                               double * Q, double * D, double *C, double * B);
extern int l_upr1fact_hess_(int *, int *);
extern int z_upr1fact_qr_(int *vec, int *id, int (*fun)(int *, int *),
                          int *N_ptr, int *P, double *Q, double *D, double *C,
                          double *B, int *M, double complex *V, int *ITS,
                          int *info);
extern int z_upr1fact_deflationcheck_(int *vec, int *N_ptr, int *P, double *Q,
                                      double *D, double *C, double *B,
                                      int *M_ptr, double complex *V,
                                      int *zero_ptr);
extern int z_upr1fact_singleshift_(int *P, double *Q, double *D, double *C,
                                   double *B, double complex *shift);
extern int z_upr1fact_buildbulge_(int *P, double *Q, double *D, double *C,
                                  double *B, double complex *shift, double *G);
extern int z_upr1utri_unimodscale_(int *row, double *D, double *C, double *B,
                                   double complex *scl);
extern int z_upr1fact_singlestep_(int *vec, int (*fun)(int *, int *),
                                  int *N_ptr, int *P, double *Q, double *D,
                                  double *C, double *B, int *M_ptr,
                                  double complex *V, int *ITCNT);
extern int z_upr1utri_decompress_(int *diag, int *N_ptr, double *D, double *C,
                                 double *B, double complex *T);
extern int z_upr1utri_rot3swap_(int *dir, double *D, double *C, double *B,
                                double *G);
extern int z_upr1fact_startchase_(int *vec, int *N_ptr, int *P, double *Q,
                                  double *D, double *C, double *B, int *M_ptr,
                                  double complex *V, int *ITCNT, double *G);
extern int z_upr1fact_chasedown_(int *vec, int *P, double *Q, double *D,
                                 double *C, double *B, int *M_ptr,
                                 double complex *V, double *misfit);
extern int z_upr1fact_endchase_(int *vec, int *N_ptr, int *P, double *Q,
                                double *D, double *C, double *B, int *M_ptr,
                                double complex *V, double *G, int *flag);
extern int d_rot2_vec2gen_(double *A, double *B, double *C, double *S,
                           double *nrm);
extern int z_rot3_vec3gen_(double *AR, double *AI, double *B, double *CR,
                           double *CI, double *S, double *nrm);
extern int z_rot3_vec4gen_(double *AR, double *AI, double *BR, double *BI,
                           double *CR, double *CI, double *S, double *nrm);
extern int z_rot3_fusion_(int *flag, double *G1, double *G2);
extern int u_fixedseed_initialize_(int *info);

void d_rot2_vec2gen(double *A_ptr, double *B_ptr, double *C_ptr, double *S_ptr,
                    double *nrm_ptr)
{
    const double A = *A_ptr;
    const double B = *B_ptr;
    if ((A == 0) && (B == 0)) {
        *C_ptr = 1;
        *S_ptr = 0;
        *nrm_ptr = 0;
    } else if (FABS(A) >= FABS(B)) {
        const double S = B / A;
        const double nrm = copysign(SQRT(1 + S*S), A);
        const double C = 1 / nrm;
        *C_ptr = C;
        *S_ptr = S*C;
        *nrm_ptr = A*nrm;
    } else {
        const double C = A / B;
        const double nrm = copysign(SQRT(1 + C*C), B);
        const double S = 1 / nrm;
        *S_ptr = S;
        *C_ptr = C*S;
        *nrm_ptr = B*nrm;
    }
}

void z_rot3_vec3gen(double *AR, double *AI, double *B, double *CR,
                    double *CI, double *S, double *nrm)
{
    double tar, tai, tb;

    const double nar = FABS(*AR);
    const double nai = FABS(*AI);
    const double nb = FABS(*B);

    if (nar == 0 && nai == 0 && nb == 0) {
        *CR = 0;
        *CI = 0;
        *S = 0;
        *nrm = 0;
    } else if (nar >= nb && nar >= nai) {
        tb = *B / *AR;
        tai = *AI / *AR;
        *nrm = copysign(SQRT(1 + tb*tb + tai*tai), *AR);
        *CR = 1 / *nrm;
        *CI = tai*(*CR);
        *S = tb*(*CR);
        *nrm = (*AR)*(*nrm);
    } else if (nai >= nb && nai >= nar) {
        tb = *B / *AI;
        tar = *AR / *AI;
        *nrm = copysign(SQRT(1 + tb*tb + tar*tar), *AI);
        *CI = 1 / *nrm;
        *CR = tar*(*CI);
        *S = tb*(*CI);
        *nrm = (*AI)*(*nrm);
    } else {
        tar = *AR / *B;
        tai = *AI / *B;
        *nrm = copysign(SQRT(1 + tai*tai + tar*tar), *B);
        *S = 1 / *nrm;
        *CR = tar*(*S);
        *CI = tai*(*S);
        *nrm = (*B)*(*nrm);
    }
}

void z_rot3_vec4gen(double *AR, double *AI, double *BR, double *BI,
                    double *CR, double *CI, double *S, double *nrm)
{
    double pAr, pAi, pBr, pBi;
    d_rot2_vec2gen(BR, BI, &pBr, &pBi, S);
    d_rot2_vec2gen(AR, AI, &pAr, &pAi, CR);
    // TODO: The conditions here look weird. Check.
    if (FABS(*S)<=INFINITY && FABS(*CR)>INFINITY && pAr==0)
        z_rot3_vec3gen(AR, AI, S, CR, CI, S, nrm);
    else if (FABS(*S)<=INFINITY && FABS(*CR)>INFINITY && pAi==0)
        z_rot3_vec3gen(AR, AI, S, CR, CI, S, nrm);
    else {
        double temp1 = pAr*pBr + pAi*pBi;
        double temp2 = -pAr*pBi + pAi*pBr;
        d_rot2_vec2gen(&temp1, &temp2, &pAr, &pAi, CI);
        temp1 = (*CR)*pAr;
        temp2 = (*CR)*pAi;
        z_rot3_vec3gen(&temp1, &temp2, S, CR, CI, S, nrm);
    }
}

void z_compmat_compress(int *N_ptr, int *P, double complex *coeffs,
                        double *Q, double *D, double *C, double *B)
{
    int i;
    double phr, phi, nrm, beta;
    double complex temp;
    const int N = *N_ptr;

    DO(i, 1, 3*(N-1))
        VEC_EL(Q, i) = 0;
    DO(i, 1, N-1)
        VEC_EL(Q, 3*i) = 1;

    DO(i, 1, N-2)
        if (VEC_EL(P, N-i-1)) {
            temp = VEC_EL(coeffs, N-i);
            if ((N-i-1)%2 == 1)
                temp = -temp;
            for (int j=N-i; j>=2; j--)
                VEC_EL(coeffs, j) = VEC_EL(coeffs, j-1);
            VEC_EL(coeffs, N-i) = temp;
        }

    DO(i, 1, 2*N)
        VEC_EL(D, i) = 0;
    DO(i, 1, N)
        VEC_EL(D, 2*i-1) = 1;
    DO(i, 1, 3*N) {
        VEC_EL(B, i) = 0;
        VEC_EL(C, i) = 0;
    }

    double *dc = (double*)&VEC_EL(coeffs, N);
    d_rot2_vec2gen(&dc[0], &dc[1], &phr, &phi, &beta);

    VEC_EL(D, 2*N-1) = phr;
    VEC_EL(D, 2*N) = phi;

    double min_one = -1;
    d_rot2_vec2gen(&beta, &min_one, &VEC_EL(C, 3*N-2), &VEC_EL(C, 3*N), &nrm);

    VEC_EL(B, 3*N-2) = VEC_EL(C, 3*N);
    VEC_EL(B, 3*N) = VEC_EL(C, 3*N-2);

    temp = nrm;
    DO(i, 1, N-1) {
        dc = (double*)&VEC_EL(coeffs, N-i);
        double *dt = (double*)&temp;
        z_rot3_vec4gen_(&dc[0], &dc[1], &dt[0], &dt[1], &VEC_EL(C, 3*(N-i)-2),
                        &VEC_EL(C, 3*(N-i)-1), &VEC_EL(C, 3*(N-i)), &nrm);

        VEC_EL(B, 3*(N-i)-2) = VEC_EL(C, 3*(N-i)-2);
        VEC_EL(B, 3*(N-i)-1) = -VEC_EL(C, 3*(N-i)-1);
        VEC_EL(B, 3*(N-i)) = -VEC_EL(C, 3*(N-i));

        temp = (VEC_EL(C, 3*(N-i)-2) - I*VEC_EL(C, 3*(N-i)-1))*VEC_EL(coeffs, N-i)
            + VEC_EL(C, 3*(N-i))*temp;
    }
}
void z_upr1utri_unimodscale(int *row, double *D, double *C, double *B,
                            double complex *scl_ptr)
{
    double nrm;
    double complex temp;
    double *temp_re = (double*)&temp;
    double *temp_im = temp_re+1;

    const double complex scl = *scl_ptr;

    temp = scl*(VEC_EL(D, 1) + I*VEC_EL(D, 2));
    d_rot2_vec2gen(temp_re, temp_im, &VEC_EL(D, 1), &VEC_EL(D, 2), &nrm);
    if (! *row) {
        temp = scl*(VEC_EL(B, 1) + I*VEC_EL(B, 2));
        z_rot3_vec3gen(temp_re, temp_im, &VEC_EL(B, 3), &VEC_EL(B, 1),
                       &VEC_EL(B, 2), &VEC_EL(B, 3), &nrm);
        temp = CONJ(scl)*(VEC_EL(C, 1) + I*VEC_EL(C, 2));
        z_rot3_vec3gen(temp_re, temp_im, &VEC_EL(C, 3), &VEC_EL(C, 1),
                       &VEC_EL(C, 2), &VEC_EL(C, 3), &nrm);
    }
}

void z_upr1fact_deflationcheck(int *vec_ptr, int *N_ptr, int *P, double *Q,
                                double *D, double *C, double *B,
                                int *M_ptr, double complex *V,
                                int *zero_ptr)
{
    int i, k;
    const double tol = EPSILON;
    double  qr, qi, nrm;
    const int vec = *vec_ptr;
    const int N = *N_ptr;
    const int M = *M_ptr;
    int zero = *zero_ptr;
    double complex scl;

    int t = 1;
    int f = 0;

    DO(i, 1, N-1) {
        nrm = FABS(VEC_EL(Q, 3*(N-i)));
        if (nrm < tol) {
            zero = N-i;
            if (zero < 0)
                zero = 0;

            qr = VEC_EL(Q, 3*zero-2);
            qi = VEC_EL(Q, 3*zero-1);

            VEC_EL(Q, 3*zero-2) = 1;
            VEC_EL(Q, 3*zero-1) = 0;
            VEC_EL(Q, 3*zero) = 0;

            if (zero == 1) {
                scl = qr + I*qi;
                z_upr1utri_unimodscale(&t, &VEC_EL(D, 2*zero-1),
                                       &VEC_EL(C, 3*zero-2),
                                       &VEC_EL(B, 3*zero-2),
                                       &scl);
            } else if (! VEC_EL(P, zero-1)) {
                scl = qr + I*qi;
                z_upr1utri_unimodscale(&t, &VEC_EL(D, 2*zero-1),
                                       &VEC_EL(C, 3*zero-2),
                                       &VEC_EL(B, 3*zero-2),
                                       &scl);
            } else {
                scl = qr + I*qi;
                z_upr1utri_unimodscale(&f, &VEC_EL(D, 2*zero-1),
                                       &VEC_EL(C, 3*zero-2),
                                       &VEC_EL(B, 3*zero-2),
                                       &scl);
                if (vec) {
                    DO(k, 1, M)
                        MAT_EL(V, M, k, zero) *= qr + I*qi;
                }
            }

            if (zero == N-1) {
                scl = qr - I*qi;
                z_upr1utri_unimodscale(&t, &VEC_EL(D, 2*zero+1),
                                       &VEC_EL(C, 3*zero+1),
                                       &VEC_EL(B, 3*zero+1),
                                       &scl);
            } else if (! VEC_EL(P, zero)) {
                scl = qr - I*qi;
                z_upr1utri_unimodscale(&f, &VEC_EL(D, 2*zero+1),
                                       &VEC_EL(C, 3*zero+1),
                                       &VEC_EL(B, 3*zero+1),
                                       &scl);
                if (vec) {
                    DO(k, 1, M)
                        MAT_EL(V, M, k, zero+1) *= qr - I*qi;
                }
            } else {
                scl = qr - I*qi;
                z_upr1utri_unimodscale(&t, &VEC_EL(D, 2*zero+1),
                                       &VEC_EL(C, 3*zero+1),
                                       &VEC_EL(B, 3*zero+1),
                                       &scl);
            }
        }
    }

    *zero_ptr = zero;
}
void z_upr1fact_startchase(int *vec_ptr, int *N_ptr, int *P, double *Q,
                           double *D, double *C, double *B, int *M_ptr,
                           double complex *V, int *ITCNT_ptr, double *G)
{
    int ir1, ir2, id1, id2;
    int tp[2] = {0};
    double Ginv[3];
    double tq[6] = {0}, td[6] = {0}, tc[9] = {0}, tb[9] = {0};
    double complex shift = 0;
    int i, m;

    const int vec = *vec_ptr;
    const int N = *N_ptr;
    const int M = *M_ptr;
    const int ITCNT = *ITCNT_ptr;

    if ((ITCNT%20 == 0) && (ITCNT > 0)) {

        VEC_EL(G, 1) = 0.5;//RANDOM_NUMBER();
        VEC_EL(G, 2) = 0.25;//RANDOM_NUMBER();
        shift = VEC_EL(G, 1) + I*VEC_EL(G, 2);

    } else {

        if (N < 3) {

            // note: tp, tq, td, tc, and tb were initialized to zero (see above)
            VEC_EL(tq, 1) = 1;
            DO(i, 4, 6)
                VEC_EL(tq, i) = VEC_EL(Q, i-3);

            VEC_EL(td, 1) = 1;
            DO(i, 3, 6)
                VEC_EL(td, i) = VEC_EL(D, i-2);

            VEC_EL(tc, 3) = 1;
            DO(i, 4, 9)
                VEC_EL(tc, i) = VEC_EL(C, i-3);

            VEC_EL(tb, 3) = -1;
            DO(i, 4, 9)
                VEC_EL(tb, i) = VEC_EL(B, i-3);

        } else {

            if (N == 3) {
                VEC_EL(tp, 1) = 0;
                VEC_EL(tp, 2) = VEC_EL(P, N-2);
            } else {
                VEC_EL(tp, 1) = VEC_EL(P, N-3);
                VEC_EL(tp, 2) = VEC_EL(P, N-2);
            }

            ir2 = 3*N; ir1 = ir2 - 8;
            id2 = 2*N; id1 = id2 - 5;
            DO (i, ir1, ir2-3)
                VEC_EL(tq, i-ir1+1) = VEC_EL(Q, i);
            DO (i, id1, id2)
                VEC_EL(td, i-id1+1) = VEC_EL(D, i);
            DO (i, ir1, ir2) {
                VEC_EL(tc, i-ir1+1) = VEC_EL(C, i);
                VEC_EL(tb, i-ir1+1) = VEC_EL(B, i);
            }

        }

        z_upr1fact_singleshift_(tp, tq, td, tc, tb, &shift);
    }

    z_upr1fact_buildbulge_(P, Q, D, C, B, &shift, G);
    VEC_EL(Ginv, 1) = VEC_EL(G, 1);
    VEC_EL(Ginv, 2) = -VEC_EL(G, 2);
    VEC_EL(Ginv, 3) = -VEC_EL(G, 3);

    if (vec) {
        const double complex A11 = VEC_EL(G, 1) + I*VEC_EL(G, 2);
        const double complex A21 = VEC_EL(G, 3);
        const double complex A12 = -A21;
        const double complex A22 = conj(A11);

        DO (m, 1, M) {
            const double complex Vm1 = MAT_EL(V, M, m, 1)*A11
                + MAT_EL(V, M, m, 2)*A21;
            const double complex Vm2 = MAT_EL(V, M, m, 1)*A12
                + MAT_EL(V, M, m, 2)*A22;
            MAT_EL(V, M, m, 1) = Vm1;
            MAT_EL(V, M, m, 2) = Vm2;
        }
    }

    int f = 0;
    int t = 1;

    if (! VEC_EL(P, 1)) {

        z_rot3_fusion_(&f, Ginv, Q);
        z_upr1utri_rot3swap_(&f, D, C, B, G);
        double complex tmp = VEC_EL(Ginv, 1) + I*VEC_EL(Ginv, 2);
        if (vec) {
            DO (m, 1, M) {
                MAT_EL(V, M, m, 1) *= tmp;
                MAT_EL(V, M, m, 2) *= conj(tmp);
            }
        }
        z_upr1utri_unimodscale(&f, D, C, B, &tmp);
        tmp = VEC_EL(Ginv, 1) - I*VEC_EL(Ginv, 2);
        z_upr1utri_unimodscale(&f, &VEC_EL(D, 3), &VEC_EL(C, 4), &VEC_EL(B, 4),
                               &tmp);

    } else {

        z_upr1utri_rot3swap_(&f, D, C, B, G);
        z_rot3_fusion_(&t, Q, G);
        double complex tmp = VEC_EL(G, 1) + I*VEC_EL(G, 2);
        z_upr1utri_unimodscale(&t, D, C, B, &tmp);
        tmp = VEC_EL(G, 1) - I*VEC_EL(G, 2);
        z_upr1utri_unimodscale(&t, &VEC_EL(D, 3), &VEC_EL(C, 4), &VEC_EL(B, 4),
                               &tmp);
        DO (i, 1, 3)
            VEC_EL(G, i) = VEC_EL(Ginv, i);

    }
}

void z_upr1fact_singlestep(int *vec, int (*fun)(int *, int *),
                           int *N_ptr, int *P, double *Q, double *D,
                           double *C, double *B, int *M_ptr,
                           double complex *V, int *ITCNT)
{
    int i, ir1, id1, final_flag;
    double misfit[3] = {0};
    const int N = *N_ptr;
    const int M = *M_ptr;

    if (N < 3)
        final_flag = 0;
    else
        final_flag = (*fun)(N_ptr, P);

    z_upr1fact_startchase_(vec, N_ptr, P, Q, D, C, B, M_ptr,
                          &MAT_EL(V, M, 1, 1), ITCNT, misfit);

    DO(i, 1, N-3) {
        ir1 = 3*i+1;
        // ir2 = 3*(i+2); // not needed
        id1 = 2*i+1;
        // id2 = 2*(i+2); // not needed
        z_upr1fact_chasedown_(vec, &VEC_EL(P, i), &VEC_EL(Q, ir1-3),
                              &VEC_EL(D, id1), &VEC_EL(C, ir1),
                              &VEC_EL(B, ir1), M_ptr,
                              &MAT_EL(V, M, 1, i+1), misfit);
    }
    z_upr1fact_endchase_(vec, N_ptr, P, Q, D, C, B, M_ptr, V, misfit,
                         &final_flag);
}

void z_upr1fact_qr(int *vec, int *id, int (*fun)(int *, int *),
                   int *N_ptr, int *P, double *Q, double *D, double *C,
                   double *B, int *M_ptr, double complex *V, int *ITS,
                   int *info)
{
    int i, k;
    int STR, STP, ZERO, ITMAX, ITCNT;
    const int N = *N_ptr;
    const int M = *M_ptr;
    int temp;

    *info = 0;

    // NOTE: The original code performs several checks
    // here that have not been ported.

    DO(i, 1, N-1)
        VEC_EL(ITS, i) = 0;

    if (*vec && *id) {
        DO(i, 1, M)
            DO(k, 1, N)
                MAT_EL(V, M, i, k) = 0;
        DO(i, 1, N)
            MAT_EL(V, M, i, i) = 1;
    }

    STR = 1;
    STP = N-1;
    ZERO = 0;
    ITMAX = 20*N;
    ITCNT = 0;

    DO(k, 1, ITMAX) {
        if (STP <= 0)
            return;
        temp = STP-STR+2;
        z_upr1fact_deflationcheck(vec, &temp, &VEC_EL(P, STR),
                                  &VEC_EL(Q, 3*STR-2), &VEC_EL(D, 2*STR-1),
                                  &VEC_EL(C, 3*STR-2), &VEC_EL(B, 3*STR-2),
                                  M_ptr, &MAT_EL(V, M, 1, STR), &ZERO);
        if (STP == STR+ZERO-1) {
            VEC_EL(ITS, STR+ZERO-1) = ITCNT;
            ITCNT = 0;
            STP--;
            ZERO = 0;
            STR = 1;
        } else {
            if (ZERO > 0) {
                STR += ZERO;
                ZERO = 0;
                VEC_EL(ITS, STR+ZERO-1) = ITCNT;
                ITCNT = 0;
            }
            temp = STP-STR+2;
            z_upr1fact_singlestep(vec, fun, &temp, &VEC_EL(P, STR),
                                  &VEC_EL(Q, 3*STR-2), &VEC_EL(D, 2*STR-1),
                                  &VEC_EL(C, 3*STR-2), &VEC_EL(B, 3*STR-2),
                                  M_ptr, &MAT_EL(V, M, 1, STR), &ITCNT);
            if (ITCNT == -1)
                ITCNT = 1;
            else
                ITCNT++;
        }

        if (k == ITMAX) {
            *info = 1;
            VEC_EL(ITS, STR+STP-1) = ITCNT;
        }
    }
}

int l_upr1fact_hess(int *N, int *P)
{
    (void) N; // suppress unused variable compiler warning
    (void) P;
    return 0;
}

INT z_poly_roots_modified(int *N_ptr, double complex const * const coeffs,
                          double complex * const roots, int *info)
{
    int i;
    int *P = NULL;
    int *ITS = NULL;
    double *Q = NULL, *D1 = NULL, *C1 = NULL, *B1 = NULL;
    double complex *V = NULL, *W = NULL;
    double complex sclc;
    INT ret_code = SUCCESS;

    int f = 0; // false in Fortran
    int t = 1; // true in Fortran

    const int N = *N_ptr;
    P = malloc( (N-2) * sizeof(int) );
    ITS = malloc( (N-1) * sizeof(int) );
    Q = malloc( 3*(N-1) * sizeof(double) );
    D1 = malloc( 2*(N+1) * sizeof(double) );
    C1 = malloc( 3*N * sizeof(double) );
    B1 = malloc( 3*N * sizeof(double) );
    V = malloc( N*N * sizeof(double complex) );
    W = malloc( N * sizeof(double complex) );
    if (P == NULL || ITS == NULL || Q == NULL || D1 == NULL || C1 == NULL
        || B1 == NULL || V == NULL || W == NULL ) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    DO(i, 1, N-2)
          VEC_EL(P, i) = 0;

    sclc = VEC_EL(coeffs, 1);
    VEC_EL(V, N) = VEC_EL(coeffs, N+1)/sclc;
    if (N%2 == 1)
        VEC_EL(V, N) = -VEC_EL(V, N);
    DO(i, 1, N-1)
        VEC_EL(V, i) = -VEC_EL(coeffs, N+1-i)/sclc;

    z_compmat_compress(N_ptr, P, V, Q, D1, C1, B1);
    z_upr1fact_qr(&f, &f, &l_upr1fact_hess, N_ptr, P, Q, D1, C1, B1, N_ptr, V,
                  ITS, info);
    if (*info != 0) {
        ret_code = E_SUBROUTINE(*info);
        *info = 1;
        goto leave_fun;
    }
    z_upr1utri_decompress_(&t, N_ptr, D1, C1, B1, roots);

leave_fun:
    free(P);
    free(ITS);
    free(Q);
    free(D1);
    free(C1);
    free(B1);
    free(V);
    free(W);

    return ret_code;
}


// Fast computation of polynomial roots. See the header file for details.
INT poly_roots_fasteigen(const UINT deg,
    COMPLEX const * const p, COMPLEX * const roots)
{
    INT int_deg, info;

    // Check inputs
    if (p == NULL)
        return E_INVALID_ARGUMENT(p);
    if (roots == NULL)
        return E_INVALID_ARGUMENT(roots);

    u_fixedseed_initialize_(&info);
    if (info != 0)
        return E_SUBROUTINE(FNFT_EC_OTHER);

    int_deg = (int)deg;
    z_poly_roots_modified(&int_deg, p, roots, &info);

    if (info == 0)
        return SUCCESS;
    else
        return E_SUBROUTINE(FNFT_EC_OTHER);
}

// Fast computation of polynomial roots. See the header file for details.
INT poly_roots_fasteigen_(const UINT deg,
                          COMPLEX const * const p, COMPLEX * const roots)
{
    INT int_deg, info;
    double threshold = INFINITY;
    // This threshold was used in the original routine. Set to INFINITY to
    // enforce QR. Set to 0 to enforce QZ.

    // Check inputs
    if (p == NULL)
        return E_INVALID_ARGUMENT(p);
    if (roots == NULL)
        return E_INVALID_ARGUMENT(roots);

    u_fixedseed_initialize_(&info);
    if (info != 0)
        return E_SUBROUTINE(FNFT_EC_OTHER);

    // Call Fortran root finding routine
    int_deg = (int)deg;
    z_poly_roots_modified_(&int_deg, p, roots, &threshold, &info);

    if (info == 0)
        return SUCCESS;
    else
        return E_SUBROUTINE(FNFT_EC_OTHER);
}
