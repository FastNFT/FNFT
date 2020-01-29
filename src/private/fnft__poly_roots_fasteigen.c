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

#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"

// Interface to the EISCOR root finding routine
extern int z_poly_roots_modified_(INT *N, double complex const * const coeffs,
    double complex * const roots, double *threshold, INT *info);

extern int z_compmat_compress_(int *N_ptr, int *P, double complex *coeffs,
                               double * Q, double * D, double *C, double * B);

extern int l_upr1fact_hess_(int *, int *);

extern int z_upr1fact_qr_(int *vec, int *id, int (*fun)(int *, int *),
                          int *N_ptr, int *P, double *Q, double *D, double *C,
                          double *B, int *M, double complex *V, int *ITS,
                          int *info);

extern int z_upr1utri_decompress_(int *diag, int *N_ptr, double *D, double *C,
                                 double *B, double complex *T);

INT z_poly_roots_modified(int *N_ptr, double complex const * const coeffs,
                          double complex * const roots, int *info)
{
    int ii;
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
    V = malloc( N * sizeof(double complex) );
    W = malloc( N * sizeof(double complex) );
    if (P == NULL || ITS == NULL || Q == NULL || D1 == NULL || C1 == NULL
        || B1 == NULL || V == NULL || W == NULL ) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    for (ii=0; ii<N-2; ii++)
        P[ii] = 0;

    sclc = coeffs[0];
    V[N-1] = coeffs[N]/sclc;
    if (N%2 == 1)
        V[N-1] = -V[N-1];
    for (ii=0; ii<N-1; ii++)
        V[ii] = -coeffs[N-1-ii]/sclc;

    z_compmat_compress_(N_ptr, P, V, Q, D1, C1, B1);
    z_upr1fact_qr_(&f, &f, &l_upr1fact_hess_, N_ptr, P, Q, D1, C1, B1, N_ptr, V,
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
    double threshold = 1e8;
    // This threshold was used in the original routine. Set to INFINITY to
    // enforce QR. Set to 0 to enforce QZ.

  // Check inputs
  if (p == NULL)
    return E_INVALID_ARGUMENT(p);
  if (roots == NULL)
    return E_INVALID_ARGUMENT(roots);

    // Call Fortran root finding routine
    int_deg = (int)deg;
    //z_poly_roots_modified_(&int_deg, p, roots, &threshold, &info);
    z_poly_roots_modified(&int_deg, p, roots, &info);

    if (info == 0)
        return SUCCESS;
    else
        return E_SUBROUTINE(FNFT_EC_OTHER);
}
