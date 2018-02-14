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
* Shrinivas Chimmalgi (TU Delft) 2017.
* Peter J Prins (TU Delft) 2017.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"
#include "fnft__nse_fscatter.h"
#include "fnft__nse_scatter.h"
#include "fnft__nse_discretization.h"

/**
 * Returns the length (in number of elements) for "result" in nse_fscatter
 * or 0 if either the discretization is unknown or D=0.
 */
UINT nse_fscatter_numel(UINT D,
        nse_discretization_t discretization)
{
    // 2x2 matrix of degree+1 elements
    return 4*D*(nse_discretization_degree(discretization) + 1);
}

/**
 * result needs to be pre-allocated with size 4*(deg+1)*D*sizeof(COMPLEX)
 */
INT nse_fscatter(const UINT D, COMPLEX const * const q,
    const REAL eps_t, const INT kappa,
    COMPLEX * const result, UINT * const deg_ptr,
    INT * const W_ptr, nse_discretization_t discretization)
{
    INT i, ret_code;
    UINT len;
    COMPLEX *p, *p11, *p12, *p21, *p22;
    REAL scl;
    REAL Q_abs;
    COMPLEX q_arg;
    COMPLEX qt, rt, B11, B12, B21, B22;
    
    // Check inputs
    if (D == 0)
        return E_INVALID_ARGUMENT(D);
    if (q == NULL)
        return E_INVALID_ARGUMENT(q);
    if (eps_t <= 0.0)
        return E_INVALID_ARGUMENT(eps_t);
    if (abs(kappa) != 1)
        return E_INVALID_ARGUMENT(kappa);
    if (result == NULL)
        return E_INVALID_ARGUMENT(result);
    if (deg_ptr == NULL)
        return E_INVALID_ARGUMENT(deg_ptr);
    
    // Allocate buffers
    len = nse_fscatter_numel(D, discretization);
    if (len == 0) { // size D>0, this means unknown discretization
        return E_INVALID_ARGUMENT(opts->discretization);
    }
    p = malloc(len*sizeof(COMPLEX));
    // degree 1 polynomials
    if (p == NULL)
        return E_NOMEM;
    
    switch (discretization) {
        
        case nse_discretization_2SPLIT2_MODAL: // Modified Ablowitz-Ladik discretization
            // Set the individual scattering matrices up
            p11 = p;
            p12 = p11 + 2*D;
            p21 = p12 + 2*D;
            p22 = p21 + 2*D;
            for (i=D-1; i>=0; i--) {
                // compute the scaling factor scl = 1/sqrt(1+kappa*|eps_t*q[i]|^2)
                scl = eps_t*CABS(q[i]);
                if (kappa == -1) {
                    if (scl >= 1.0) {
                        ret_code = E_OTHER("kappa == -1 but eps_t*|q[i]|>=1 ... decrease step size");
                        goto release_mem;
                    }
                    scl = 1.0 / ( SQRT(1.0 + scl) * SQRT(1.0 - scl) );
                } else
                    scl = 1.0 / HYPOT(1.0, scl);
                
                // construct the scattering matrix for the i-th sample
                p11[0] = 0.0;
                p11[1] = scl;
                p12[0] = scl*eps_t*q[i];
                p12[1] = 0.0;
                p21[0] = 0.0;
                p21[1] = -kappa*scl*eps_t*CONJ(q[i]);
                p22[0] = scl;
                p22[1] = 0.0;
                p11 += 2;
                p12 += 2;
                p21 += 2;
                p22 += 2;
            }
            
            break;
            
        case nse_discretization_2SPLIT4B: // Fifth order splitting dual scheme add citation
            
            
            // Set the individual scattering matrices up
            p11 = p;
            p12 = p11 + 3*D;
            p21 = p12 + 3*D;
            p22 = p21 + 3*D;
            
            for (i = D-1; i >= 0; i--) {
                // compute few values
                qt = q[i];
                rt = -kappa*CONJ(q[i]);
                if (qt == 0){
                    B11 = 1; B12 = 0; B21 = 0; B22 = 1;
                }else{
                    B11 = CCOSH(0.25*eps_t*CSQRT(qt)*CSQRT(rt));
                    B12 = (CSINH(0.25*eps_t*CSQRT(rt*qt))*CSQRT(rt*qt))/rt;
                    B21 = -kappa * CONJ(B12);
                    B22 = B11;
                }
                // construct the scattering matrix for the i-th sample
                p11[0] = B11*B11*B11*B11 + (2*B12*B21*B11*B11)/3 - (B12*B12*B21*B21)/3;
                p11[1] = (8*B12*B21*B11*B11)/3 + (8*B12*B21*B22*B11)/3;
                p11[2] = (4*B12*B12*B21*B21)/3 + B12*B21*B22*B22 - (2*B11*B12*B21*B22)/3 - (B11*B11*B12*B21)/3;
                p12[0] = -(B12*(- 3*B11*B11*B11 + B22*B11*B11 - 3*B12*B21*B11 + B12*B21*B22))/3;
                p12[1] = (B12*(4*B11*B11*B22 + 4*B11*B22*B22 + 4*B12*B21*B11 + 4*B12*B21*B22))/3;
                p12[2] = -(B12*(- 3*B22*B22*B22 + B11*B22*B22 - 3*B12*B21*B22 + B11*B12*B21))/3;
                p21[0] = -(B21*(- 3*B11*B11*B11 + B22*B11*B11 - 3*B12*B21*B11 + B12*B21*B22))/3;
                p21[1] = (B21*(4*B11*B11*B22 + 4*B11*B22*B22 + 4*B12*B21*B11 + 4*B12*B21*B22))/3;
                p21[2] = -(B21*(- 3*B22*B22*B22 + B11*B22*B22 - 3*B12*B21*B22 + B11*B12*B21))/3;
                p22[0] = B11*B11*B12*B21 - (2*B22*B11*B12*B21)/3 + (4*B12*B12*B21*B21)/3 - (B22*B22*B12*B21)/3;
                p22[1] = (8*B12*B21*B22*B22)/3 + (8*B11*B12*B21*B22)/3;
                p22[2] = B22*B22*B22*B22 + (2*B12*B21*B22*B22)/3 - (B12*B12*B21*B21)/3;
                p11 += 3;
                p12 += 3;
                p21 += 3;
                p22 += 3;
            }
            
            break;
            
        case nse_discretization_2SPLIT4A: // Fifth order splitting scheme add citation
            
            // Set the individual scattering matrices up
            p11 = p;
            p12 = p11 + 5*D;
            p21 = p12 + 5*D;
            p22 = p21 + 5*D;
            for (i = D-1; i >= 0; i--) {
                // compute few values
                qt = q[i];
                rt = -kappa*CONJ(q[i]);
                if (qt == 0){
                    B11 = 1; B12 = 0; B21 = 0; B22 = 1;
                }else{
                    B11 = CCOSH(0.5*eps_t*CSQRT(qt)*CSQRT(rt));
                    B12 = (CSINH(0.5*eps_t*CSQRT(rt*qt))*CSQRT(rt*qt))/rt;
                    B21 = -kappa * CONJ(B12);
                    B22 = B11;
                }
                
                // construct the scattering matrix for the i-th sample
                p11[0] = B11*B11-(B12*B21)/3;
                p11[1] = 0.0;
                p11[2] = (4*B12*B21)/3;
                p11[3] = 0.0;
                p11[4] = 0.0;
                p12[0] = 0.0;
                p12[1] = (4*B11*B12)/3;
                p12[2] = -(B12*B22)/3-(B11*B12)/3;
                p12[3] = (4*B12*B22)/3;
                p12[4] = 0.0;
                p21[0] = 0.0;
                p21[1] = (4*B11*B21)/3;
                p21[2] = -(B21*B22)/3-(B11*B21)/3;
                p21[3] = (4*B21*B22)/3;
                p21[4] = 0.0;
                p22[0] = 0.0;
                p22[1] = 0.0;
                p22[2] = (4*B12*B21)/3;
                p22[3] = 0.0;
                p22[4] = B22*B22-(B12*B21)/3;
                p11 += 5;
                p12 += 5;
                p21 += 5;
                p22 += 5;
            }
            
            break;
            
        case nse_discretization_2SPLIT2A: // Boffetta-Osborne plus Strang splitting
            
            // Set the individual scattering matrices up
            p11 = p;
            p12 = p11 + 2*D;
            p21 = p12 + 2*D;
            p22 = p21 + 2*D;
            
            
            if (kappa == -1) { // defocusing case
                
                for (i=D-1; i>=0; i--) {
                    
                    // pre-compute a few values
                    Q_abs = eps_t * CABS(q[i]);
                    q_arg = CEXP( I * CARG(q[i]) );
                    
                    // construct the scattering matrix for the i-th sample
                    p11[0] = 0.0;
                    p11[1] = COSH(Q_abs);
                    p12[0] = SINH(Q_abs);
                    p21[1] = CONJ(q_arg)*p12[0];
                    p12[0] *= q_arg;
                    p12[1] = 0.0;
                    p21[0] = 0.0;
                    // p21[1] was already set earlier
                    p22[0] = p11[1];
                    p22[1] = 0.0;
                    p11 += 2;
                    p12 += 2;
                    p21 += 2;
                    p22 += 2;
                }
                
            } else { // focusing case
                
                for (i=D-1; i>=0; i--) {
                    
                    // pre-compute a few values
                    Q_abs = eps_t * CABS(q[i]);
                    q_arg = CEXP( I * CARG(q[i]) );
                    
                    // construct the scattering matrix for the i-th sample
                    p11[0] = 0.0;
                    p11[1] = COS(Q_abs); // since cosh(ix)=cos(x)
                    p12[0] = I*SIN(Q_abs); // since sinh(ix)=i*sin(x)
                    p21[1] = I*CONJ(q_arg)*p12[0];
                    p12[0] *= q_arg * (-I);
                    p12[1] = 0.0;
                    p21[0] = 0.0;
                    // p21[1] was already set earlier
                    p22[0] = p11[1];
                    p22[1] = 0.0;
                    p11 += 2;
                    p12 += 2;
                    p21 += 2;
                    p22 += 2;
                }
            }
            
            break;
            
        default: // Unknown discretization
            
            ret_code = E_INVALID_ARGUMENT(discretization);
            goto release_mem;
    }
    
    // Multiply the individual scattering matrices
    *deg_ptr = nse_discretization_degree(discretization);
    if (*deg_ptr == 0) {
        ret_code = E_INVALID_ARGUMENT(discretization);
        goto release_mem;
    }
    ret_code = poly_fmult2x2(deg_ptr, D, p, result, W_ptr);
    if (ret_code != SUCCESS)
        ret_code = E_SUBROUTINE(ret_code);
    
    release_mem:
        free(p);
        return ret_code;
}
