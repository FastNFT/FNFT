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
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>
#include "fnft__nse_fscatter.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"

static INT nse_fscatter_test_defocusing_2split4B()
{
    UINT i, D = 8, deg;
    INT W;
    REAL scl;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8];
    COMPLEX result[3*4*8];
    COMPLEX result_exact[68] = {
           1.00578402514809 -  5.83210136741117e-05*I,
         0.0155954521142025 -  0.000467093980630542*I,
         0.0124758235986958 -    0.0012200695988487*I,
        0.00854978166679528 -   0.00180605696970404*I,
        0.00794998805226457 -   0.00330086725833245*I,
        0.00494602046863013 -   0.00325922937957204*I,
        0.00457437001988662 -   0.00468450876004249*I,
          0.003323129618157 -   0.00362219897727821*I,
        0.00383263814195228 -   0.00415838147145455*I,
        0.00300832770761044 -    0.0025232535398731*I,
        0.00359443831265403 -    0.0020580058659592*I,
        0.00236478550332468 -  0.000776610384360161*I,
        0.00220294672839155 +  2.88422244171861e-05*I,
        0.00103138166556304 +  0.000356635975839707*I,
       0.000429677977956195 +  0.000702698094317228*I,
       0.000141186123348109 +  0.000308760666557993*I,
       1.76161843819631e-05 +  3.85241069693263e-05*I,
        0.00471000774614879 +   0.00321951527165378*I,
         0.0189122705005425 +    0.0129253478144194*I,
        0.00143721720296166 +   0.00959129465377047*I,
        -0.0141785792092596 +    0.0248899171756221*I,
        -0.0121752886109488 +    0.0152445101105612*I,
         -0.034473296266092 +    0.0347437886898533*I,
        -0.0146134526780112 +    0.0196037409478134*I,
        -0.0229483461879859 +    0.0415158853611465*I,
        -0.0033916553772744 +    0.0221807639839301*I,
        0.00984625500059017 +    0.0445112101146303*I,
             0.011015403441 +    0.0226759354165135*I,
          0.033523817026556 +    0.0433990886521313*I,
         0.0151100594624086 +     0.020990642771451*I,
         0.0261833595590324 +    0.0382376846087666*I,
        0.00518512881431605 +    0.0172863054149115*I,
       -0.00509713173069031 +    0.0295477309013482*I,
       -0.00126908276981291 +   0.00735925061620365*I,
       -0.00126908276981291 -   0.00735925061620365*I,
       -0.00509713173069031 -    0.0295477309013482*I,
        0.00518512881431605 -    0.0172863054149115*I,
         0.0261833595590324 -    0.0382376846087666*I,
         0.0151100594624086 -     0.020990642771451*I,
          0.033523817026556 -    0.0433990886521313*I,
             0.011015403441 -    0.0226759354165135*I,
        0.00984625500059017 -    0.0445112101146303*I,
        -0.0033916553772744 -    0.0221807639839301*I,
        -0.0229483461879859 -    0.0415158853611465*I,
        -0.0146134526780112 -    0.0196037409478134*I,
         -0.034473296266092 -    0.0347437886898533*I,
        -0.0121752886109488 -    0.0152445101105612*I,
        -0.0141785792092596 -    0.0248899171756221*I,
        0.00143721720296166 -   0.00959129465377047*I,
         0.0189122705005425 -    0.0129253478144194*I,
        0.00471000774614879 -   0.00321951527165378*I,
       1.76161843819631e-05 -  3.85241069693263e-05*I,
       0.000141186123348109 -  0.000308760666557993*I,
       0.000429677977956195 -  0.000702698094317228*I,
        0.00103138166556304 -  0.000356635975839707*I,
        0.00220294672839155 -   2.8842224417186e-05*I,
        0.00236478550332468 +  0.000776610384360161*I,
        0.00359443831265403 +    0.0020580058659592*I,
        0.00300832770761044 +   0.00252325353987311*I,
        0.00383263814195228 +   0.00415838147145455*I,
          0.003323129618157 +   0.00362219897727822*I,
        0.00457437001988662 +   0.00468450876004249*I,
        0.00494602046863012 +   0.00325922937957204*I,
        0.00794998805226457 +   0.00330086725833245*I,
        0.00854978166679528 +   0.00180605696970404*I,
         0.0124758235986958 +    0.0012200695988487*I,
         0.0155954521142025 +  0.000467093980630541*I,
           1.00578402514809 +  5.83210136741118e-05*I};
           
    /* Matlab code to generate result_exact:
    e = 0.13;
    kappa = -1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    r=-kappa*conj(q);
    ms=3;% 
    S=zeros(2,2*ms,D);
    for n=1:D
        if q(n) == 0
             B11 = 1; B12 = 0; B21 = 0; B22 = 1;
        else
             B11=cosh(0.25*e*q(n)^(1/2)*r(n)^(1/2));
             B12= (sinh(0.25*e*(q(n)*r(n))^(1/2))*(q(n)*r(n))^(1/2))/r(n);
             B21= -kappa*conj(B12);
             B22=(B11);
        end
        % Scattering matrix at n
        S(:,:,n)=[[B11^4 + (2*B12*B21*B11^2)/3 - (B12^2*B21^2)/3,(8*B12*B21*B11^2)/3 + (8*B12*B21*B22*B11)/3,(4*B12^2*B21^2)/3 + B12*B21*B22^2 - (2*B11*B12*B21*B22)/3 - (B11^2*B12*B21)/3],...
            [-(B12*(- 3*B11^3 + B22*B11^2 - 3*B12*B21*B11 + B12*B21*B22))/3,(B12*(4*B11^2*B22 + 4*B11*B22^2 + 4*B12*B21*B11 + 4*B12*B21*B22))/3,-(B12*(- 3*B22^3 + B11*B22^2 - 3*B12*B21*B22 + B11*B12*B21))/3];...
            [-(B21*(- 3*B11^3 + B22*B11^2 - 3*B12*B21*B11 + B12*B21*B22))/3,(B21*(4*B11^2*B22 + 4*B11*B22^2 + 4*B12*B21*B11 + 4*B12*B21*B22))/3,-(B21*(- 3*B22^3 + B11*B22^2 - 3*B12*B21*B22 + B11*B12*B21))/3],...
            [ B11^2*B12*B21 - (2*B22*B11*B12*B21)/3 + (4*B12^2*B21^2)/3 - (B22^2*B12*B21)/3, (8*B12*B21*B22^2)/3 + (8*B11*B12*B21*B22)/3, B22^4 + (2*B12*B21*B22^2)/3 - (B12^2*B21^2)/3]];
    end
    % Multiply scattering matrices
    temp=S(:,:,D);
    m=ms;
    for n=D-1:-1:1
        temp=[conv(temp(1,1:m),S(1,1:ms,n))+conv(temp(1,m+1:end),S(2,1:ms,n)),...
            conv(temp(1,1:m),S(1,ms+1:end,n))+conv(temp(1,m+1:end),S(2,ms+1:end,n));...
            conv(temp(2,1:m),S(1,1:ms,n))+conv(temp(2,m+1:end),S(2,1:ms,n)),...
            conv(temp(2,1:m),S(1,ms+1:end,n))+conv(temp(2,m+1:end),S(2,ms+1:end,n))];
     m=length(temp)/2;
    end
    format long g; result_exact = [temp(1,1:m) temp(1,m+1:end) temp(2,1:m) temp(2,m+1:end)].'
    */

    for (i=0; i<D; i++)
        q[i] = 0.4*cos(i+1) + 0.5*I*sin(0.3*(i+1));

    // without normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, NULL, nse_discretization_2SPLIT4B);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, -1, result, &deg, &W, nse_discretization_2SPLIT4B);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W != 0) // it really should be zero in this case!
        return E_TEST_FAILED;
    scl = POW(2.0, W);
    for (i=0; i<4*(deg+1); i++)
        result[i] *= scl;
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    return SUCCESS;
}

INT main()
{
    if (nse_fscatter_test_defocusing_2split4B() != SUCCESS)
        return EXIT_FAILURE;

    return SUCCESS;
}
