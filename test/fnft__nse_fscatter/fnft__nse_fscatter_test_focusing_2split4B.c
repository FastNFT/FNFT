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

static INT nse_fscatter_test_focusing_2split4B()
{
    UINT i, D = 8, deg;
    INT W;
    REAL scl;
    INT ret_code;
    const REAL eps_t = 0.13;
    COMPLEX q[8];
    COMPLEX result[3*4*8];
    COMPLEX result_exact[68] = {
          0.994247576912421 +  5.77391551056385e-05*I,
         -0.015430195607459 +  0.000461392025535083*I,
        -0.0121488006133603 +   0.00119974633455092*I,
       -0.00818373882110172 +   0.00176808026557291*I,
       -0.00759531054331097 +    0.0032419342835477*I,
       -0.00461764064622704 +   0.00317948346909759*I,
       -0.00428799881282084 +   0.00460177722522967*I,
       -0.00307137405587193 +   0.00354074310921268*I,
       -0.00361551687953187 +   0.00410399515855735*I,
       -0.00282129722687509 +   0.00249458343394973*I,
       -0.00345092101380092 +   0.00206293179565897*I,
        -0.0022620945739206 +   0.00080405342304031*I,
       -0.00214534379123276 +  6.90912975446876e-06*I,
       -0.00100657703177414 -  0.000323522309188338*I,
      -0.000420836189635619 -  0.000684649970948553*I,
      -0.000139172594699345 -  0.000304090350539432*I,
      -1.74283363456518e-05 -  3.80815220504591e-05*I,
        0.00465538328234486 +   0.00318351880542936*I,
         0.0185500633088764 +     0.012687279800172*I,
        0.00071949215768491 +   0.00904870991946628*I,
        -0.0146663485952178 +    0.0240524428991485*I,
        -0.0121900955597852 +    0.0139728142014193*I,
        -0.0341610379123053 +    0.0331615329471449*I,
        -0.0138742265701376 +    0.0175851650470893*I,
        -0.0223713500442326 +    0.0392889017194499*I,
       -0.00302183839113178 +    0.0196573663888199*I,
        0.00982150070922305 +    0.0419716031721442*I,
         0.0105456591749782 +    0.0200682713858596*I,
         0.0330491701937163 +    0.0410290563612796*I,
         0.0146000884630374 +    0.0188339341742728*I,
         0.0260853256405769 +    0.0365883674001791*I,
        0.00535913442350011 +    0.0160588970706718*I,
       -0.00499129584927267 +    0.0289944669938156*I,
       -0.00125296552447965 +   0.00727601002158376*I,
        0.00125296552447965 +   0.00727601002158376*I,
        0.00499129584927267 +    0.0289944669938156*I,
       -0.00535913442350011 +    0.0160588970706718*I,
        -0.0260853256405769 +    0.0365883674001791*I,
        -0.0146000884630374 +    0.0188339341742728*I,
        -0.0330491701937163 +    0.0410290563612796*I,
        -0.0105456591749782 +    0.0200682713858596*I,
       -0.00982150070922306 +    0.0419716031721442*I,
        0.00302183839113178 +    0.0196573663888199*I,
         0.0223713500442326 +    0.0392889017194499*I,
         0.0138742265701376 +    0.0175851650470893*I,
         0.0341610379123053 +    0.0331615329471449*I,
         0.0121900955597852 +    0.0139728142014193*I,
         0.0146663485952178 +    0.0240524428991485*I,
       -0.00071949215768491 +   0.00904870991946628*I,
        -0.0185500633088764 +     0.012687279800172*I,
       -0.00465538328234486 +   0.00318351880542936*I,
      -1.74283363456518e-05 +  3.80815220504591e-05*I,
      -0.000139172594699345 +  0.000304090350539432*I,
      -0.000420836189635619 +  0.000684649970948553*I,
       -0.00100657703177414 +  0.000323522309188338*I,
       -0.00214534379123276 -  6.90912975446876e-06*I,
        -0.0022620945739206 -   0.00080405342304031*I,
       -0.00345092101380092 -   0.00206293179565897*I,
       -0.00282129722687509 -   0.00249458343394973*I,
       -0.00361551687953187 -   0.00410399515855735*I,
       -0.00307137405587193 -   0.00354074310921268*I,
       -0.00428799881282084 -   0.00460177722522967*I,
       -0.00461764064622704 -   0.00317948346909759*I,
       -0.00759531054331097 -   0.00324193428354771*I,
       -0.00818373882110172 -   0.00176808026557291*I,
        -0.0121488006133603 -   0.00119974633455092*I,
         -0.015430195607459 -  0.000461392025535083*I,
          0.994247576912421 -  5.77391551056386e-05*I};
          
    /* Matlab code to generate result_exact:
    e = 0.13;
    kappa = 1; D=8; q = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
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
        q[i] = 0.4*COS(i+1) + 0.5*I*SIN(0.3*(i+1));

    // without normalization
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, NULL, nse_discretization_2SPLIT4B);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (misc_rel_err(4*(deg+1), result, result_exact) > 100*EPSILON)
        return E_TEST_FAILED;

    // with normalization
    ret_code = nse_fscatter(D, q, eps_t, +1, result, &deg, &W,nse_discretization_2SPLIT4B);
    if (ret_code != SUCCESS)
        return E_SUBROUTINE(ret_code);
    if (W == 0)
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
    if (nse_fscatter_test_focusing_2split4B() != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
