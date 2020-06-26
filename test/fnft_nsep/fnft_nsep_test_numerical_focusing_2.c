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
 * Sander Wahls (TU Delft) 2020.
 * Shrinivas Chimmalgi (TU Delft) 2020.
 */

#define FNFT_ENABLE_SHORT_NAMES
#include "fnft.h"
#include "fnft_nsep.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"
#ifdef DEBUG
#include <stdio.h>
#endif

// The test signal was build using the following Matlab code.

/*
% Code to build specific genus-2 solutions for focusing NSE
% Built using Eq.20 of A. O. Smirnov et. al, "On a periodic solution of the
% focusing nonlinear Schrödinger equation," arXiv:1407.7974v1, 2014. Spines
% are shown in Figure 1.
clear
clc
%%
l0 = 0;
a = 1; b = 2; c = 5;
Z1 = pi/4; Z2 = -pi/4; %phases
K1 = -l0; K2 = 0;%wavenumber/frequency
lambda = [l0+a*1i,l0+b*1i,l0+c*1i];
%% Numerical
f = @(t) 1./sqrt((t-a^2).*(b^2-t).*(c^2-t));
Ap = integral(f,a^2,b^2);
f = @(t) 1./sqrt(t.*(a^2-t).*(b^2-t).*(c^2-t));
Am = integral(f,0,a^2);
f = @(t) 1./sqrt((t-a^2).*(t-b^2).*(c^2-t));
Bp = integral(f,b^2,c^2);
f = @(t) 1./sqrt(t.*(t-a^2).*(b^2-t).*(c^2-t));
Bm = integral(f,a^2,b^2);
f = @(t) 1./sqrt(t.*(t-a^2).*(t-b^2).*(t-c^2));
Bm1 = integral(f,c^2,Inf);
%%
bm = Bm/Am;
bp = Bp/Ap;
kappa1 = 4/Am;
k = 2/Ap;
kappa2 = 8*l0/Ap;
del = Bm1/Am;
f = @(t) t./sqrt(t.*(a^2-t).*(b^2-t).*(c^2-t));
Dm = 0.5*integral(f,0,a^2);
f = @(t) (t./sqrt(t.*(t-a^2).*(t-b^2).*(t-c^2)))-(1./t);
Fm = 0.5*integral(f,c^2,Inf);
K0 = 1i*c*exp(Dm*del-Fm);
%%
T = [0, Ap];  % T(1) is the beginning of the period, T(2) is the end

D = 2^8+1;        % number of samples
kappa = +1;     % focusing nonlinear Schroedinger equation

%%% Setup the signal %%%
X = linspace(T(1),T(2),D);


%% Building the signal

t = 0;
q = zeros(1,D);
for j = 1:1:length(X)
    x = X(j);
    u1 = kappa1*t+2*Z1; u2 = k*x+kappa2*t+2*Z2;
    q(j) = -2i*K0*(H(u1+1i*del,u2+1,bp,bm))*(exp(2i*K1*x+2i*K2*t))/(H(u1,u2,bp,bm));
end

%%
function s = H(u1,u2,bp,bm)

s = v3(u1,2i*bm)*v3(u2,2i*bp)+v2(u1,2i*bm)*v3(u2,2i*bp)+...
    v3(u1,2i*bm)*v2(u2,2i*bp)-v2(u1,2i*bm)*v2(u2,2i*bp);
end


function s = v2(u,b)
low_lim = 1e-50;
up_lim = 1e100;
h = exp(pi*1i*b);
s = 0;
m = 1;
temp = 1;
while abs(temp)>low_lim && abs(temp)<up_lim
    temp = 2*(h^((m-0.5)^2))*cos((2*m-1)*pi*u);
    s = s+temp;
    m = m+1;
end
end

function s = v3(u,b)
low_lim = 1e-50;
up_lim = 1e100;
h = exp(pi*1i*b);
s = 1;
m = 1;
temp = 1;
while abs(temp)>low_lim && abs(temp)<up_lim
    temp = 2*(h^(m^2))*cos(2*m*pi*u);
    s = s+temp;
    m = m+1;
end
end
*/

INT fnft_nsep_test_numerical_focusing_2()
{
    INT ret_code = SUCCESS;

    COMPLEX q[257] = {
      4.342615484209547e+00 + 9.303042232295938e-01*I,
      4.376264025488326e+00 + 9.299262942389809e-01*I,
      4.410337956208852e+00 + 9.295435874185952e-01*I,
      4.444823341542983e+00 + 9.291562592796946e-01*I,
      4.479705659973817e+00 + 9.287644729230267e-01*I,
      4.514969791258486e+00 + 9.283683981740255e-01*I,
      4.550600004981735e+00 + 9.279682117113738e-01*I,
      4.586579949772434e+00 + 9.275640971881201e-01*I,
      4.622892643258073e+00 + 9.271562453445065e-01*I,
      4.659520462834854e+00 + 9.267448541116349e-01*I,
      4.696445137333549e+00 + 9.263301287050749e-01*I,
      4.733647739663573e+00 + 9.259122817074811e-01*I,
      4.771108680519800e+00 + 9.254915331392770e-01*I,
      4.808807703238547e+00 + 9.250681105164315e-01*I,
      4.846723879890605e+00 + 9.246422488943393e-01*I,
      4.884835608700487e+00 + 9.242141908968083e-01*I,
      4.923120612881922e+00 + 9.237841867291386e-01*I,
      4.961555940980000e+00 + 9.233524941742778e-01*I,
      5.000117968810484e+00 + 9.229193785710418e-01*I,
      5.038782403086276e+00 + 9.224851127733836e-01*I,
      5.077524286820024e+00 + 9.220499770897107e-01*I,
      5.116318006590354e+00 + 9.216142592012754e-01*I,
      5.155137301757020e+00 + 9.211782540586723e-01*I,
      5.193955275707504e+00 + 9.207422637555210e-01*I,
      5.232744409214154e+00 + 9.203065973784420e-01*I,
      5.271476575976768e+00 + 9.198715708324896e-01*I,
      5.310123060420705e+00 + 9.194375066412483e-01*I,
      5.348654577814929e+00 + 9.190047337208730e-01*I,
      5.387041296768062e+00 + 9.185735871274235e-01*I,
      5.425252864153208e+00 + 9.181444077769146e-01*I,
      5.463258432504553e+00 + 9.177175421376083e-01*I,
      5.501026689919683e+00 + 9.172933418941610e-01*I,
      5.538525892492248e+00 + 9.168721635833506e-01*I,
      5.575723899289113e+00 + 9.164543682012247e-01*I,
      5.612588209875022e+00 + 9.160403207816366e-01*I,
      5.649086004376043e+00 + 9.156303899462677e-01*I,
      5.685184186060393e+00 + 9.152249474263707e-01*I,
      5.720849426402195e+00 + 9.148243675566323e-01*I,
      5.756048212579819e+00 + 9.144290267416860e-01*I,
      5.790746897346238e+00 + 9.140393028959856e-01*I,
      5.824911751194041e+00 + 9.136555748579066e-01*I,
      5.858509016722668e+00 + 9.132782217791092e-01*I,
      5.891504965100041e+00 + 9.129076224903836e-01*I,
      5.923865954495263e+00 + 9.125441548453513e-01*I,
      5.955558490343460e+00 + 9.121881950435923e-01*I,
      5.986549287288405e+00 + 9.118401169349266e-01*I,
      6.016805332633221e+00 + 9.115002913067554e-01*I,
      6.046293951114585e+00 + 9.111690851565423e-01*I,
      6.074982870801374e+00 + 9.108468609516602e-01*I,
      6.102840289904887e+00 + 9.105339758790012e-01*I,
      6.129834944274746e+00 + 9.102307810868878e-01*I,
      6.155936175342400e+00 + 9.099376209219515e-01*I,
      6.181113998263188e+00 + 9.096548321637872e-01*I,
      6.205339169997965e+00 + 9.093827432602801e-01*I,
      6.228583257066869e+00 + 9.091216735666229e-01*I,
      6.250818702700659e+00 + 9.088719325910907e-01*I,
      6.272018893109740e+00 + 9.086338192507312e-01*I,
      6.292158222587105e+00 + 9.084076211401523e-01*I,
      6.311212157159504e+00 + 9.081936138166108e-01*I,
      6.329157296501142e+00 + 9.079920601046205e-01*I,
      6.345971433825905e+00 + 9.078032094232632e-01*I,
      6.361633613478119e+00 + 9.076272971393483e-01*I,
      6.376124185947528e+00 + 9.074645439495032e-01*I,
      6.389424860042221e+00 + 9.073151552941839e-01*I,
      6.401518751962999e+00 + 9.071793208064917e-01*I,
      6.412390431034529e+00 + 9.070572137985315e-01*I,
      6.422025961862502e+00 + 9.069489907879216e-01*I,
      6.430412942701484e+00 + 9.068547910668546e-01*I,
      6.437540539835633e+00 + 9.067747363159493e-01*I,
      6.443399517793282e+00 + 9.067089302648872e-01*I,
      6.447982265236975e+00 + 9.066574584016261e-01*I,
      6.451282816392218e+00 + 9.066203877317203e-01*I,
      6.453296867901169e+00 + 9.065977665890229e-01*I,
      6.454021791011417e+00 + 9.065896244987890e-01*I,
      6.453456639034603e+00 + 9.065959720939012e-01*I,
      6.451602150034862e+00 + 9.066168010846724e-01*I,
      6.448460744732753e+00 + 9.066520842823929e-01*I,
      6.444036519635776e+00 + 9.067017756764800e-01*I,
      6.438335235432470e+00 + 9.067658105648364e-01*I,
      6.431364300712131e+00 + 9.068441057367032e-01*I,
      6.423132751097013e+00 + 9.069365597070388e-01*I,
      6.413651223897923e+00 + 9.070430530011800e-01*I,
      6.402931928427093e+00 + 9.071634484882770e-01*I,
      6.390988612124187e+00 + 9.072975917617561e-01*I,
      6.377836522671912e+00 + 9.074453115648218e-01*I,
      6.363492366296875e+00 + 9.076064202588093e-01*I,
      6.347974262468863e+00 + 9.077807143319865e-01*I,
      6.331301695227565e+00 + 9.079679749462348e-01*I,
      6.313495461379732e+00 + 9.081679685188825e-01*I,
      6.294577615821948e+00 + 9.083804473368198e-01*I,
      6.274571414254164e+00 + 9.086051501999205e-01*I,
      6.253501253557436e+00 + 9.088418030907010e-01*I,
      6.231392610115260e+00 + 9.090901198670717e-01*I,
      6.208271976362121e+00 + 9.093498029750036e-01*I,
      6.184166795844837e+00 + 9.096205441778972e-01*I,
      6.159105397082531e+00 + 9.099020252994446e-01*I,
      6.133116926509263e+00 + 9.101939189767945e-01*I,
      6.106231280779902e+00 + 9.104958894208708e-01*I,
      6.078479038714388e+00 + 9.108075931807504e-01*I,
      6.049891393148962e+00 + 9.111286799090882e-01*I,
      6.020500082954314e+00 + 9.114587931256648e-01*I,
      5.990337325471102e+00 + 9.117975709762496e-01*I,
      5.959435749602277e+00 + 9.121446469840826e-01*I,
      5.927828329789723e+00 + 9.124996507914294e-01*I,
      5.895548321089684e+00 + 9.128622088887894e-01*I,
      5.862629195547806e+00 + 9.132319453295114e-01*I,
      5.829104580060231e+00 + 9.136084824277150e-01*I,
      5.795008195892247e+00 + 9.139914414375987e-01*I,
      5.760373800010790e+00 + 9.143804432123717e-01*I,
      5.725235128371701e+00 + 9.147751088412324e-01*I,
      5.689625841286889e+00 + 9.151750602629868e-01*I,
      5.653579470981233e+00 + 9.155799208550702e-01*I,
      5.617129371433514e+00 + 9.159893159969184e-01*I,
      5.580308670580579e+00 + 9.164028736067927e-01*I,
      5.543150224949136e+00 + 9.168202246513411e-01*I,
      5.505686576765246e+00 + 9.172410036273303e-01*I,
      5.467949913577651e+00 + 9.176648490151412e-01*I,
      5.429972030417931e+00 + 9.180914037037756e-01*I,
      5.391784294507761e+00 + 9.185203153872494e-01*I,
      5.353417612511642e+00 + 9.189512369324010e-01*I,
      5.314902400322270e+00 + 9.193838267182515e-01*I,
      5.276268555355282e+00 + 9.198177489471798e-01*I,
      5.237545431320409e+00 + 9.202526739282845e-01*I,
      5.198761815427210e+00 + 9.206882783334024e-01*I,
      5.159945907975498e+00 + 9.211242454263412e-01*I,
      5.121125304273254e+00 + 9.215602652659719e-01*I,
      5.082326978818367e+00 + 9.219960348838920e-01*I,
      5.043577271674784e+00 + 9.224312584374477e-01*I,
      5.004901876968716e+00 + 9.228656473389376e-01*I,
      4.966325833426277e+00 + 9.232989203618907e-01*I,
      4.927873516870455e+00 + 9.237308037253410e-01*I,
      4.889568634592369e+00 + 9.241610311570436e-01*I,
      4.851434221509606e+00 + 9.245893439366246e-01*I,
      4.813492638022846e+00 + 9.250154909196573e-01*I,
      4.775765569480754e+00 + 9.254392285436681e-01*I,
      4.738274027162805e+00 + 9.258603208170997e-01*I,
      4.701038350689523e+00 + 9.262785392922410e-01*I,
      4.664078211770017e+00 + 9.266936630231343e-01*I,
      4.627412619197610e+00 + 9.271054785094668e-01*I,
      4.591059925005391e+00 + 9.275137796274351e-01*I,
      4.555037831695127e+00 + 9.279183675485492e-01*I,
      4.519363400454791e+00 + 9.283190506473392e-01*I,
      4.484053060281895e+00 + 9.287156443988833e-01*I,
      4.449122617932246e+00 + 9.291079712670671e-01*I,
      4.414587268616156e+00 + 9.294958605844486e-01*I,
      4.380461607366730e+00 + 9.298791484245699e-01*I,
      4.346759641007713e+00 + 9.302576774675440e-01*I,
      4.313494800651175e+00 + 9.306312968596800e-01*I,
      4.280679954658289e+00 + 9.309998620679147e-01*I,
      4.248327421999467e+00 + 9.313632347297551e-01*I,
      4.216448985953124e+00 + 9.317212824994189e-01*I,
      4.185055908085447e+00 + 9.320738788908189e-01*I,
      4.154158942456490e+00 + 9.324209031180030e-01*I,
      4.123768350001064e+00 + 9.327622399336365e-01*I,
      4.093893913035754e+00 + 9.330977794660600e-01*I,
      4.064544949846460e+00 + 9.334274170554513e-01*I,
      4.035730329313571e+00 + 9.337510530895565e-01*I,
      4.007458485534871e+00 + 9.340685928394532e-01*I,
      3.979737432408872e+00 + 9.343799462957540e-01*I,
      3.952574778143970e+00 + 9.346850280056460e-01*I,
      3.925977739661346e+00 + 9.349837569111207e-01*I,
      3.899953156862023e+00 + 9.352760561887322e-01*I,
      3.874507506730804e+00 + 9.355618530911864e-01*I,
      3.849646917252107e+00 + 9.358410787910431e-01*I,
      3.825377181114811e+00 + 9.361136682267874e-01*I,
      3.801703769185345e+00 + 9.363795599515056e-01*I,
      3.778631843730093e+00 + 9.366386959843742e-01*I,
      3.756166271370099e+00 + 9.368910216651580e-01*I,
      3.734311635752727e+00 + 9.371364855118864e-01*I,
      3.713072249926558e+00 + 9.373750390818615e-01*I,
      3.692452168407382e+00 + 9.376066368361397e-01*I,
      3.672455198924429e+00 + 9.378312360075985e-01*I,
      3.653084913837451e+00 + 9.380487964727061e-01*I,
      3.634344661216380e+00 + 9.382592806270780e-01*I,
      3.616237575576482e+00 + 9.384626532649029e-01*I,
      3.598766588262906e+00 + 9.386588814623074e-01*I,
      3.581934437479559e+00 + 9.388479344647153e-01*I,
      3.565743677958046e+00 + 9.390297835782480e-01*I,
      3.550196690263226e+00 + 9.392044020652098e-01*I,
      3.535295689732661e+00 + 9.393717650436813e-01*I,
      3.521042735047880e+00 + 9.395318493912506e-01*I,
      3.507439736435942e+00 + 9.396846336528958e-01*I,
      3.494488463500285e+00 + 9.398300979530319e-01*I,
      3.482190552680333e+00 + 9.399682239117275e-01*I,
      3.470547514339688e+00 + 9.400989945650906e-01*I,
      3.459560739483085e+00 + 9.402223942898291e-01*I,
      3.449231506102544e+00 + 9.403384087319701e-01*I,
      3.439560985153462e+00 + 9.404470247397401e-01*I,
      3.430550246161467e+00 + 9.405482303005903e-01*I,
      3.422200262461110e+00 + 9.406420144823583e-01*I,
      3.414511916067498e+00 + 9.407283673785491e-01*I,
      3.407486002182108e+00 + 9.408072800577298e-01*I,
      3.401123233334015e+00 + 9.408787445170161e-01*I,
      3.395424243157793e+00 + 9.409427536396396e-01*I,
      3.390389589809335e+00 + 9.409993011565825e-01*I,
      3.386019759020786e+00 + 9.410483816122669e-01*I,
      3.382315166795724e+00 + 9.410899903342808e-01*I,
      3.379276161745633e+00 + 9.411241234071391e-01*I,
      3.376903027068595e+00 + 9.411507776500572e-01*I,
      3.375195982171058e+00 + 9.411699505987362e-01*I,
      3.374155183933335e+00 + 9.411816404911494e-01*I,
      3.373780727619415e+00 + 9.411858462573246e-01*I,
      3.374072647431474e+00 + 9.411825675131141e-01*I,
      3.375030916709330e+00 + 9.411718045579558e-01*I,
      3.376655447774964e+00 + 9.411535583766202e-01*I,
      3.378946091421977e+00 + 9.411278306449422e-01*I,
      3.381902636049820e+00 + 9.410946237395482e-01*I,
      3.385524806442364e+00 + 9.410539407515712e-01*I,
      3.389812262190313e+00 + 9.410057855043736e-01*I,
      3.394764595756742e+00 + 9.409501625752718e-01*I,
      3.400381330184999e+00 + 9.408870773212801e-01*I,
      3.406661916448007e+00 + 9.408165359088850e-01*I,
      3.413605730437950e+00 + 9.407385453478535e-01*I,
      3.421212069595215e+00 + 9.406531135290970e-01*I,
      3.429480149175401e+00 + 9.405602492665969e-01*I,
      3.438409098153167e+00 + 9.404599623434153e-01*I,
      3.447997954761640e+00 + 9.403522635617908e-01*I,
      3.458245661666161e+00 + 9.402371647973471e-01*I,
      3.469151060771124e+00 + 9.401146790574207e-01*I,
      3.480712887658766e+00 + 9.399848205435233e-01*I,
      3.492929765658880e+00 + 9.398476047179465e-01*I,
      3.505800199548527e+00 + 9.397030483745269e-01*I,
      3.519322568881039e+00 + 9.395511697135712e-01*I,
      3.533495120943831e+00 + 9.393919884209522e-01*I,
      3.548315963344794e+00 + 9.392255257513753e-01*I,
      3.563783056227402e+00 + 9.390518046158169e-01*I,
      3.579894204115019e+00 + 9.388708496731248e-01*I,
      3.596647047385367e+00 + 9.386826874257751e-01*I,
      3.614039053376569e+00 + 9.384873463197658e-01*I,
      3.632067507126807e+00 + 9.382848568486268e-01*I,
      3.650729501750219e+00 + 9.380752516615131e-01*I,
      3.670021928452394e+00 + 9.378585656753473e-01*I,
      3.689941466189617e+00 + 9.376348361909653e-01*I,
      3.710484570976832e+00 + 9.374041030132041e-01*I,
      3.731647464850285e+00 + 9.371664085748694e-01*I,
      3.753426124491830e+00 + 9.369217980645047e-01*I,
      3.775816269522950e+00 + 9.366703195578683e-01*I,
      3.798813350477841e+00 + 9.364120241530152e-01*I,
      3.822412536466128e+00 + 9.361469661088657e-01*I,
      3.846608702537258e+00 + 9.358752029871243e-01*I,
      3.871396416760052e+00 + 9.355967957973977e-01*I,
      3.896769927032560e+00 + 9.353118091453425e-01*I,
      3.922723147639052e+00 + 9.350203113836502e-01*I,
      3.949249645572780e+00 + 9.347223747656697e-01*I,
      3.976342626645072e+00 + 9.344180756014204e-01*I,
      4.003994921403389e+00 + 9.341074944157588e-01*I,
      4.032198970883042e+00 + 9.337907161084105e-01*I,
      4.060946812219545e+00 + 9.334678301155649e-01*I,
      4.090230064150947e+00 + 9.331389305727107e-01*I,
      4.120039912441849e+00 + 9.328041164783476e-01*I,
      4.150367095263449e+00 + 9.324634918581911e-01*I,
      4.181201888566501e+00 + 9.321171659294610e-01*I,
      4.212534091486893e+00 + 9.317652532648002e-01*I,
      4.244353011826203e+00 + 9.314078739553496e-01*I,
      4.276647451652630e+00 + 9.310451537724775e-01*I,
      4.309405693070493e+00 + 9.306772243276101e-01*I,
      4.342615484209547e+00 + 9.303042232295938e-01*I
    };
    const REAL T[2] = {0, 6.628591515456010e-01};

    const COMPLEX mainspec_exact[6] = {-5*I, -2*I, -1*I, 1*I, 2*I, 5*I};
    const UINT K_exact = sizeof(mainspec_exact) / sizeof(mainspec_exact[0]);

    COMPLEX *mainspec = NULL;
    COMPLEX *spines = NULL;

    const UINT D = sizeof(q) / sizeof(q[0]);
    UINT K = D;
    UINT M = 0;
    UINT points_per_spine = 250;
    UINT K_spine = K*points_per_spine;
    mainspec = malloc(K * sizeof(COMPLEX));
    spines = malloc(K_spine * sizeof(COMPLEX));
    if (mainspec == NULL || spines == NULL) {
        ret_code = E_NOMEM;
        goto leave_fun;
    }

    // Compute main and auxiliary spectrum
    fnft_nsep_opts_t opts = fnft_nsep_default_opts();
    opts.filtering = fnft_nsep_filt_MANUAL;
    opts.bounding_box[0] = -1;
    opts.bounding_box[1] = 1;
    opts.bounding_box[2] = -10;
    opts.bounding_box[3] = 10;
    ret_code = fnft_nsep(D, q, T, &K, mainspec, &M, NULL, NULL, +1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check that the main and auxiliary spectrum are correct
    REAL dist = misc_hausdorff_dist(K_exact, mainspec_exact, K, mainspec);
    if (dist > 4e-4) {
#ifdef DEBUG
        printf("dist = %e\n", dist);
#endif
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    // Compute spines
    opts.points_per_spine = points_per_spine;
    ret_code = fnft_nsep(D, q, T, &K_spine, spines, &M, NULL, NULL, +1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check that all found points are on one of the three spines. The ok_flags
    // are turned on when we find a point "close" to the center of a spines.
    const REAL tol = 0.00004;
    INT ok_flag[3] = {0, 0, 0};
    for (UINT k=0; k<K_spine; k++) {

        const COMPLEX lam = spines[k];

        // Spines are imaginary in this example
        if (FABS(CREAL(lam)) > 1000*EPSILON)
            return E_TEST_FAILED;

        const REAL lam_i = CIMAG(lam);

        // Spine [-5*I, -2*I]
        if (lam_i>-4.5 && lam_i<-2.5)
            ok_flag[0] = 1;
        if (lam_i>=-5-tol && lam_i<=-2+tol)
            continue;

        // Spine [-I, I]
        if (FABS(lam_i) < 0.5)
            ok_flag[1] = 1;
        if (FABS(lam_i) <= 1+tol)
            continue;

        // Spine [2*I, 5*I]
        if (lam_i>2.5 && lam_i<4.5)
            ok_flag[2] = 1;
        if (lam_i>=2-tol && lam_i<=5+tol)
            continue;

        // lam was in none of the spines
#ifdef DEBUG
        printf("lam = %e + %ej\n", CREAL(lam), CIMAG(lam));
#endif
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

    // Check if points "close" to the center of each spine were found
    if (!ok_flag[0]) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (!ok_flag[1]) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    if (!ok_flag[2]) {
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    free(mainspec);
    free(spines);
    return SUCCESS;
}

int main()
{
    if (fnft_nsep_test_numerical_focusing_2() == SUCCESS)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
