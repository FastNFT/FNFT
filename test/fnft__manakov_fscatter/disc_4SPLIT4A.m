%% Matlab file to get the polynomial coefficients for the 4SPLIT4A
% CF24 method with 4th order splitting
clear
close all

syms s1 s2 l kappa h z

AE = [1/z^8,0,0;0,z^8,0;0,0,z^8];
BE = sym('BE',[3,3]);


% we have CE = expm(A+B)
% Let G_approx be approximation of CE after application of splitting
% scheme
% AE_8 = expm(A*h/8)
% BE_2 = expm(2B*h/4) = expm(B*h/2)
AE_8 = [1/z,0,0;0,z,0;0,0,z];
BE_2 = sym('BE_2',[3,3]);


G_approx = ((4/3)*AE_8*(BE_2)*(AE_8*AE_8)*(BE_2)*AE_8-...
        (1/3)*(AE_8*AE_8)*BE*(AE_8*AE_8));

% Dividing throughout by z^-4
G_approx_pos = expand(G_approx*z^4);

% We can now look at coefficients of individual polynomials. Here c gives
% the coefficients of the element specified by (i,j) in CE_approx_pos(i,j)
% and t gives which power of z the coefficients belong to
[c11,t11] = coeffs(G_approx_pos(1,1),z);
[c12,t12] = coeffs(G_approx_pos(1,2),z);
[c13,t13] = coeffs(G_approx_pos(1,3),z);

[c21,t21] = coeffs(G_approx_pos(2,1),z);
[c22,t22] = coeffs(G_approx_pos(2,2),z);
[c23,t23] = coeffs(G_approx_pos(2,3),z);

[c31,t31] = coeffs(G_approx_pos(3,1),z);
[c32,t32] = coeffs(G_approx_pos(3,2),z);
[c33,t33] = coeffs(G_approx_pos(3,3),z);

% now look at all the c's to see what the matrices for each coefficient
% should be:

matz0 = [c11(2)     0           0;
        0           0           0;
        0           0           0];
    
matz2 = [0          c12(3)      c13(3);
        c21(3)      0           0;
        c31(3)      0           0];

matz4 = [c11(1)     c12(2)      c13(2);
        c21(2)      c22(2)      c23(2);
        c31(2)      c32(2)      c33(2)];    
    
matz6 = [0          c12(1)      c13(1);
        c21(1)      0           0;
        c31(1)      0           0];    
    
matz8 = [0          0           0;
        0           c22(1)      c23(1);
        0           c32(1)      c33(1)];    
    
    
%% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
c1 = 0.5-sqrt(3)/6;
c2 = 0.5+sqrt(3)/6;
a1  = 0.25 + sqrt(3)/6;
a2  = 0.25 - sqrt(3)/6;

% getting the non-equidistant samples
[q1_c1, q2_c1] =  bandlimited_interpolation(eps_t,...
                                    [q1; q2],c1*eps_t);
[q1_c2, q2_c2] = bandlimited_interpolation(eps_t,...
                                    [q1; q2],c2*eps_t);

%%
% Interlacing the samples:
qeff = zeros(2,2*D);
reff = zeros(2,2*D);

qeff(1,2:2:end)=a2*q1_c1+a1*q1_c2;
qeff(2,2:2:end)=a2*q2_c1+a1*q2_c2;
reff(1,2:2:end)=-kappa*(a2*conj(q1_c1)+a1*conj(q1_c2));
reff(2,2:2:end)=-kappa*(a2*conj(q2_c1)+a1*conj(q2_c2));

qeff(1,1:2:end)=a1*q1_c1+a2*q1_c2;
qeff(2,1:2:end)=a1*q2_c1+a2*q2_c2;
reff(1,1:2:end)=-kappa*(a1*conj(q1_c1)+a2*conj(q1_c2));
reff(2,1:2:end)=-kappa*(a1*conj(q2_c1)+a2*conj(q2_c2));

for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:2*D
        % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_4 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, qeff(1,n), qeff(2,n);...
             reff(1,n), 0, 0;
             reff(2,n), 0, 0];
        eB = expm(B*eps_t);
        eB_1_2 = expm(B*eps_t/2);
        U = (4/3)*eA_4*eB_1_2*eA_4*eA_4*eB_1_2*eA_4+...
            -(1/3)*eA_4*eA_4*eB*eA_4*eA_4;
        S = U*S;
    end
    S=S*z^(2*D*4);
    result_exact(i) = S(1,1);
    result_exact(5+i) = S(1,2);
    result_exact(10+i) = S(1,3);
    result_exact(15+i) = S(2,1);
    result_exact(20+i) = S(2,2);
    result_exact(25+i) = S(2,3);
    result_exact(30+i) = S(3,1);
    result_exact(35+i) = S(3,2);
    result_exact(40+i) = S(3,3);
end

%%
function [q1s, q2s] = bandlimited_interpolation(eps_t, qn, ts)
% implements bandlimited interpolation (interpolation using the time shift
% property of the Fourier Transform)
% 
% inputs:
% eps_t sampling period
% qn    values of the samples of q taken at the points in vector tn
% ts    desired shift of the samples
% 
% output: qs    the shifted vector of samples

Qn = [fft(qn(1,:)); fft(qn(2,:))];
N = length(qn(1,:));
Np = floor(N/2);
Nn = -floor((N-1)/2);
Qn = [Qn(1,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t));...
      Qn(2,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t))];
q1s = ifft(Qn(1,:));
q2s = ifft(Qn(2,:));
end





