%% Matlab file to get the polynomial coefficients for the 4SPLIT4A
%  BO method with 4th order splitting
%clear
close all

syms s1 s2 l kappa h z

AE = [1/z^2,0,0;0,z^2 0;0,0,z^2];
BE = sym('BE',[3,3]);


% we have CE = expm(A+B)
% Let G_approx be approximation of CE after application of splitting
% scheme, choose z = exp(j*lambda*h/2)
% AE_2 = expm(A*h/2)
% BE_2 = expm(2B*h/4) = expm(B*h/2)
% BE_4 = expm(B*h/4)
AE_2 = [1/z,0,0;0,z,0;0,0,z];
BE_2 = sym('BE_2',[3,3]);
BE_4 = sym('BE_4',[3,3]);

G_approx = ((4/3)*BE_4*AE_2*BE_2*AE_2*BE_4)+...
            -((1/3)*BE_2*AE_2*AE_2*BE_2);

% Dividing throughout by z^-4
G_approx_pos = expand(G_approx*z^2);

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
    
%% Exact answer for test file:
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
r1 = -kappa*conj(q1);
r2 = -kappa*conj(q2);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);


%%
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 20 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_2 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB_1_2 = expm(B*eps_t/2);
        eB_1_4 = expm(B*eps_t/4);
        U = ((4/3)*eB_1_4*eA_2*eB_1_2*eA_2*eB_1_4)+...
            -((1/3)*eB_1_2*eA_2*eA_2*eB_1_2);
        S = U*S;
    end
    S=S*z^(D*2);
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
