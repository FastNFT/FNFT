%% Matlab file to get the polynomial coefficients for the 4SPLIT4A
% BO method with 6th order splitting (Eq. 22 in P. J. Prins and S. Wahls,
% Higher order exponential
% splittings for the fast non-linear Fourier transform of the KdV
% equation

clear
close all

syms s1 s2 l kappa h z

AE = [1/z^12,0,0;0,z^12,0;0,0,z^12];
BE = sym('BE',[3,3]);


% we have CE = expm(A+B)
% Let G_approx be approximation of CE after application of splitting
% scheme
% AE_12 = expm(A*h/12)
% then expm(A*h/6) = AE_12*AE_12, expm(A*h/4) = AE_12*AE_12*AE_12,
% expm(A*h/3) = AE_12^4, expm(A*h/2) = AE_12^6
% BE_p = expm(B*h/p)
AE_12 = [1/z,0,0;0,z,0;0,0,z];
BE = sym('BE',[3,3]);
BE_3 = sym('BE_3',[3,3]);
BE_2 = sym('BE_2',[3,3]);


G_approx = (81/40) * AE_12^2*(BE_3*AE_12^4)^2*BE_3*AE_12^2+...
            -(16/15)*AE_12^3*BE_2*AE_12^6*BE_2*AE_12^3+...
            (1/24)*AE_12^6*BE*AE_12^6;

% Dividing throughout by z^-12
G_approx_pos = expand(G_approx*z^12);

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
% TODO: update
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        eA_12 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB = expm(B*eps_t);
        eB_2 = expm(B*eps_t/2);
        eB_3 = expm(B*eps_t/3);
        U = (81/40) * eA_12^2*(eB_3*eA_12^4)^2*eB_3*eA_12^2+...
            -(16/15)*eA_12^3*eB_2*eA_12^6*eB_2*eA_12^3+...
            (1/24)*eA_12^6*eB*eA_12^6;
        S = U*S;
    end
    S=S*z^(D*12);
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
