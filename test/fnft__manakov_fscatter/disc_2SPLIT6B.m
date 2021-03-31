%% Matlab file to get the polynomial coefficients for the 4SPLIT4A
%  BO method with 6th order splitting (Eq. 22 in P. J. Prins and S. Wahls,
% Higher order exponential
% splittings for the fast non-linear Fourier transform of the KdV
% equation

clear
close all

syms s1 s2 l kappa h z

AE = [1/z^6,0,0;0,z^6,0;0,0,z^6];
BE = sym('BE',[3,3]);


% we have CE = expm(A+B)
% Let G_approx be approximation of CE after application of splitting
% scheme
% AE_6 = expm(A*h/6)
% then expm(A*h/3) = AE_6*AE_6, expm(A*h/2) = AE_6^3,
% expm(A*h) = AE_12^6
% BE_p = expm(B*h/p)
AE_6 = [1/z,0,0;0,z,0;0,0,z];
BE = sym('BE',[3,3]);
BE_6 = sym('BE_6',[3,3]);
BE_4 = sym('BE_4',[3,3]);
BE_3 = sym('BE_3',[3,3]);
BE_2 = sym('BE_2',[3,3]);

G_approx = (81/40) * BE_6*(AE_6^2*BE_3)^2*AE_6^2*BE_6+...
            -(16/15)*BE_4*AE_6^3*BE_2*AE_6^3*BE_4+...
            (1/24)*BE_2*AE_6^6*BE_2;

% Dividing throughout by z^-6
G_approx_pos = expand(G_approx*z^6);

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
q1 = 1.5*0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 1.5*2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
r1 = -kappa*conj(q1);
r2 = -kappa*conj(q2);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);


%%
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        eA_6 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n);...
             r1(n), 0, 0;
             r2(n), 0, 0];
        eB = expm(B*eps_t);
        eB_2 = expm(B*eps_t/2);
        eB_3 = expm(B*eps_t/3);
        eB_4 = expm(B*eps_t/4);
        eB_6 = expm(B*eps_t/6);
        U = (81/40) * eB_6*(eA_6^2*eB_3)^2*eA_6^2*eB_6+...
            -(16/15)*eB_4*eA_6^3*eB_2*eA_6^3*eB_4+...
            (1/24)*eB_2*eA_6^6*eB_2;
        S = U*S;
    end
    S=S*z^(D*6);
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

%% check p results
p11 = zeros(D*13,1); p12 = zeros(D*13,1); p13 = zeros(D*13,1);
p21 = zeros(D*13,1); p22 = zeros(D*13,1); p23 = zeros(D*13,1);
p31 = zeros(D*13,1); p32 = zeros(D*13,1); p33 = zeros(D*13,1);
 i=1;
for j = D:-1:1

B = [0      q1(j)    q2(j);
     r1(j)   0       0;
     r2(j)   0       0];

e_1_2B = expm(B*eps_t/2);
e_1_3B = expm(B*eps_t/3);
e_1_4B = expm(B*eps_t/4);
e_1_6B = expm(B*eps_t/6);

BE_21_1 = e_1_2B(1,1); BE_21_2 = e_1_2B(1,2); BE_21_3 = e_1_2B(1,3);
BE_22_1 = e_1_2B(2,1); BE_22_2 = e_1_2B(2,2); BE_22_3 = e_1_2B(2,3);
BE_23_1 = e_1_2B(3,1); BE_23_2 = e_1_2B(3,2); BE_23_3 = e_1_2B(3,3);

BE_31_1 = e_1_3B(1,1); BE_31_2 = e_1_3B(1,2); BE_31_3 = e_1_3B(1,3);
BE_32_1 = e_1_3B(2,1); BE_32_2 = e_1_3B(2,2); BE_32_3 = e_1_3B(2,3);
BE_33_1 = e_1_3B(3,1); BE_33_2 = e_1_3B(3,2); BE_33_3 = e_1_3B(3,3);

BE_41_1 = e_1_4B(1,1); BE_41_2 = e_1_4B(1,2); BE_41_3 = e_1_4B(1,3);
BE_42_1 = e_1_4B(2,1); BE_42_2 = e_1_4B(2,2); BE_42_3 = e_1_4B(2,3);
BE_43_1 = e_1_4B(3,1); BE_43_2 = e_1_4B(3,2); BE_43_3 = e_1_4B(3,3);

BE_61_1 = e_1_6B(1,1); BE_61_2 = e_1_6B(1,2); BE_61_3 = e_1_6B(1,3);
BE_62_1 = e_1_6B(2,1); BE_62_2 = e_1_6B(2,2); BE_62_3 = e_1_6B(2,3);
BE_63_1 = e_1_6B(3,1); BE_63_2 = e_1_6B(3,2); BE_63_3 = e_1_6B(3,3);

% coefficients for z^12
p11(13*(i-1)+1) = (BE_21_2*BE_22_1)/24 + (BE_21_3*BE_23_1)/24 - (16*BE_22_2*BE_41_2*BE_42_1)/15 - (16*BE_22_3*BE_41_2*BE_43_1)/15 - (16*BE_23_2*BE_41_3*BE_42_1)/15 - (16*BE_23_3*BE_41_3*BE_43_1)/15 + (81*BE_32_2^2*BE_61_2*BE_62_1)/40 + (81*BE_33_3^2*BE_61_3*BE_63_1)/40 + (81*BE_32_2*BE_32_3*BE_61_2*BE_63_1)/40 + (81*BE_32_2*BE_33_2*BE_61_3*BE_62_1)/40 + (81*BE_32_3*BE_33_2*BE_61_2*BE_62_1)/40 + (81*BE_32_3*BE_33_2*BE_61_3*BE_63_1)/40 + (81*BE_32_3*BE_33_3*BE_61_2*BE_63_1)/40 + (81*BE_33_2*BE_33_3*BE_61_3*BE_62_1)/40;
p12(13*(i-1)+1) = (BE_21_2*BE_22_2)/24 + (BE_21_3*BE_23_2)/24 - (16*BE_22_2*BE_41_2*BE_42_2)/15 - (16*BE_22_3*BE_41_2*BE_43_2)/15 - (16*BE_23_2*BE_41_3*BE_42_2)/15 - (16*BE_23_3*BE_41_3*BE_43_2)/15 + (81*BE_32_2^2*BE_61_2*BE_62_2)/40 + (81*BE_33_3^2*BE_61_3*BE_63_2)/40 + (81*BE_32_2*BE_32_3*BE_61_2*BE_63_2)/40 + (81*BE_32_2*BE_33_2*BE_61_3*BE_62_2)/40 + (81*BE_32_3*BE_33_2*BE_61_2*BE_62_2)/40 + (81*BE_32_3*BE_33_2*BE_61_3*BE_63_2)/40 + (81*BE_32_3*BE_33_3*BE_61_2*BE_63_2)/40 + (81*BE_33_2*BE_33_3*BE_61_3*BE_62_2)/40;
p13(13*(i-1)+1) = (BE_21_2*BE_22_3)/24 + (BE_21_3*BE_23_3)/24 - (16*BE_22_2*BE_41_2*BE_42_3)/15 - (16*BE_22_3*BE_41_2*BE_43_3)/15 - (16*BE_23_2*BE_41_3*BE_42_3)/15 - (16*BE_23_3*BE_41_3*BE_43_3)/15 + (81*BE_32_2^2*BE_61_2*BE_62_3)/40 + (81*BE_33_3^2*BE_61_3*BE_63_3)/40 + (81*BE_32_2*BE_32_3*BE_61_2*BE_63_3)/40 + (81*BE_32_2*BE_33_2*BE_61_3*BE_62_3)/40 + (81*BE_32_3*BE_33_2*BE_61_2*BE_62_3)/40 + (81*BE_32_3*BE_33_2*BE_61_3*BE_63_3)/40 + (81*BE_32_3*BE_33_3*BE_61_2*BE_63_3)/40 + (81*BE_33_2*BE_33_3*BE_61_3*BE_62_3)/40;
p21(13*(i-1)+1) = (BE_22_1*BE_22_2)/24 + (BE_22_3*BE_23_1)/24 - (16*BE_22_2*BE_42_1*BE_42_2)/15 - (16*BE_22_3*BE_42_2*BE_43_1)/15 - (16*BE_23_2*BE_42_1*BE_42_3)/15 - (16*BE_23_3*BE_42_3*BE_43_1)/15 + (81*BE_32_2^2*BE_62_1*BE_62_2)/40 + (81*BE_33_3^2*BE_62_3*BE_63_1)/40 + (81*BE_32_2*BE_32_3*BE_62_2*BE_63_1)/40 + (81*BE_32_2*BE_33_2*BE_62_1*BE_62_3)/40 + (81*BE_32_3*BE_33_2*BE_62_1*BE_62_2)/40 + (81*BE_32_3*BE_33_2*BE_62_3*BE_63_1)/40 + (81*BE_32_3*BE_33_3*BE_62_2*BE_63_1)/40 + (81*BE_33_2*BE_33_3*BE_62_1*BE_62_3)/40;
p22(13*(i-1)+1) = BE_22_2^2/24 - (16*BE_22_2*BE_42_2^2)/15 + (81*BE_32_2^2*BE_62_2^2)/40 + (BE_22_3*BE_23_2)/24 - (16*BE_22_3*BE_42_2*BE_43_2)/15 - (16*BE_23_2*BE_42_2*BE_42_3)/15 - (16*BE_23_3*BE_42_3*BE_43_2)/15 + (81*BE_32_3*BE_33_2*BE_62_2^2)/40 + (81*BE_33_3^2*BE_62_3*BE_63_2)/40 + (81*BE_32_2*BE_32_3*BE_62_2*BE_63_2)/40 + (81*BE_32_2*BE_33_2*BE_62_2*BE_62_3)/40 + (81*BE_32_3*BE_33_2*BE_62_3*BE_63_2)/40 + (81*BE_32_3*BE_33_3*BE_62_2*BE_63_2)/40 + (81*BE_33_2*BE_33_3*BE_62_2*BE_62_3)/40;
p23(13*(i-1)+1) = (BE_22_2*BE_22_3)/24 - (16*BE_23_2*BE_42_3^2)/15 + (BE_22_3*BE_23_3)/24 - (16*BE_22_2*BE_42_2*BE_42_3)/15 - (16*BE_22_3*BE_42_2*BE_43_3)/15 - (16*BE_23_3*BE_42_3*BE_43_3)/15 + (81*BE_32_2*BE_33_2*BE_62_3^2)/40 + (81*BE_33_2*BE_33_3*BE_62_3^2)/40 + (81*BE_32_2^2*BE_62_2*BE_62_3)/40 + (81*BE_33_3^2*BE_62_3*BE_63_3)/40 + (81*BE_32_2*BE_32_3*BE_62_2*BE_63_3)/40 + (81*BE_32_3*BE_33_2*BE_62_2*BE_62_3)/40 + (81*BE_32_3*BE_33_2*BE_62_3*BE_63_3)/40 + (81*BE_32_3*BE_33_3*BE_62_2*BE_63_3)/40;
p31(13*(i-1)+1) = (BE_22_1*BE_23_2)/24 + (BE_23_1*BE_23_3)/24 - (16*BE_22_2*BE_42_1*BE_43_2)/15 - (16*BE_22_3*BE_43_1*BE_43_2)/15 - (16*BE_23_2*BE_42_1*BE_43_3)/15 - (16*BE_23_3*BE_43_1*BE_43_3)/15 + (81*BE_32_2^2*BE_62_1*BE_63_2)/40 + (81*BE_33_3^2*BE_63_1*BE_63_3)/40 + (81*BE_32_2*BE_32_3*BE_63_1*BE_63_2)/40 + (81*BE_32_2*BE_33_2*BE_62_1*BE_63_3)/40 + (81*BE_32_3*BE_33_2*BE_62_1*BE_63_2)/40 + (81*BE_32_3*BE_33_2*BE_63_1*BE_63_3)/40 + (81*BE_32_3*BE_33_3*BE_63_1*BE_63_2)/40 + (81*BE_33_2*BE_33_3*BE_62_1*BE_63_3)/40;
p32(13*(i-1)+1) = (BE_22_2*BE_23_2)/24 - (16*BE_22_3*BE_43_2^2)/15 + (BE_23_2*BE_23_3)/24 - (16*BE_22_2*BE_42_2*BE_43_2)/15 - (16*BE_23_2*BE_42_2*BE_43_3)/15 - (16*BE_23_3*BE_43_2*BE_43_3)/15 + (81*BE_32_2*BE_32_3*BE_63_2^2)/40 + (81*BE_32_3*BE_33_3*BE_63_2^2)/40 + (81*BE_32_2^2*BE_62_2*BE_63_2)/40 + (81*BE_33_3^2*BE_63_2*BE_63_3)/40 + (81*BE_32_2*BE_33_2*BE_62_2*BE_63_3)/40 + (81*BE_32_3*BE_33_2*BE_62_2*BE_63_2)/40 + (81*BE_32_3*BE_33_2*BE_63_2*BE_63_3)/40 + (81*BE_33_2*BE_33_3*BE_62_2*BE_63_3)/40;
p33(13*(i-1)+1) = BE_23_3^2/24 - (16*BE_23_3*BE_43_3^2)/15 + (81*BE_33_3^2*BE_63_3^2)/40 + (BE_22_3*BE_23_2)/24 - (16*BE_22_2*BE_42_3*BE_43_2)/15 - (16*BE_22_3*BE_43_2*BE_43_3)/15 - (16*BE_23_2*BE_42_3*BE_43_3)/15 + (81*BE_32_3*BE_33_2*BE_63_3^2)/40 + (81*BE_32_2^2*BE_62_3*BE_63_2)/40 + (81*BE_32_2*BE_32_3*BE_63_2*BE_63_3)/40 + (81*BE_32_2*BE_33_2*BE_62_3*BE_63_3)/40 + (81*BE_32_3*BE_33_2*BE_62_3*BE_63_2)/40 + (81*BE_32_3*BE_33_3*BE_63_2*BE_63_3)/40 + (81*BE_33_2*BE_33_3*BE_62_3*BE_63_3)/40;

% coefs for z^8
p31(13*(i-1)+5) = (81*BE_31_2*BE_32_3*BE_63_1^2)/40 + (81*BE_31_3*BE_33_3*BE_63_1^2)/40 + (81*BE_31_2*BE_32_1*BE_62_1*BE_63_2)/40 + (81*BE_31_2*BE_32_2*BE_62_1*BE_63_1)/40 + (81*BE_32_1*BE_32_2*BE_61_1*BE_63_2)/40 + (81*BE_31_2*BE_33_1*BE_62_1*BE_63_3)/40 + (81*BE_31_3*BE_32_1*BE_63_1*BE_63_2)/40 + (81*BE_31_3*BE_33_2*BE_62_1*BE_63_1)/40 + (81*BE_32_1*BE_33_2*BE_61_1*BE_63_3)/40 + (81*BE_32_3*BE_33_1*BE_61_1*BE_63_2)/40 + (81*BE_31_3*BE_33_1*BE_63_1*BE_63_3)/40 + (81*BE_33_1*BE_33_3*BE_61_1*BE_63_3)/40;

% coefs for z^6

% coefs for z^4

i=i+1;
end

p = [p11; p12; p13; p21; p22; p23; p31; p32; p33];
