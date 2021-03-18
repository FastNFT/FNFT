%% Matlab file to get the polynomial coefficients for the 2SPLIT3B
% also to generate results for the test file
% BO (2nd order base method) with 3rd order splitting (Eq. 19)

%% computing the coefficients in sybmbolic form
syms q1 q2 l kappa h z

A = [-1i*l,0,0;0,1i*l,0;0,0,1i*l];
B = [0,q1,q2;-kappa*conj(q1),0,0;-kappa*conj(q2),0,0];

% Let AE = expm(A*h) and BE = expm(B*h)
% Set z = exp(j lambda h/3) (m=1/3) such that
AE = [1/z^3,0,0;0,z^3,0;0,0,z^3];
BE = sym('BE',[3,3]);
% AE_3 = expm(A*h/3)
AE_3 = [1/z,0,0;0,z,0;0,0,z];
% BE_1_3 = expm(B*h/3)
BE_1_3 = sym('BE_1_3',[3,3]);
% BE_2_3 = expm(B*h*2/3)
BE_2_3 = sym('BE_2_3',[3,3]);


% CE=expm((A+B)*h)
% Let CE_approx be approximation of CE after application of splitting
% scheme

CE_approx = (9/8)*BE_1_3*AE_3*AE_3*BE_2_3*AE_3 -...
            (1/8)*BE*AE;
CE_approx_pos = expand(CE_approx*z^3);

% We can now look at coefficients of individual polynomials. Here c gives
% the coefficients of the element specified by (i,j) in CE_approx_pos(i,j)
% and t gives which power of z the coefficients belong to
[c11,t11] = coeffs(CE_approx_pos(1,1),z);
[c12,t12] = coeffs(CE_approx_pos(1,2),z);
[c13,t13] = coeffs(CE_approx_pos(1,3),z);

[c21,t21] = coeffs(CE_approx_pos(2,1),z);
[c22,t22] = coeffs(CE_approx_pos(2,2),z);
[c23,t23] = coeffs(CE_approx_pos(2,3),z);

[c31,t31] = coeffs(CE_approx_pos(3,1),z);
[c32,t32] = coeffs(CE_approx_pos(3,2),z);
[c33,t33] = coeffs(CE_approx_pos(3,3),z);

% Organizing all coefficients in corresponding matrices
matz0 = [c11(2),    0,     0;
         c21(2),    0,     0;
         c31(2),    0,     0];
         
    
matz2 = [0,         c12(2),  c13(2);
         0,         c22(2),  c23(2);
         0,         c32(2),  c33(2)];

matz4 = [c11(1),    0,  0;
         c21(1),    0,  0;
         c31(1),    0,  0];

matz6 = [0,         c12(1),  c13(1);
         0,         c22(1),  c23(1);
         0,         c32(1),  c33(1)];
     
     
%% Computing values for the test file
eps_t = 0.13;
kappa = 1;
D=8;
q = (0.41*cos(1:2*D)+0.59j*sin(0.28*(1:2*D)))*50;
r = -conj(q);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 19 in P. J. Prins and S. Wahls, Higher order exponential
        % splittings for the fast non-linear Fourier transform of the KdV
        % equation, to appear in Proc. of IEEE ICASSP 2018, Calgary, April 2018.
        eA_3 = [1/z 0 0; 0 z 0; 0 0 z];
        B = [0, q(2*n-1), q(2*n); r(2*n-1), 0, 0; r(2*n), 0, 0];
        U = (9/8)*expm(B*eps_t/3)*eA_3^2*expm(B*eps_t*2/3)*eA_3 - ...
            (1/8)*expm(B*eps_t)*eA_3^3;
        S = U*S;
    end
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


