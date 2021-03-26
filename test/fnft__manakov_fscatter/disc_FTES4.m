%% Matlab file to get the polynomial coefficients for FTES4
% also to generate results for the test file

%% computing the coefficients in sybmbolic form
syms q1 q2 l kappa h z

% Let AE = expm(A*h) and BE = expm(B*h)
% Set z = exp(j lambda h/4) such that
AE = [1/z^4,0,0;0,z^4,0;0,0,z^4];
BE = sym('BE',[3,3]);
% CE=expm((A+B)*h)
% Let CE_approx be approximation of CE after application of splitting
% scheme
% AE_4 = expm(A*h/4)
AE_4 = [1/z,0,0;0,z,0;0,0,z];
% BE_4 = expm(B*h/4)
BE_4 = sym('BE_4',[3,3]);

CE_approx = (4/3)*AE_4*(BE_4*BE_4)*(AE_4*AE_4)*(BE_4*BE_4)*AE_4-...
        (1/3)*(AE_4*AE_4)*BE*(AE_4*AE_4);

% Dividing throughout by z^-4
CE_approx_pos = expand(CE_approx*z^4);

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

matz0 = [c11(2),    c12(2),     c13(2);
         0,         0,          0;
         0,         0,          0];
         
    
matz2 = [0,         0,          0;
         c21(2),    c22(2),     c23(2);
         c31(2),    c32(2),     c33(2)];

matz4 = [c11(1),    c12(1),     c13(1);
         0,         0,          0;
         0,         0,          0];

matz6 = [0,         0,          0;
         c21(1),    c22(1),     c23(1);
         c31(1),    c32(1),     c33(1)];

%% Computing values for the test file
% TODO: update
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
        eA_3 = [z^-1 0 0; 0 z 0; 0 0 z];
        B = [0, q(2*n-1), q(2*n); r(2*n-1), 0, 0; r(2*n), 0, 0];
        U = (9/8)*eA_3*expm(B*eps_t*2/3)*eA_3^2*expm(B*eps_t/3) - (1/8)*eA_3^3*expm(B*eps_t);
        S = U*S;
    end
    S=S*z^(D*3);
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


