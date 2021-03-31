%% Matlab file to get the polynomial coefficients for the 2SPLIT3S
% also to generate results for the test file
% BO (2nd order base method) with 3rd order splitting
% symmethic weighted sequential splitting (Strang splitting)

%% computing the coefficients in sybmbolic form
syms q1 q2 l kappa h z

A = [-1i*l,0,0;0,1i*l,0;0,0,1i*l];
B = [0,q1,q2;-kappa*conj(q1),0,0;-kappa*conj(q2),0,0];

% Let AE = expm(A*h) and BE = expm(B*h)
% Set z = exp(j lambda h/2) (m=1/2) such that
AE = [1/z^2,0,0;0,z^2,0;0,0,z^2];
BE = sym('BE',[3,3]);
% AE_2 = expm(A*h/2)
AE_2 = [1/z,0,0;0,z,0;0,0,z];
% BE_2 = expm(B*h/2)
BE_2 = sym('BE_2',[3,3]);

% CE=expm((A+B)*h)
% Let CE_approx be approximation of CE after application of splitting
% scheme

CE_approx = (4/3)*(AE_2*BE*AE_2 + BE_2*AE*BE_2)/2-...
            (1/3)*(AE*BE + BE*AE)/2;
CE_approx_pos = expand(CE_approx*z^2);

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

%% Computing values for the test file
eps_t = 0.13;
kappa = 1;
D=8;
q = (0.41*cos(1:2*D)+0.59j*sin(0.28*(1:2*D)))*50;
r = -kappa*conj(q);
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        % Eq. 4.14 in Third Order Difference Methods for Hyperbolic
        % Equations S.Z. Burstein, A.A. Mirin
        % (https://doi.org/10.1016/0021-9991(70)90080-X)
        eA_2 = [z^-1 0 0; 0 z 0; 0 0 z];
        eA = [z^(-2) 0 0; 0 z^2 0; 0 0 z];
        B = [0, q(2*n-1), q(2*n); r(2*n-1), 0, 0; r(2*n), 0, 0];
        U = (4/3)*(eA_2*expm(B*eps_t)*eA_2 + expm(B*eps_t/2)*eA*expm(B*eps_t/2))/2-...
            (1/3)*(eA*expm(B*eps_t) + expm(B*eps_t)*eA)/2;
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
