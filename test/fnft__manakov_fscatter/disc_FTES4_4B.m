%% Matlab file to get the polynomial coefficients for FTES4
% also to generate results for the test file

syms q1 q2 l kappa h z...

A = [-1i*l,0,0;0,1i*l,0;0,0,1i*l];
B = [0,q1,q2;-kappa*conj(q1),0,0;-kappa*conj(q2),0,0];
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
BE_2 = sym('BE_2',[3,3]);
BE_4 = sym('BE_4',[3,3]);
E1 = sym('E1',[3,3]);
E2 = sym('E2',[3,3]);

CE_approx = E1*((4/3)*BE_4*(AE_4*AE_4)*(BE_2)*(AE_4*AE_4)*BE_4-...
        (1/3)*BE_2*AE_4*AE_4*AE_4*AE_4*BE_2)*E2;
    
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

%% Computing values for the test file
eps_t = 0.13;
kappa = 1;
D=512;
q1 = 0.92*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
q2 = 2.13*sech((-(D-1)*eps_t/2):eps_t:((D-1)*eps_t/2));
zlist = exp(1i*[0,pi/4,9*pi/14,4*pi/3,-pi/5]);
result_exact = zeros(45,1);
for i = 1:1:5
    z = zlist(i);
    S = eye(3);
    for n=1:D
        if n == 1
    Q1 = [0                                 q1(2)-q1(1) q2(2)-q2(1);
          -kappa*(conj(q1(2))-conj(q1(1)))  0           0;
          -kappa*(conj(q2(2))-conj(q2(1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(2)-2*q1(1)+q1(1)       q2(2)-2*q2(1)+q2(1);
          -kappa*(conj(q1(2))-2*conj(q1(1))+conj(q1(1)))    0 0;
          -kappa*(conj(q2(2))-2*conj(q2(1))+conj(q2(1)))    0 0]/(eps_t^2);
        elseif n == D
    Q1 = [0                                 q1(n)-q1(n-1) q2(n)-q2(n-1);
          -kappa*(conj(q1(n))-conj(q1(n-1)))  0           0;
          -kappa*(conj(q2(n))-conj(q2(n-1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(n)-2*q1(n)+q1(n-1)       q2(n)-2*q2(n)+q2(n-1);
          -kappa*(conj(q1(n))-2*conj(q1(n))+conj(q1(n-1)))    0 0;
          -kappa*(conj(q2(n))-2*conj(q2(n))+conj(q2(n-1)))    0 0]/(eps_t^2);
        else
    Q1 = [0                                 q1(n+1)-q1(n-1) q2(n+1)-q2(n-1);
          -kappa*(conj(q1(n+1))-conj(q1(n-1)))  0           0;
          -kappa*(conj(q2(n+1))-conj(q2(n-1)))  0           0]/(2*eps_t);
    Q2 = [0         q1(n+1)-2*q1(n)+q1(n-1)       q2(n+1)-2*q2(n)+q2(n-1);
          -kappa*(conj(q1(n+1))-2*conj(q1(n))+conj(q1(n-1)))    0 0;
          -kappa*(conj(q2(n+1))-2*conj(q2(n))+conj(q2(n-1)))    0 0]/(eps_t^2);
        end
  
E1 = expm(eps_t^2*Q1/12 + eps_t^3*Q2/48);
E2 = expm(-eps_t^2*Q1/12 + eps_t^3*Q2/48);
        eA_4 = [z^-1 0 0; 0 z 0; 0 0 z];
        B = [0, q1(n), q2(n); -kappa*conj(q1(n)), 0, 0; -kappa*conj(q2(n)), 0, 0];
        eB_1_2 = expm(B*eps_t/2);
        eB_1_4 = expm(B*eps_t/4);
        U = E1*((4/3)*eB_1_4*eA_4*eA_4*eB_1_2*eA_4*eA_4*eB_1_4+...
            -(1/3)*eB_1_2*eA_4*eA_4*eA_4*eA_4*eB_1_2)*E2;
        S = U*S;
    end
    S=S*z^(D*4);
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


