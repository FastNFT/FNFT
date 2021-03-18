%% Matlab file to get the polynomial coefficients for the FCF24_4
% CF24 method with 4th order splitting

syms s1 s2 l kappa h z

A = [-1i*l,0,0;0,1i*l,0;0,0,1i*l];
B = [0,s1,s2;-kappa*conj(s1),0,0;-kappa*conj(s2),0,0];
% s1, s2 are some combination of the q1, q2's. See overleaf
% Let AE = expm(A*h) and BE = expm(B*h)
% because we are only splitting one of the 2 matrix exponentials used in
% this method, this splitting is now basically the same as the splitting
% for TES4
% Set z = exp(j*xi*1/4(a1+a2)*h) = exp(j*xi*1/8*h) such that
AE = [1/z^8,0,0;0,z^8,0;0,0,z^8];
BE = sym('BE',[3,3]);


% we have CE = expm(A+B)
% Let G_approx be approximation of CE after application of splitting
% scheme
% AE_8 = expm(A*h/8)
AE_8 = [1/z,0,0;0,z,0;0,0,z];
BE_4 = sym('BE_4',[3,3]);


G_approx = ((4/3)*AE_8*(BE_4*BE_4)*(AE_8*AE_8)*(BE_4*BE_4)*AE_8-...
        (1/3)*(AE_8*AE_8)*BE*(AE_8*AE_8))

% Dividing throughout by z^-4
G_approx_pos = expand(G_approx*z^4)

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
    