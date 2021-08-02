%% This script determines symbolic expression for the exact solution of 
% the Manakov equation with rectangle potential
% q1 = A1 for L1 < t < L2
% q2 = A2 for L1 < t < L2
%% Parameters:
% A1            amplitude of first entry of potential
% A2            amplitude of second entry of potential
% L             L = [L1 L2], where L1 is the left boundary of the rectangle
%               and L2 the right boundary
% kappa         -1 defocussing, +1 focussing

syms lam t kappa A1 A2 L1 L2
P0 = [-1j*lam, 0, 0;
          0, 1j*lam, 0;
          0, 0, 1j*lam];
PL = [-1j*lam, A1, A2;
      -kappa*conj(A1), 1j*lam, 0;
      -kappa*conj(A2), 0, 1j*lam];
  
vL1 = [exp(-1j*lam*L1);  0;  0];     % value of v at L1
vt = expm(P0*t)*inv(expm(P0*L2))*expm(PL*L2)*inv(expm(PL*L1))*vL1;
vt_simp = simplify(vt);

a = vt(1)*exp(1j*lam*t);
b1 = vt(2)*exp(-1j*lam*t);
b2 = vt(3)*exp(-1j*lam*t);

