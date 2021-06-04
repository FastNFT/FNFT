%% This function determines the exact solution of the Manakov equation
% with rectangle potential:
% q1 = A1 for L1 < t < L2
% q2 = A2 for L1 < t < L2
%% Inputs:
% A1            amplitude of first entry of potential
% A2            amplitude of second entry of potential
% XI_vector     vector cotaining all xi's we want to evaluate
% L             L = [L1 L2], where L1 is the left boundary of the rectangle
%               and L2 the right boundary
% kappa         -1 defocussing, +1 focussing
function [a, b1, b2] = Manakov_rectangle_exact(A1, A2, XI_vector, kappa, L)
a=zeros(1,length(XI_vector));
b1=zeros(1,length(XI_vector));
b2=zeros(1,length(XI_vector));
for i=1:length(XI_vector)
    L1 = L(1); L2 = L(2);
    lam = XI_vector(i);
    P0 = [-1j*lam, 0, 0;
          0, 1j*lam, 0;
          0, 0, 1j*lam];
    PL = [-1j*lam, A1, A2;
          -kappa*conj(A1), 1j*lam, 0;
          -kappa*conj(A2), 0, 1j*lam];
    
      % For the first part -inf<t<L1, we have v = [exp(j*lam*t); 0; 0]
      % From this we get the BC at L1:
      v_L1 = [exp(-1j*lam*L1); 0; 0];
      
      % L1<t<L2, v=exp(PL*t)*c2 where we get the constant c2 from v_L1
      c2 = expm(PL*L1)\v_L1;
      v_L2 = expm(PL*L2)*c2;
      
      % L2<t, v = exp(P0*t)*c2 where we get c3 from v_L2
      c3 = expm(P0*L2)\v_L2;
%      v_fun = @(t)expm(PL*t)*c3;
%       v_end = expm(P0*t)*c3;
%       ab(1) = v_end(1)*exp(1j*lam*t);
%       ab(2) = v_end(2)*exp(-1j*lam*t);
%       ab(3) = v_end(3)*exp(-1j*lam*t);
% expm(P0*t)*[epx(j*lam*t); epx(-j*lam*t); epx(-j*lam*t)] cancel against
% each other (becomes [1; 1; 1], so we do not need ot specify at which t
% and we just get c3 as the answer
    a(i) = c3(1);
    b1(i) = c3(2);
    b2(i) = c3(3);
end
end

