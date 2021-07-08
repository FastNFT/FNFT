%% This function determines the exact solution of the Manakov equation
% with rectangle potential:
% q1 = A1 for L1 < t < L2
% q2 = A2 for L1 < t < L2
%% Inputs:
% S             amplitude of first entry of potential
% xi            parameter defining the velocity of the soliton
% eta           parameter defining the amplitude of the soliton
% x             This is the x coordinate at which we get the potential q(x,t) = q(t)
% lam_vector    points in the frequency domain where we determined the a
%               coefficient
function [a] = Manakov_soliton_exact(xi, eta, lam_vector)
zeta = xi+1i*eta;
a = ones(size(lam_vector));
a = a.*((lam_vector - zeta)./(lam_vector - conj(zeta)));
end

