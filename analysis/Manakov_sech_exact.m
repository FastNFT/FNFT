%% This function determines the exact solution of the Manakov equation
% equipped with secant hyperbolic potential
%% Inputs:
% A1            amplitude of first entry of potential: q1 = A1 sech(t)
% A2            amplitude of first entry of potential: q2 = A2 sech(t)
% XI_vector     vector containing all xi's we want to evaluate
% t_vector      vector containing all time instances where we take a sample
% kappa         -1 defocussing, +1 focussing
function [a, b1, b2] = Manakov_sech_exact(A1, A2, XI_vector, kappa)
qo = sqrt(kappa)*sqrt(abs(A1)^2+abs(A2)^2);
%qo = sqrt(kappa)*sqrt(A1^2+A2^2);
a_fun = @(xi) double((gamma(sym(0.5-1i*(xi)))).^2./(gamma(sym(0.5-1i*(xi)+qo)).*gamma(sym(0.5-1i*(xi)-qo))));
a = a_fun(XI_vector);
b1 = -kappa*(sin(pi*qo)*A1/qo)*sech(XI_vector*pi);
b2 = -kappa*(sin(pi*qo)*A2/qo)*sech(XI_vector*pi);
end

