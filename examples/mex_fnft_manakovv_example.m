% This file is part of FNFT.  
%                                                                  
% FNFT is free software; you can redistribute it and/or
% modify it under the terms of the version 2 of the GNU General
% Public License as published by the Free Software Foundation.
%
% FNFT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%                                                                      
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Contributors:
% Sander Wahls (TU Delft) 2017-2018.
% Lianne de Vries (TU Delft student) 2021.

% This examples demonstrates how the nonlinear Fourier transform with
% respect to the nonlinear Schroedinger equation with vanishing boundary
% conditions can be computed using mex_fnft_nsev. The signal is a sech
% pulse as discussed in the paper https://doi.org/10.1143/PTPS.55.284

clear all;
close all;

%%% Setup parameters %%%

T = [-5, 5];  % location of the 1st and last sample in the time domain
D = 512;       % number of samples
XI = [-7/4, 8/4];   % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation

%%% Setup the signal %%%

eps_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1):eps_t:T(2);
q1 = 0.8*sech(t);               % signal samples
q2 = 5.2*sech(t);
% q1 = sqrt(2)*sech(t);               % signal samples
% q2 = sqrt(2)*sech(t);

% xi = 4.869;     % parameter defining the velocity of the soliton. real number
% eta = 0.56;     % parameter defining the amplitude of the soliton. real number
% x = 0.1;          % x coordinate at which we get the potential q(x,t) = q(t)
% S = [6, 1+5i];  % vector defining the polarization of both components of q
% c = S/norm(S);      % normaliza vector c for polarization
%     x0 = log(norm(S)^2)/(4*eta);
%     q1 = conj(-2*eta*sech(2*eta*(t-x0)+8*xi*eta*x).*exp(2*1i*xi*t+4i*(xi^2-eta^2)*x)*c(1));
%     q2 = conj(-2*eta*sech(2*eta*(t-x0)+8*xi*eta*x).*exp(2*1i*xi*t+4i*(xi^2-eta^2)*x)*c(2));

%%% Compute the nonlinear Fourier transform %%%

[contspec, bound_states] = mex_fnft_manakovv(complex(q1), complex(q2), T, XI, kappa);
% mex_fnft_manakovv has many options => run "help mex_fnft_manakovv" to learn more

%%% Plot the results %%%

ep_xi = (XI(2) - XI(1)) / (D - 1);
xi = XI(1):ep_xi:XI(2); % grid in the nonlinear frequency domain

figure;
plot(t, real(q1), t, imag(q1));
title('Time-domain');
xlabel('t');
ylabel('q1(t)');
legend('Real part', 'Imaginary part');

figure;
plot(t, real(q2), t, imag(q2));
title('Time-domain');
xlabel('t');
ylabel('q2(t)');
legend('Real part', 'Imaginary part');

% Note: plotting the results for the default case where 
% contspec = [rho1; rho2]. If different inputs are chosen for cstype or M
% plots should be adjusted
figure;
plot(xi, real(contspec(1:length(contspec)/2)), xi, imag(contspec(1:length(contspec)/2)));
title('Continuous spectrum, first element');
xlabel('\xi');
ylabel('r(\xi)');
legend('Real part', 'Imaginary part');

figure;
plot(xi, real(contspec(length(contspec)/2+1:length(contspec))), xi, imag(contspec(length(contspec)/2+1:length(contspec))));
title('Continuous spectrum, second element');
xlabel('\xi');
ylabel('r(\xi)');
legend('Real part', 'Imaginary part');
