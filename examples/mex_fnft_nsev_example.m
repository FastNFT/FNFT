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

% This examples demonstrates how the nonlinear Fourier transform with
% respect to the nonlinear Schroedinger equation with vanishing boundary
% conditions can be computed using mex_fnft_nsev. The signal is a sech
% pulse as discussed in the paper https://doi.org/10.1143/PTPS.55.284

clear all;
close all;

%%% Setup parameters %%%

T = [-10, 10];  % location of the 1st and last sample in the time domain
D = 4096;       % number of samples
XI = [-5, 5];   % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation

%%% Setup the signal %%%

ep_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1):ep_t:T(2);
q = 3.2j*sech(t);               % signal samples

%%% Compute the nonlinear Fourier transform %%%

[contspec, bound_states, residuals] = mex_fnft_nsev(q, T, XI, kappa);
% mex_fnft_nsev has many options => run "help mex_fnft_nsev" to learn more

%%% Plot the results %%%

ep_xi = (XI(2) - XI(1)) / (D - 1);
xi = XI(1):ep_xi:XI(2); % grid in the nonlinear frequency domain

figure;
plot(t, real(q), t, imag(q));
title('Time-domain');
xlabel('t');
ylabel('q(t)');
legend('Real part', 'Imaginary part');

figure;
plot(xi, real(contspec), xi, imag(contspec));
title('Continuous spectrum');
xlabel('\xi');
ylabel('r(\xi)');
legend('Real part', 'Imaginary part');

figure;
plot(real(bound_states), imag(bound_states), 'x');
title('Bound states');
xlabel('Re(\lambda_k)');
ylabel('Im(\lambda_k)');
ylim([0 3]);
axis equal;
