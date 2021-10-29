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
% Sander Wahls (TU Delft) 2017-2018, 2021.
% Lianne de Vries (TU Delft student) 2021.

% This examples demonstrates how the nonlinear Fourier transform with
% respect to the Manakov equation with vanishing boundary
% conditions can be computed using mex_fnft_manakovv. The signal is a sech
% pulse as discussed in 
% http://resolver.tudelft.nl/uuid:0276e693-3408-4472-9749-b754c2114183


clear all;
close all;

%%% Setup parameters %%%

T = [-7, 7];    % location of the 1st and last sample in the time domain
D = 512;        % number of samples
XI = [-7/4, 8/4]; % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation

%%% Setup the signal %%%

eps_t = (T(2) - T(1)) / (D - 1);    % time domain step size
t = T(1):eps_t:T(2);
q1 = 0.8*sech(t);                   % signal samples
q2 = 5.2*sech(t);

%%% Compute the nonlinear Fourier transform %%%

[contspec, bound_states] = mex_fnft_manakovv(complex(q1), complex(q2), ...
    T, XI, kappa);
% mex_fnft_manakovv has many options => run "help mex_fnft_manakovv" to
% learn more

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
plot(xi, real(contspec(1:length(contspec)/2)), xi, ...
    imag(contspec(1:length(contspec)/2)));
title('Continuous spectrum, first element');
xlabel('\xi');
ylabel('r(\xi)');
legend('Real part', 'Imaginary part');

figure;
plot(xi, real(contspec(length(contspec)/2+1:length(contspec))), xi, ...
    imag(contspec(length(contspec)/2+1:length(contspec))));
title('Continuous spectrum, second element');
xlabel('\xi');
ylabel('r(\xi)');
legend('Real part', 'Imaginary part');
