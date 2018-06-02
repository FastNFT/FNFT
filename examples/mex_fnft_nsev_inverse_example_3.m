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
% Sander Wahls (TU Delft) 2018.

% The fast inverse nonlinear Fourier transform routine
% mex_finft_nsev_inverse supports different methods to invert the
% continuous spectrum. The default method is usually fine, but sometimes
% an alternative method can give better results. This example illustrates
% such a situation in the focusing case.

clear all;
close all;

%%% Setup parameters %%%

T = [-20 20];   % Location of the 1st and last sample in the time domain
D = 2^10;       % Desired number of samples in the time domain
M = D;          % Number of samples in the nonlinear frequency domain.
                % Increasing M does not help in this example. In fact, it
                % deterioates the performance of the alternative method
                % to that of the default method.
XI = mex_fnft_nsev_inverse_XI(D, T, M);
                % Location of the 1st and last sample in the nonlinear
                % frequency domain -> we currently have to use the specific
                % values returned by mex_fnft_nsev_inverse_XI
XI_plot = [-4, 4];
                % Nonlinear frequency range used for plotting and error
                % computation
kappa = +1;     % Docusing nonlinear Schroedinger equation

%%% Define the desired continuous spectrum and compute the corresponding
% time domain signal numerically using two different methods %%%

eps_xi = (XI(2) - XI(1))/(M - 1);
xi = XI(1) + (0:M-1)*eps_xi;
contspec_exact = complex(200*exp(-xi.^4));

q_via_default_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa);

q_via_alternative_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa, 'csinv_b_from_a');

%%% Compute the reflection coefficents of the numerically determined time
% domain signals and determine errors %%%

[contspec_default, bs] = mex_fnft_nsev(q_via_default_method, T, ...
    XI_plot, kappa, 'M', D);
assert(length(bs) == 0);

[contspec_alternative, bs] = mex_fnft_nsev(q_via_alternative_method, T, ...
    XI_plot, kappa, 'M', D);
assert(length(bs) == 0);

eps_xi = (XI_plot(2) - XI_plot(1))/(D - 1);
xi = XI_plot(1) + (0:D-1)*eps_xi;
contspec_exact = complex(200*exp(-xi.^4));

error_default_method = ...
    norm(contspec_exact - contspec_default)/norm(contspec_exact)
error_alternative_method = ...
    norm(contspec_exact - contspec_alternative)/norm(contspec_exact)

%%% Plot results %%%

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

subplot(2,2,1);
semilogy(xi, abs(contspec_default), xi, abs(contspec_alternative), xi, abs(contspec_exact), '--k');
xlim(XI_plot);
ylim([1e-4 2]*max(abs(contspec_exact)));
xlabel('\xi');
ylabel('|reflection coefficent(\xi)|');
legend('Via default method', 'Via alternative method', 'Specification', ...
    'location', 'south');
grid on;

subplot(2,2,3);
plot(xi, angle(contspec_default), xi, angle(contspec_alternative), xi, angle(contspec_exact), '--k');
xlim(XI_plot);
ylim(1.1*pi*[-1 1]);
xlabel('\xi');
ylabel('\anglereflection coefficent(\xi)');
grid on;

eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;

subplot(2,2,2);
plot(t, abs(q_via_default_method), t, abs(q_via_alternative_method));
xlabel('t');
ylabel('|q(t)|');
grid on;

subplot(2,2,4);
plot(t, angle(q_via_default_method), t, angle(q_via_alternative_method));
ylim(1.1*pi*[-1 1]);
xlabel('t');
ylabel('\angleq(t)');
grid on;
