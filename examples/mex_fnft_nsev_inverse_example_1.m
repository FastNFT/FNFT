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

% This examples demonstrates how the inverse nonlinear Fourier transform
% with respect to the nonlinear Schroedinger equation with vanishing
% boundary conditions can be computed using mex_fnft_nsev_inverse.
%
% The signal is a truncated soliton pulse given in Rourke and Morris,
% Phys. Rev. A 46(7), 1992, https://doi.org/10.1103/PhysRevA.46.3631
%
% Note that this is a difficult test case since the time domain signal is
% discontinuous and the reflection coefficient decays slowly.

clear all;
close all;

%%% Setup parameters %%%

T = [-2 2];     % Location of the 1st and last sample in the time domain
D = 2^10;       % Desired number of samples in the time domain
M = 2*D;        % Number of samples in the nonlinear frequency domain
[XI, xi] = mex_fnft_nsev_inverse_XI(D, T, M);
                % Location of the 1st and last sample in the nonlinear
                % frequency domain, as well as the grid of all locations
kappa = +1;     % Focusing nonlinear Schroedinger equation

%%% Setup the reflection coefficent %%%

al = 2;
be = 0.55;
gam = sqrt(abs(al)^2 + be^2);
contspec_fun = @(xi) complex(al./(xi - be*1j));
contspec_exact = contspec_fun(xi);

%%% Compute the corresponding time domain signal numerically %%%

bound_states = [1j*be];
normconsts = [-1j*al/(gam + be)];
q = mex_fnft_nsev_inverse(contspec_exact, XI, bound_states, normconsts, ...
    D, T, kappa);

%%% Compute the nonlinear Fourier transform of the numerically determined
% time domain signal %%%

[contspec, bs, nc] = mex_fnft_nsev(q, T, XI, kappa, 'M', M);

%%% Compute the exact values of the time domain signal using the formula
% in the reference given above %%%

eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;
q_exact = (t <= 0).*-2.0j*gam*al/abs(al).*sech(2*gam*t + atanh(be/gam));

%%% Plot results %%%

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

subplot(2,2,1);
plot(t, abs(q), t, abs(q_exact), '--k');
xlabel('t');
ylabel('|q(t)|');
legend('Numerically', 'Exact solution');
grid on;

subplot(2,2,3);
plot(t, angle(q), t, angle(q_exact), '--k');
ylim(1.1*pi*[-1 1]);
xlabel('t');
ylabel('\angleq(t)');
grid on;

subplot(2,2,2);
plot(xi, abs(contspec), xi, abs(contspec_exact), '--k');
xlabel('\xi');
ylabel('|reflection coefficent(\xi)|');
legend('Numerically', 'Specification');
grid on;

subplot(2,2,4);
plot(xi, angle(contspec), xi, angle(contspec_exact), '--k');
ylim(1.1*pi*[-1 1]);
xlabel('\xi');
ylabel('\anglereflection coefficent(\xi)');
grid on;
