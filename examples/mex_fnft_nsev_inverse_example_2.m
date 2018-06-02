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
% such a situation in the defocusing case.
%
% The example used here is given in the paper
% Frumin et al., J. Opt. Soc. Am. B 32(2), 2015.

clear all;
close all;

%%% Setup parameters %%%

T = [-1 1];     % Location of the 1st and last sample in the time domain
D = 2^12;       % Desired number of samples in the time domain
M = D;          % Number of samples in the nonlinear frequency domain
[XI, xi] = mex_fnft_nsev_inverse_XI(D, T, M);
                % Location of the 1st and last sample in the nonlinear
                % frequency domain, as well as the grid of all locations
kappa = -1;     % Defocusing nonlinear Schroedinger equation

%%% Define the desired continuous spectrum and compute the corresponding
% time domain signal numerically using two different methods %%%

try
    eps_xi = (XI(2) - XI(1))/(M - 1);
    Q = 5;
    GAM = 1/25;
    F = 1.5;
    xi = XI(1) + (0:(M-1))*eps_xi;
    cgamma = @(z) double(gamma(sym(z)));
    d = 0.5 + 1i*(xi*GAM-F);
    fp = 0.5 - 1i*(xi*GAM+sqrt(F^2+Q^2));
    fm = 0.5 - 1i*(xi*GAM-sqrt(F^2+Q^2));
    gp = 1 - 1i*(F+sqrt(F^2+Q^2));
    gm = 1 - 1i*(F-sqrt(F^2+Q^2));
    contspec_exact = -2^(-2i*F)*Q*cgamma(d).*cgamma(fm).* ...
        cgamma(fp)./(cgamma(conj(d)).*cgamma(gm).*cgamma(gp));
catch
    error('This example requires the Symbolic Math Toolbox in order to compute the complex gamma function.');
end

q_via_default_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa);
q_via_alternative_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa, 'csinv_a_from_b_iter');

%%% Compute the exact solution of the problem and determine errors %%%

eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;
q_exact = -conj(Q/GAM*sech(t/GAM).^(1-2j*F));

error_default_method = ...
    norm(q_exact - q_via_default_method)/norm(q_exact)
error_alternative_method = ...
    norm(q_exact - q_via_alternative_method)/norm(q_exact)

%%% Plot results %%%

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;

subplot(1,2,1);
semilogy(t, abs(q_via_default_method), t, abs(q_via_alternative_method),...
    t, abs(q_exact), '--k');
xlabel('t');
ylabel('|q(t)|');
grid on;
legend('Via default method', 'Via alternative method', 'Exact solution',...
    'location', 'south');

subplot(1,2,2);
plot(t, angle(q_via_default_method), t, angle(q_via_alternative_method));
ylim(1.1*pi*[-1 1]);
xlabel('t');
ylabel('\angleq(t)');
grid on;
