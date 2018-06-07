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
% mex_finft_nsev_inverse supports two different methods to invert the
% continuous spectrum in the defocusing case. The default method is usually
% fine, but sometimes the alternative method can give better results.
% This example illustrates such a situation.
%
% The example used here is given in the paper
% Frumin et al., J. Opt. Soc. Am. B 32(2), 2015.

clear all;
close all;

%%% Setup parameters %%%

T = [-1 1];     % Location of the 1st and last sample in the time domain
D = 2^10;       % Desired number of samples in the time domain
kappa = -1;     % Defocusing nonlinear Schroedinger equation

eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;

%%% Define the desired continuous spectrum and compute the corresponding
% time domain signal numerically using the default method %%%

M = 2*D;        % Number of samples in the nonlinear frequency domain
[XI, xi] = mex_fnft_nsev_inverse_XI(D, T, M);
                % Location of the 1st and last sample in the nonlinear
                % frequency domain, as well as the grid of all locations

fprintf('Computing exact solution symbolically - please wait ...');                
[q_exact, contspec_exact] = ...
    mex_fnft_nsev_inverse_example_2_exact_solution(t, xi);
fprintf('done\n');     

q_via_default_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa);

error_default_method = ...
    norm(q_exact - q_via_default_method)/norm(q_exact)

%%% Recompute the desired continuous spectrum (different grid) and compute
% the corresponding time domain signal numerically using the alternative
% method %%%

M = D;          % The alternative method does not support oversampling
[XI, xi] = mex_fnft_nsev_inverse_XI(D, T, M);

fprintf('Computing exact solution symbolically - please wait ...');                
[q_exact, contspec_exact] = ...
    mex_fnft_nsev_inverse_example_2_exact_solution(t, xi);
fprintf('done\n');     

q_via_alternative_method = mex_fnft_nsev_inverse(contspec_exact, T, D, ...
    XI, kappa, 'csmethod_tfmatrix_contains_ab_from_iter');

error_alternative_method = ...
    norm(q_exact - q_via_alternative_method)/norm(q_exact)

%%% Plot results %%%

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

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
