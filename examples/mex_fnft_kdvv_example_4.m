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
% Shrinivas Chimmalgi (TU Delft) 2020.
% Sander Wahls (TU Delft) 2020.
% Peter J. Prins (TU Delft) 2021.

% This example shows the feature of computing bound states using the slow 
% algorithms.

clear;
close all;
clc;

%%% Setup parameters %%%

T = [0,16];     % location of the 1st and last sample in the time domain
XI = [0.5,23];  % location of the 1st and last sample in the xi-domain
D = 2^10;

%%% Setup the signal %%%

A = sym(15);
d = sym(0.5);
exp_t0 = sym(3000);
q_fun = @(t) A*sech((t-log(exp_t0))/d).^2; % signal function

%%% Exact values of the bound states %%%

bound_states_exact = [1i,3i];
normconsts_exact = [-9e6, 729e18];

%%% Compute the discrete part of the nonlinear Fourier transform %%%

t = linspace(T(1),T(2),D);
q = double(q_fun(t)); % signal samples
dxi = sqrt(max(q)) / 1000; % use approximately 1000 grid points for the
% bound state localization step
[~,bound_states_computed,normconsts_computed]=mex_fnft_kdvv(q, T, XI,... 
    'discr_CF4_2','bsloc_gridsearch_refine', 'grid_spacing', dxi, ...
    'skip_cs', 'bsloc_niter',20);

%%% Plot results %%%

plot(bound_states_exact,'o')
hold on
plot(bound_states_computed,'xk')
ylim([0,4])
xlabel('Real part');
ylabel('Imaginary part');
legend('Exact bound state','Computed bound state');