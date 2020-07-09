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

% This example shows the feature of computing bound states using the slow 
% algorithms. The example is Example 1 from doi: 10.1109/ACCESS.2019.2945480.

clear;
close all;
clc;

%%% Setup parameters %%%

T = [-32,32];   % location of the 1st and last sample in the time domain
XI = [-10,10];  % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation
D = 2^10;

%%% Setup the signal %%%

qo =  5.4;
lam0 = 3;
q_fun = @(t) qo*sech(t).*exp(-2i*t*lam0); % signal function

%%% Exact values of the bound states %%%

bound_states_exact = lam0+1i*(qo+0.5-(1:floor(qo+0.5)));
normconsts_exact = (-1).^(1:floor(qo+0.5));

%%% Compute the discrete part of the nonlinear Fourier transform %%%

bound_states_guesses = bound_states_exact+0.035*exp(1i*pi*rand(1,5)); % Add some error
% NOTE: The initial guess has to be quite close to the actual value for the
% Newton method to converge.

t = linspace(T(1),T(2),D);
q = q_fun(t); % signal samples
[~,bound_states_computed,normconsts_computed]=mex_fnft_nsev(q, T, XI,... 
    kappa,'discr_CF4_2','bsloc_newton',bound_states_guesses,'skip_cs', 'bsloc_niter',20, 'bsfilt_manual', [2,4,0,6]);

%%% Plot results %%%

plot(bound_states_exact,'o')
hold on
plot(bound_states_guesses,'*')
plot(bound_states_computed,'xk')
xlim([2.8,3.2])
xlabel('Real part');
ylabel('Imaginary part');
legend('Exact bound state','Guess of bound state','Computed bound state');