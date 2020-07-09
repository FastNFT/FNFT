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

% This example shows the feature of computing bound states using the phase 
% jump tracking algorithm. The example is Example 1 from doi: 10.1103/PhysRevE.96.063302

clear;
close all;
clc;

%%% Setup parameters %%%

T = [-350,350];   % location of the 1st and last sample in the time domain
XI = [-10,10];  % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation
D = 2^11;

%%% Exact values of the bound states %%%

K = 32; % Number of bound states
theta0=pi/3;
J=4;
L=14;
deltatheta=(pi-2*theta0)/(J-1);
theta=theta0+((1:J)-1)*deltatheta;
lactual=kron([1:L],exp(1i*theta));
bactual=exp(1i*pi*((1:J*L)-1)/(L*J));
bound_states_exact = lactual(1:K);
kap = 2*sqrt(sum(imag(bound_states_exact)));
bound_states_exact= bound_states_exact/kap;
normconsts_exact = bactual(1:K);

%%% Generating the multi-soliton %%%

q = mex_fnft_nsev_inverse([], [-1,1], bound_states_exact, normconsts_exact, D, T, kappa);

%%% Compute the discrete part of the nonlinear Fourier transform %%%

[~,bound_states_computed,normconsts_computed]=mex_fnft_nsev(q, T, XI,... 
    kappa,'discr_ES4','bsloc_pjt');

%%% Plot results %%%

plot(bound_states_exact,'o')
hold on
plot(bound_states_computed,'xk')
xlabel('Real part');
ylabel('Imaginary part');
legend('Exact bound state', 'Computed bound state');