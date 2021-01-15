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
% Peter J. Prins (TU Delft) 2021.

% This examples demonstrates how the nonlinear Fourier transform with
% respect to the Korteweg-de Vries equation with vanishing boundary
% conditions can be computed using mex_fnft_kdvv. The signal is a squared
% sech pulse.

clear all;
close all;

%%% Setup parameters %%%


D = 2^8;            % location of the 1st and last sample in the time domain
T = [-10, 10];      % number of samples
XI = [1e-10, 6];   % location of the 1st and last sample in the xi-domain

%%% Setup the signal %%%

ep_t = (T(2) - T(1)) / (D - 1);
t = T(1):ep_t:T(2);
q = 16*sech(t).^2;

%%% Compute the nonlinear Fourier transform %%%

[contspec, bound_states, norming_constants] = mex_fnft_kdvv(q, T, XI);
% mex_fnft_kdvv has many options => run "help mex_fnft_kdvv" to learn more

%%% Plot the results %%%

ep_xi = (XI(2) - XI(1)) / (D - 1);
xi = XI(1):ep_xi:XI(2);

figure;
plot(t, q);
title('Time-domain');
xlabel('t');
ylabel('q(t)');
legend('q(t)');

figure;
hmag=subplot(2,1,1);
semilogy(xi, abs(contspec));
title('Continuous spectrum');
xlabel('\xi');
ylabel('|r(\xi)|');
hang=subplot(2,1,2);
plot(xi, angle(contspec));
xlabel('\xi');
ylabel('\angle r(\xi)');
linkaxes([hmag,hang],'x');

figure;
stem(imag(bound_states),norming_constants, 'x');
title('Bound states and norming constants');
xlabel('\lambda_k/j');
ylabel('b(\lambda_k)');
axis([0,4,-2,2]);