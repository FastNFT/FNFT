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
% Sander Wahls (KIT) 2026.
% Fabian Fischer (Hiwi KIT) 2026.

% This example should illustrate, how to use the inverse kdvv in combination 
% with the forward kdvv transform. For basic usage of the inverse kdvv, please
% please have a look at mex_fnft_kdvv_inverse_example_1.m (state 04/2026).

clear all;
close all;

contspec_initial = [];
XI = [1e-6 10];

bound_states_initial = 1i*sqrt([4, 3, 2, 1] ./2);
norming_constants_initial = complex([1, -1, 1, -1].*[10, 0.1, 1, 0.00001]);

D = 1001;
T = [-10 10];

% inverse kdvv transform
q = mex_fnft_kdvv_inverse(contspec_initial, XI, bound_states_initial, norming_constants_initial, D, T);

q = double(real(q(:)'));

% compute the nonlinear Fourier transform of the output of the inverse kdvv
[contspec_computed, bound_states_computed, norming_constants_computed] = mex_fnft_kdvv(q, T, XI);


% --- Plot the results ---
t = linspace(T(1), T(2), D);

figure;
plot(t, q);
title('output of inverse kdvv: time-domain');
xlabel('t');
ylabel('q(t)');
legend('q(t)');

% plotting the initial and the computed discrete spectrum
% the difference between the initial and the computed spectrum should be very small

figure;
stem(imag(bound_states_initial),real(norming_constants_initial), 'x', 'color', 'r');
hold on;
stem(imag(bound_states_computed),real(norming_constants_computed), 'x', 'color', 'g');
hold off;
title('initial and computed bound states and norming constants');
xlabel('bound states');
ylabel('norming constants');
legend('initial', 'computed');

% plotting of continuous spectrum
% values are small because initial continuous spectrum was assumed as zero
% and the computed continuous spectrum is only the result of numerical errors

ep_xi = (XI(2) - XI(1)) / (D - 1);
xi = XI(1):ep_xi:XI(2);

figure;
hmag=subplot(2,1,1);
semilogy(xi, abs(contspec_computed));
title('computed continuous spectrum');
xlabel('\xi');
ylabel('|r(\xi)|');
hang=subplot(2,1,2);
plot(xi, angle(contspec_computed));
xlabel('\xi');
ylabel('\angle r(\xi)');
linkaxes([hmag,hang],'x');

