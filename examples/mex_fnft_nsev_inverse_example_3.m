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
% mex_finft_nsev_inverse also supports the inversion of a continuous
% spectrum that is specificied using the inverse Fourier transform B(tau)
% of b(xi), which is useful for b-modulation (see
% https://doi.org/10.1109/ECOC.2017.8346231).
%
% The test signal is a shifted sech pulse. The formula for b(xi) for an
% unshifted sech can be found in Satsuma and Yajima, Proc. Theor. Phys.
% Suppl. 55(1), 1974. The effect of a a time domain shift on b(xi) is
% given in Yousefi and Kschischang, IEEE Trans. Inform. Theor. 60(7), 2014.

clear all;
close all;

%%% Setup parameters %%%

T = [-20 20];   % Location of the 1st and last sample in the time domain.
                % Note that we need -T(1)=T(2) for the inversion of B(tau)
D = 2^8;        % Desired number of samples in the time domain
XI = [-1 1];    % Not used if the continuous spectrum is described using
                % B(tau). Any vector with X(1)<XI(2) will do.
kappa = +1;     % Focusing nonlinear Schroedinger equation

%%% Setup the signal and the continuous spectrum %%%

A = 0.45; % should be >0 and <0.5 to avoid solitons
q_exact_fun = @(t) 1j*A*sech(t-1);
b_fun = @(xivec) exp(-2j*xivec).*1j*sin(pi*A)./cosh(pi*xivec);
B_fun = @(tauvec) 1j/(2*pi)*sin(pi*A)*sech((tauvec-2)/2);

%%% Setup the grids %%%

ep_t = (T(2) - T(1))/(D-1);
t = T(1) + (0:D-1)*ep_t;
tau = 2*t;

%%% Perform the inverse nonlinear Fourier transform and compare the result
% with the known closed-form formula for q(t).

B_vals = complex(B_fun(tau));
q = mex_fnft_nsev_inverse(B_vals, T, D, XI, kappa, 'cstype_B_of_tau');

q_exact = q_exact_fun(t);
error_in_q = norm(q_exact(:) - q(:))/norm(q_exact)

%%% Plot the result %%%

figure;
subplot(2,1,1);
plot(t, abs(q_exact), t, abs(q), '--');
xlabel('t');
ylabel('|q(t)|');
legend('Formula', 'Numerically');
grid on;

subplot(2,1,2);
plot(t, angle(q_exact), t, angle(q), '--');
xlabel('t');
ylabel('\angleq(t)');
grid on;

%%% Verify numerically that B(tau) is the inverse Fourier transform of
% b(xi) %%%

B_numerically_fun1 = @(tau) ...
    integral(@(xi) b_fun(xi).*exp(1j*xi*tau), -inf, inf)/(2*pi);
B_numerically_fun = @(tauvec) arrayfun(B_numerically_fun1, tauvec);
mismatch_in_B = norm(B_vals - B_numerically_fun(tau))/norm(B_vals)