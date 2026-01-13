% This file is part of FNFT.
%
% FNFT is free software; you can redistribute it and/or
% modify it under the terms of the version 2 of the GNU General
% Public License  as published by the Free Software Foundation.
%aux
% FNFT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Contributors:
% Sander Wahls (KIT) 2025.

%% Periodic nonlinear Fourier transform of two sine waves

% This example demonstrates that the nonlinear Fourier transform (NFT) for the
% periodic KdV equation reduces to the FFT in the linear limit (i.e., not too
% large amplitudes). The signal is a sum of two sine waves. The NFT and FFT
% both detect the same frequencies and amplitudes in this case.

clear all
close all

%% Set parameters

P = 2;                  % Period of the time domain signal
D = 256;                % Number of time domain samples
A1 = 2;                 % Amplitude of the first sine wave
A2 = -1;                % Amplitude of the second sine wave
f1 = 8/P;               % Frequency of the first sine wave
f2 = 16/P;              % Frequency of the second sine wave

E = [-100 800];         % Spectral interval (increase the upper bound
                        % to cover higher frequencies)
grid_spacing = 0.01;    % Max. allowed distance between consecutive grid
                        % points on the spectral interval

%% Generate the time-domain signal

t = linspace(0, P, D+1);
t = t(1:end-1);
q = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);

%% Compute the single-sided FFT for real signals
Q = fft(q);
Q = 2 * Q(1:D/2)/D;

dt = t(2)-t(1);
assert(mod(D,2)==0)
f = (0:D/2-1)/(D*dt);

%% Compute the periodic KdV-NFT

[main_spec, ~, ~] = mex_fnft_kdvp(q, [0 P], E, 'grid_spacing', grid_spacing, 'keep_degenerate');

ampmodfreqs = mex_fnft_kdvp_ampmodfreq(main_spec);
A = ampmodfreqs(1:3:end);    % amplitudes
m = ampmodfreqs(2:3:end);    % moduli
F = ampmodfreqs(3:3:end);

%% Plot results

figure
subplot(2,1,1)
plot(t, q);
xlabel('t')
ylabel('q(t)')
title('Time-domain signal')

subplot(2,1,2)
stem(F, A, 'o')
hold on
plot(f(1:D/2), abs(Q), 'x')
hold off
grid on
xlabel('f')
ylabel('A')
legend('NFT', 'FFT')
xlim([0 20])
title('Nonlinear and linear amplitude spectrum')
