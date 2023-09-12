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
% Sander Wahls (KIT) 2023.

%% Periodic nonlinear Fourier transform of the Zabusky-Kruskal signal

% This example recreates the right panel of Figure 1 in the paper "Hidden 
% solitons in the Zabusky–Kruskal experiment: Analysis using the periodic, 
% inverse scattering transform" by I.C. Christov, Math. Comput. Simul.
% 82, 2012, 1069-1078, https://dx.doi.org/10.1016/j.matcom.2010.05.021
%
% It furthermore shows a plot with the Floquet determinant, band edges,
% etc.

clear all
close all

%% Set parameters (page 1071 in the paper)

del = 0.022;    % Nonlinearity parameter for the partly normalized KdV 
L = 2;          % Period [cm]

%% Compute the parameters of the corresponding dimensional KdV equation (page 1071)

g = 981;                        % gravity [cm/s^2]
be = del^2;
h = (6*del^2/sqrt(g))^(2/5);    % water depth [cm]
c0 = (6*del^2*g^2)^(1/5);
al = 3/2*(6*del^2/g^3)^(-1/5);
lam = al/(6*be);

%% Generate the signal

D = 256;
x = linspace(0, L, D+1);    % note the D+1
x = x(1:D);                 % exclude x=L because it already belongs to the
                            % next period
X = [x(1) x(end)];
a = 1/al;                   % below Eq. 7
u = a*sin(2*pi/L*x);        % Eq. 7
q = lam*u;

%% Compute various representations of the nonlinear Fourier spectrum

E = [-2000 2500];   % spectral interval
grid_spacing = 1;   % Max. allowed distance between consecutive grid
                    % points on the spectral interval

ampmodfreqs = mex_fnft_kdvp(q, [0 L], E, 'grid_spacing', grid_spacing, 'mstype_amplitudes_moduli_freqs');
A = ampmodfreqs(1:3:end)/(a*lam);    % amplitudes
m = ampmodfreqs(2:3:end);            % moduli
F = ampmodfreqs(3:3:end)/pi;         % nonlinear frequencies (not used here)

%% Recreate Fig. 1 (right panel) from the paper

figure
subplot(1,2,1)
i_ref = find(m>=0.99, 1, 'last');
plot(2*pi*(1:length(A))/L, A, '-o')
hold on
plot(2*pi*(1:length(A))/L, m, '-s')
plot(2*pi/L*[i_ref i_ref], ylim, '--k');
xlabel('2\pij/L [1/cm]')
xlim(2*pi/L*[1 length(A)])
grid minor
legend('Amplitudes [cm]', 'Moduli');

subplot(1,2,2)
plot(F, A, '-o')
hold on
plot(F, m, '-s')
xlim([min(F) max(F)])
xlabel('Wave numbers [1/cm]')
grid minor
legend('Amplitudes [cm]', 'Moduli');