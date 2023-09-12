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
% Sander Wahls (KIT) 2023.

%% Periodic nonlinear Fourier transform of a Korteweg-de Vries soliton

% This example recreates Figure 3 from the paper "The solitons of Zabusky
% and Kruskal revisited: Perspective in terms of the periodic spectral
% transform" by Osborne and Bergamasco, Physica D: Nonlin. Phen. 18(1-3),
% Jan. 1986, p. 26-46, https://doi.org/10.1016/0167-2789(86)90160-0
%
% It furthermore shows a plot with the soliton amplitudes.

clear all
close all

%% Set parameters

g = 981;    % Gravity [cm/s^2]
h = 5;      % Water depth [cm]
u1 = 4;     % Soliton amplitude [cm]
L = 100;    % Period [cm]

%% Compute KdV and normalization parameters (Page 27 of the paper)

c0 = sqrt(g*h);
al = 3*c0/(2*h);
be = c0*h^2/6;
lam = al/(6*be);

%% Generate the signal (Eq. 52 in the paper)

D = 256;                    % number of samples
x = linspace(0, L, D+1);    % note the D+1
x = x(1:D);                 % exclude x=L because it already belongs to the
                            % next period
L1 = sqrt(12*be/(al*u1));
q = lam*u1*sech((x-L/2)/L1).^2;

%% Compute various representations of the nonlinear Fourier spectrum

E = [-0.03 0.015];         % Spectral interval
R = 10000;                 % Number of grid points for the Floquet diagram
grid_spacing = 0.0001;     % Max. allowed distance between consecutive grid
                           % points on the spectral interval

floq_det =              mex_fnft_kdvp(q, [0 L], E, 'mstype_floquet', R);
[main_spec, aux_spec] = mex_fnft_kdvp(q, [0 L], E, 'grid_spacing', grid_spacing);
bands =                 mex_fnft_kdvp(q, [0 L], E, 'grid_spacing', grid_spacing, 'mstype_openbands');
ampmodfreq =            mex_fnft_kdvp(q, [0 L], E, 'grid_spacing', grid_spacing, 'mstype_amplitudes_moduli_freqs');

%% Compute the soliton amplitudes

A = ampmodfreq(1:3:end)/lam;    % amplitudes
m = ampmodfreq(2:3:end);        % moduli
K = ampmodfreq(3:3:end)/pi;     % nonlinear wave numbers
A_sol = A(K<0);                 % soliton amplitudes
K_sol = K(K<0);                 % soliton wave numbers

%% Recreate Fig. 3 from the paper

figure
subplot(2,1,1)
plot(x, q/lam)
xlabel('x [cm]')
ylabel(['\eta(x) [cm]'])
grid on

ix = abs(floq_det) > 1;
floq_scaled = real(floq_det);
floq_scaled(ix) = sign(real(floq_det(ix))).*(1+log(abs(floq_det(ix))));
Es = linspace(E(1), E(2), R);
B = [bands(1:2:end) ; bands(2:2:end)];

subplot(2,1,2)
plot(Es, floq_scaled)
xlim(E);
hold on
plot(B, 0*B, '-g', 'LineWidth', 1.2)
plot(main_spec(1:2:end), main_spec(2:2:end), 'ro')
plot(aux_spec, 0*aux_spec, 'rs')
plot([E(1) E(2)], [1 1], '-k')
plot([E(1) E(2)], [-1 -1], '-k')
hold off
grid on
xlabel('E');
legend('Floquet det. (scaled)', 'Open bands', 'Main spectrum', 'Aux. spectrum');

%% Plot the soliton amplitudes

figure;
stem(K_sol, A_sol, '-o');
xlabel('Wave number [1/cm]')
ylabel('Amplitude [cm]')
title('Soliton spectrum')
grid on
