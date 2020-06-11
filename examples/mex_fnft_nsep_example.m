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
% Shrinivas Chimmalgi (TU Delft) 2020.

% This examples demonstrates how the nonlinear Fourier transform with
% respect to the nonlinear Schroedinger equation with periodic boundary
% conditions can be computed using mex_fnft_nsep. The signal is the same
% plane wave as in the paper https://doi.org/10.1109/TIT.2015.2485944

clear all;
close all;

%%% Setup parameters %%%

T = [0, 2*pi];  % T(1) is the beginning of the period, T(2) is the end
D = 2^8+1;        % number of samples
kappa = +1;     % focusing nonlinear Schroedinger equation  

%%% Setup the signal %%%

t = linspace(T(1),T(2),D);
% Note: in the previous version of FNFT The location of the 1st sample is T(1), but the location of the
% sample if T(2) - (T(2)-T(1)/D, i.e. t = T(1) + (0:D-1)*(T(2) - T(1))/D.
q = 3*exp(3j*t);

%%% Compute the nonlinear Fourier transform %%%

[main_spec, aux_spec] = mex_fnft_nsep(q, T, kappa);

%%% Plot results %%%

figure;
plot(t, real(q), t, imag(q));
title('Time-domain');
xlabel('t');
ylabel('q(t)');
legend('Real part', 'Imaginary part');

figure;
plot(real(main_spec), imag(main_spec), 'x');
title('Main Spectrum');
xlabel('Re(\lambda_k)');
ylabel('Im(\lambda_k)');
xlim([-4 4]);
ylim([-4 4]);
axis equal;
grid on;

figure;
plot(real(aux_spec), imag(aux_spec), '+r');
title('Auxiliary Spectrum');
xlabel('Re(\lambda_k)');
ylabel('Im(\lambda_k)');
xlim([-4 4]);
ylim([-4 4]);
axis equal;
grid on;
