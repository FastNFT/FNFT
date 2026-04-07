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

% This example should illustrate how to use the inverse kdvv transform

clear all;
close all;

% desired height of solitions:
desired_solitions = [4, 3, 2, 1];

% resulting bound states out of desired heights of solitions:
% - bound states have to be true and positive imaginary
bound_states = 1i*sqrt(desired_solitions ./2);

% desired norming constants
% defines how much the solitions are shifted towards each other
% - signs have to alternate regards to the order of the bound states
% - sign of the normconst for the biggest eigenvalue has to be positive
norming_constants = complex([1, -1, 1, -1].*[0.00001, 0.1, 1, 10]);

% Number of samples of the output of the inverse kdvv
D = 1001;

% Area, for which the output has to be computed
% Can be asymmetric
T = [-10 10];

% defining the continuous spectrum:
% out of function, just for seek of completeness (state 04/2026)
contspec = [];
XI = [1e-6 1];

% calls function of c-library FNFT
q = mex_fnft_kdvv_inverse(contspec, XI, bound_states, norming_constants, D, T);
% the result consists of solitions with desired height (if the solitions are far
% away enough to each other)

% --- Plot the results ---
t = linspace(T(1), T(2), D);

% the output of the inverse kdvv is defined as a complex number (state 04/2026),
% however only the imaginary part should be zero and can be neglected
q = double(real(q(:)'));

figure;
plot(t, q);
title('output of inverse kdvv: time-domain');
xlabel('t');
ylabel('q(t)');
legend('q(t)');

figure;
stem(imag(bound_states),real(norming_constants), 'x');
title('Bound states and norming constants');
xlabel('bound states');
ylabel('norming constants');