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

close all;

contspec = [];
XI = [0 1];
% bound_states = 1i*[1, 2, 3, 4] /2;
bound_states = 1i*[1] /2;
% normconsts = complex([-1, 1, -1, 1].*[10, 0.1, 1, 0.00001]);
normconsts = complex([1].*[10]);
D = 1001;
T = [-20 20];
q = mex_fnft_kdvv_inverse(contspec, XI, bound_states, normconsts, D, T);

t = linspace(T(1), T(2), D);
plot(t, q)
xlabel('t')
ylabel('q(t)')
