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
% Fabian Fischer (Hiwi KIT) 2026.

% This example shows a fast way to try different values (bound states, norming constants)
% for the inverse kdvv transform, because this script sorts and considers the alternating 
% signs of the norming constants automatically.

clear all;
close all;

% === modify here ====================================================================================
% a solition in the output of inverse kdvv is defined by the elements of both arrays at the same index
soliton_amplitudes = [9, 2, 4, 3, 7];
solition_shifts = [10000, 0.1, 1, 0.00001, 2];
% ====================================================================================================


% Check if both arrays contains the same number of elements
if length(solition_shifts) ~= length(soliton_amplitudes)
    error('Number of desired solitions is not equal to number of values given for shifting!');
end

% --- determining sorted bound states and norming constants out of desired values ---

number_values = length(soliton_amplitudes);

tmp_bound_states = 1i*sqrt(soliton_amplitudes ./2);

% sorting
[bound_states, indices] = sort(tmp_bound_states, "descend");
tmp_norming_constants = solition_shifts(indices);

% alternating sequence 1 and -1 as elements and with 1 as first element
tmp_seq = (-1).^( 0:(number_values-1) );

norming_constants = complex(tmp_seq.*tmp_norming_constants);

% Continuous spectrum functionality has not yet been implemented.
% Hand over empty arrays
contspec = [];
XI = [];

% inverse kdvv
D = 1001;
T = [-12 8];
q = mex_fnft_kdvv_inverse(contspec, XI, bound_states, norming_constants, D, T);


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