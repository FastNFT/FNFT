% MEX_FNFT_KDVV_INVERSE Fast inverse nonlinear Fourier transform for the
% Korteweg-de  Vries equation with vanishing boundaries.
%
%   q = MEX_FNFT_KDVVV_INVERSE(contspec, XI, bound_states, ...
%                             norming_constants, D, T);
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_kdvv_inverse.
%
% INPUTS
%   contspec        Complex row vector of length M>=D, contains the samples
%                   of the reflection coefficient, the b-scattering 
%                   coefficient or the inverse Fourier transform of the 
%                   b-scattering coefficient on an equidistant grid.
%                   Pass [] if the continuous spectrum is zero.
%                   (i.e., a multi-soliton is desired)
%                   Note: continuous spectrum related functionality has not 
%                   yet been implemented (state 04/2026)! 
%                   If not empty array is handed over, error is returned!
%                   It will be implicitly assumed to \f$ 0 \f$ for all \f$ \xi \f$.
%   XI              Real 1x2 vector, contains the position of the first and the last
%                   sample of the continuous spectrum.
%                   Note: continuous spectrum related functionality has not 
%                   yet been implemented (state 04/2026)! 
%                   If not empty array is handed over, error is returned!
%   bound_states    Complex row vector. Bound states have to be positive, 
%                   purely imaginary numbers (lie on the upper half of the imaginary axis). 
%                   The bound states have to be in descending order. To add a solition 
%                   with height \f$ h_i \f$ the bound state have to be 
%                   \f$ \gamma_i = \sqrt{ h_i/2} \f$.
%   norming_constants Complex row vector, same length as bound_states.
%                   Contains the corresponding norming constants.
%                   The signs of the norming constants have to alternate regards to 
%                   the order of the bound states. The sign of the normconst for the 
%                   biggest eigenvalue has to be positive
%   D               Real scalar, number of time domain samples (i.e. number of samples 
%                   of the resulting signal q).
%   T               Real 1x2 vector, contains the location of the first and
%                   the last sample in q
%
% OUTPUTS
%   q               Complex row vector of length D

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