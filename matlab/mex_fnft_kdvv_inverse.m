% MEX_FNFT_KDVV_INVERSE Fast inverse nonlinear Fourier transform for the
% Korteweg-de  Vries equation with vanishing boundaries.
%
%   q = MEX_FNFT_KDVVV_INVERSE(contspec, XI, bound_states, ...
%                             normconsts_or_residuals, D, T);
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_kdvv_inverse.
%
% INPUTS
%   contspec        Complex row vector of length M>=D, contains the samples
%                   of the reflection coefficient, the b-scattering 
%                   coefficient or the inverse Fourier transform of the 
%                   b-scattering coefficient on an equidistant grid.
%                   Pass [] if the continuous spectrum is zero 
%                   (i.e., a multi-soliton is desired)
%   XI              Real 1x2 vector, contains the location of the first and
%                   the last sample in contspec
%   bound_states    Complex row vector, contains the desired bound states.
%                   Pass [] if the discrete spectrum is empty.
%   normconsts_or_residues Complex row vector, same length as bound_states.
%                   Contains the corresponding norming constants (default) or,
%                   if the corresponding option is passed, residues. Pass []
%                   if the discrete spectrum is empty.
%   D               Real scalar, number of time domain samples; must be a
%                   positive power of two
%   T               Real 1x2 vector, contains the location of the first and
%                   the last sample in q
%   kappa           +1.0 or -1.0
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