% MEX_FNFT_KDVP_AMPMODFREQ Determines the amplitudes, moduli and frequencies
% of the hyperelliptic modes associated to the nonlinear Fourier transform
% for the Korteweg-de Vries equation with periodic boundaries.
%
%   ampmodfreq = MEX_FNFT_KDVP_AMPMODFREQ(main_spec)
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_kdvp_ampmodfreq.
%
% INPUTS
%   main_spec       Main spectrum computed with fnft_kdvp.
%
% OUTPUTS
%   ampmodfreq      Real row vector

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
% Sander Wahls (KIT) 2025.
