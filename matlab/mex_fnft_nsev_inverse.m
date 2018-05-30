% MEX_FNFT_NSEV_INVERSE Fast inverse nonlinear Fourier transform for the
% nonlinear Schroedinger equation with vanishing boundaries.
%
%   q = MEX_FNFT_NSEV(contspec, T, D, XI, kappa);
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_nsev_inverse.
%
% INPUTS
%   contspec        Complex row vector of length M>D, contains the samples
%                   of the reflection coefficient on an equidistant grid
%                   with end points given by XI below
%   T               Real 1x2 vector, contains the location of the first and
%                   the last sample in q
%   D               Real scalar, number of time domain samples; must be a
%                   positive power of two
%   XI              Real 1x2 vector, contains the location of the first and
%                   the last sample in contspec; must be computed with
%                   MEX_FNFT_NSEV_INVERSE_XI
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
% Sander Wahls (TU Delft) 2017-2018.
