% MEX_FNFT_NSEV_INVERSE Provides the range in the nonlinear frequency
% domain that has to be used with MEX_FNFT_NSEV_INVERSE.
%
%   q = MEX_FNFT_NSEV_XI(contspec, T, D, XI, kappa);
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_nsev_inverse.
%
% INPUTS
%   T               Real 1x2 vector
%   D               Real scalar, must be a positive power of two
%   XI              Real 1x2 vector
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
