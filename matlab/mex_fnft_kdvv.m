% MEX_FNFT_KDVV Fast nonlinear Fourier transform for the Korteweg-de Vries
% equation with vanishing boundaries.
%
%   contspec = MEX_FNFT_KDVV(q, T, XI, kappa)
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_kdvv.
%
% INPUTS
%   q               Complex row vector of length D=2^n
%   T               Real 1x2 vector
%   XI              Real 1x2 vector
%   kappa           +1.0 or -1.0
%
% OUTPUTS
%   contspec        Complex row vector of length D

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
% Sander Wahls (TU Delft) 2018.
