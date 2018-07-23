% MEX_FNFT_VERSION Provides the version of the used FNFT C library.
%
%   MEX_FNFT_VERSION();
%   [major, minor, patch, suffix] = MEX_FNFT_VERSION();
%
% DESCRIPTION
%   Provides the version of the underlying FNFT library. When called
%   without return value, the version is simply printed to the screen.
%   When called with four return values, it returns the constants
%   that are used to build the version string.
%
% OUTPUTS
%   major              Value of the constant FNFT_VERSION_MAJOR
%   minor              Value of the constant FNFT_VERSION_MINOR
%   patch              Value of the constant FNFT_VERSION_PATCH
%   suffix             Value of the constant FNFT_VERSION_SUFFIX

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
