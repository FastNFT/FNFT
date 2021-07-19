% MEX_FNFT_MANAKOVV Fast nonlinear Fourier transform for the nonlinear
% Schroedinger equation with vanishing boundaries.
%
%   contspec = MEX_FNFT_MANAKOVV(q1, q2, T, XI, kappa)
%   contspec = MEX_FNFT_MANAKOVV(q1, q2, T, XI, kappa,
%       OPTIONAL INPUTS)
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_manakovv.
%
% INPUTS
%   q1              Complex row vector of length D>1
%   q2              Complex row vector of length D>1
%   T               Real 1x2 vector
%   XI              Real 1x2 vector
%   kappa           +1.0 or -1.0
%
% OPTIONAL INPUTS
%   It is possible to provide additional inputs. These come either in the
%   form of a single string or of a string followed by a value.
%
%   'M'             The length of the vector contspec. Followed by the
%                   desired value, a positive integer. If not provided, M
%                   is set to D
%   'RE'            Use Richardson extrapolation to improve accuracy. The
%                   approximations of the nonlinear Fourier spectrum are
%                   calcuated using all given samples and again with half of
%                   the samples. The two approximations are combined
%                   through the idea of Richardson extrpolation to
%                   hopefully obtain a more accurate approximation. Note
%                   that in certain situations such as discontinuous signals,
%                   enabling this option may result in worse accuracy than
%                   without it.
%   'cstype_ab'     Returns values of a(xi), b1(xi) and b2(xi) individually 
%                   instead of the values of reflection coefficient 
%                   b1(xi)/a(xi) and b2(xi)/a(xi)
%   'cstype_both'   Returns both the reflection coefficients and a(xi),
%                   b1(xi) and b2(xi)
%   'quiet'         Turns off messages generated by then FNFT C library.
%                   (To turn off the messages generated by the mex
%                   interface functions, use MATLAB's warning and error
%                   commands instead.)
%
%       The following options specify the use of fast discretization schemes,
%       characterised by the Order of the base method (O), Polynomial degree per
%       potential sample (P), and order of occuracy of the splitting scheme (S).
%   'discr_2split3A' O = 2; P = 3;   S=3.
%   'discr_2split3B' O = 2; P = 3;   S=3.
%   'discr_2split4A' O = 2; P = 4;   S=4.
%   'discr_2split4B' O = 2; P = 2;   S=4.
%   'discr_2split6B' O = 2; P = 6;   S=6.
%   'discr_4split4A' O = 4; P = 4;   S=4.
%   'discr_4split4B' O = 4; P = 2;   S=4.
%   'discr_4split6B' O = 4; P = 7;   S=6.
%   'discr_FTES4_suzuki' O = 4; P = 2;   S=4.
%   TODO: check O, P, S
% 
%       The following options specify slow discretization schemes.
%   'discr_BO'      Use the second-order method by Boffetta-Osborne.
%                   Requires one matrix exponential per sample.
%   'discr_CF4_2'   Use fourth-order commutator-free exponential integrator
%                   which requires two matrix exponentials per sample.
%
% OUTPUTS
%   contspec        Complex row vector of length 2*M (default), length 3*M
%                   if 'cstype_ab' is set and 5*M if cstype_both is set

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
% Contributor:
% Lianne de Vries (TU Delft) 2021.