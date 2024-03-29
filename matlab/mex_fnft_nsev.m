% MEX_FNFT_NSEV Fast nonlinear Fourier transform for the nonlinear
% Schroedinger equation with vanishing boundaries.
%
%   [contspec, bound_states, normconsts] = MEX_FNFT_NSEV(q, T, XI, kappa)
%   [contspec, bound_states, normconsts] = MEX_FNFT_NSEV(q, T, XI, kappa,
%       OPTIONAL INPUTS)
%
% DESCRIPTION
%   Provides an interface to the C routine fnft_nsev.
%
% INPUTS
%   q               Complex row vector of length D>1
%   T               Real 1x2 vector
%   XI              Real 1x2 vector
%   kappa           +1.0 or -1.0
%
% OPTIONAL INPUTS
%   It is possible to provide additional inputs. These come either in the
%   form of a single string or of a string followed by a value.
%
%   'M'             The length of the vector contspec. Followed by the
%                   desired value, a positive integer.
%   'bsloc_fasteigen'   Use fast eigenvalue method to locate bound states.
%                   This method is very reliable, but requires O(D^2)
%                   FLOPs. This input is not followed by a value.
%   'bsloc_newton'  Use Newton's method to locate bound states. This method
%                   is reliable if good intial guesses for the bound states
%                   are known. Followed by a complex row vector of K
%                   initial guesses. It requires O(niter KD) FLOPs.
%   'bsloc_subsamp_refine' Use a mixed method to locate bound states. First
%                   get initial guesses for the bound states by applying
%                   the 'fasteigen' method to a subsampled version of the
%                   signal. Then refine using 'Newton' based on the full
%                   signal. This method is reliable if D is not too low.
%                   It requires O(D log^2 D + niter K D) FLOPs if Dsub (see
%                   below) is set by the algorithm. Not followed by a
%                   value.
%   'bsloc_niter'   Number of iterations to be carried by Newton's method.
%                   Followed by a non-negative integer.
%   'bsloc_tol'     Tolerance used for stopping Newton's method. Followed
%                   by a non-negative real number.
%   'bsloc_Dsub'    The desired number of samples for the subsampled signal
%                   in the 'subsamp_refine' method. Less samples in the
%                   subsampled stage result in faster execution time, but
%                   might also lead to a loss of precision. Followed
%                   by scalar value. Note that the routine treates Dsub as
%                   as an indication. The actually used value might differ.
%   'bsfilt_none'   Do not filter bound states at all.
%   'bsfilt_basic'  Basic bound state filtering. Removes duplicates and
%                   bound states in the lower half plane.
%   'bsfilt_full'   Full bound state filtering. Additionally removes
%		            physically implausible bound states.
%   'bsfilt_manual' Manual filtering. Remove bound states outside of a box
%                   provided by the user. Must be followed by a real vector
%                   of the form [min_real, max_real, min_imag, max_imag]
%                   that specifies the box.
%   'RE'            Use Richardson extrpolation to improve accuracy. The
%                   approximations of the nonlinear Fourier spectrum are
%                   calcuated using all given samples and again with half of
%                   the samples. The two approximations are combined
%                   through the idea of Richardson extrpolation to
%                   hopefully obtain a more accurate approximation. Note
%                   that in certain situations such as discontinuous signals,
%                   enabling this option may result in worse accuracy than
%                   without it.
%   'dstype_residues'   Return residues instead of norming constants
%   'cstype_ab'     Returns values of a(xi) and b(xi) individually instead
%                   of the values of reflection coefficient b(xi)/a(xi)
%   'skip_cs'       Skip computation of the continuous spectrum.
%   'skip_bs'       Skip computation of the bound states. Implies 'skip_nc'.
%   'skip_normalization' Skips normalization steps during the computation of
%                   scattering data. Faster but more prone to numerical overflow.
%   'skip_nc'       Skip computation of the norming constants.
%   'quiet'         Turns off messages generated by then FNFT C library.
%                   (To turn off the messages generated by the mex
%                   interface functions, use MATLAB's warning and error
%                   commands instead.)
%
%       The following options specify the use of fast discretization schemes,
%       characterised by the Order of the base method (O), Polynomial degree per
%       potential sample (P), and order of occuracy of the splitting scheme (S).
%   'discr_modal'    O = 2; P = 1;   S=2.
%   'discr_2split1A' O = 2; P = 1;   S=1.
%   'discr_2split1B' O = 2; P = 1;   S=1.
%   'discr_2split2A' O = 2; P = 1;   S=2.
%   'discr_2split2B' O = 2; P = 1;   S=2.
%   'discr_2split2S' O = 2; P = 1;   S=2.
%   'discr_2split3A' O = 2; P = 3;   S=3.
%   'discr_2split3B' O = 2; P = 3;   S=3.
%   'discr_2split3S' O = 2; P = 2;   S=3.
%   'discr_2split4A' O = 2; P = 4;   S=4.
%   'discr_2split4B' O = 2; P = 2;   S=4.
%   'discr_2split5A' O = 2; P = 15;  S=5.
%   'discr_2split5B' O = 2; P = 15;  S=5.
%   'discr_2split6A' O = 2; P = 12;  S=6.
%   'discr_2split6B' O = 2; P = 6;   S=6.
%   'discr_2split7A' O = 2; P = 105; S=7.
%   'discr_2split7B' O = 2; P = 105; S=7.
%   'discr_2split8A' O = 2; P = 24;  S=8.
%   'discr_2split8B' O = 2; P = 12;  S=8.
%   'discr_4split4A' O = 4; P = 4;   S=4.
%   'discr_4split4B' O = 4; P = 2;   S=4.
%
%       The following options specify slow discretization schemes.
%   'discr_BO'      Use the second-order method by Boffetta-Osborne.
%                   Requires one matrix exponential per sample.
%   'discr_CF4_2'   Use fourth-order commutator-free exponential integrator
%                   which requires two matrix exponentials per sample.
%   'discr_CF4_3'   Use fourth-order commutator-free exponential integrator
%                   which requires three matrix exponentials per sample.
%   'discr_CF5_3'   Use fifth-order commutator-free exponential integrator
%                   which requires three matrix exponentials per sample.
%   'discr_CF6_4'   Use sixth-order commutator-free exponential integrator
%                   which requires four matrix exponentials per sample.
%   'discr_ES4'     Use fourth-order exponential integrator
%                   which requires one matrix exponential per sample.
%   'discr_TES4'    Use fourth-order exponential integrator
%                   which requires three matrix exponentials per sample.
%
% OUTPUTS
%   contspec        Complex row vector of length M, the default is M=D
%                   (or [] if skipped)
%   bound_states    Complex row vector of length K (or [] if skipped)
%   normconsts      Complex row vector of length K (or [] if skipped)


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
% Sander Wahls (TU Delft) 2017-2018, 2022.
% Shrinivas Chimmalgi (TU Delft) 2019-2020.
% Peter J. Prins (TU Delft) 2021.
% Sander Wahls (KIT) 2023.
