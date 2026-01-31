% This file is part of FNFT.
%
% FNFT is free software; you can redistribute it and/or
% modify it under the terms of the version 2 of the GNU General
% Public License  as published by the Free Software Foundation.
%aux
% FNFT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Contributors:
% Sander Wahls (KIT) 2026

% This script builds the mex interface using GNU Octave.

clear

if ~exist('OCTAVE_VERSION', 'builtin')
  error('This script is indented for GNU Octave only.');
end

mex_src_files = dir('mex_fnft_*.c');
if length(mex_src_files) == 0
  error('No FNFT mex source files found. Please cd into FNFT''s matlab folder before running this script.')
end

if isunix % Linux
  lib_file = '../lib/libfnft.so';
elseif ispc % Windows
  lib_file = '..\build\libfnft.dll';
elseif ismac
  error('Mac support is not available');
end
if ~isfile(lib_file)
  error(sprintf('Could not locate the FNFT library file "%s". Please build FNFT first.', lib_file));
end

fprintf('Building mex files:\n');
for f = mex_src_files'
  fprintf('  %s\n', f.name);
  mkoctfile('--mex', f.name, '-DSKIP_MATRIX_H', '-I../include',
  '-I../include/private', '-I../include/3rd_party/kiss_fft', '-L../lib', '-lfnft');
end

if isunix
  f = dir(lib_file);
  lib_full = fullfile(f.folder, f.name);
  fprintf('\nIn order to make sure that the mex files can find the FNFT share library file,\n');
  fprintf('Octave must be told where to find this file. One simple way of doing this under\n');
  fprintf('Linux is to start Octave from the shell using the command\n\n');
  fprintf('  LD_PRELOAD="%s" octave --gui\n\n', lib_full);
end
