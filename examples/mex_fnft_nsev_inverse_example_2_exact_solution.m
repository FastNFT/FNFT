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

% The formulas used here are given in the paper
% Frumin et al., J. Opt. Soc. Am. B 32(2), 2015.

function [q_exact, contspec_exact] = ...
    mex_fnft_nsev_inverse_example_2_exact_solution(t, xi)
try
    Q = 5;
    GAM = 1/25;
    F = 1.5;
    
    q_exact = -conj(Q/GAM*sech(t/GAM).^(1-2j*F));
    
    cgamma = @(z) double(gamma(sym(z)));
    d = 0.5 + 1i*(xi*GAM-F);
    fp = 0.5 - 1i*(xi*GAM+sqrt(F^2+Q^2));
    fm = 0.5 - 1i*(xi*GAM-sqrt(F^2+Q^2));
    gp = 1 - 1i*(F+sqrt(F^2+Q^2));
    gm = 1 - 1i*(F-sqrt(F^2+Q^2));
    contspec_exact = -2^(-2i*F)*Q*cgamma(d).*cgamma(fm).* ...
        cgamma(fp)./(cgamma(conj(d)).*cgamma(gm).*cgamma(gp));
catch
    error('This function requires the Symbolic Math Toolbox in order to compute the complex gamma function.');
end
end