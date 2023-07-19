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
% Sander Wahls (TU Delft) 2017-2018, 2020.
% Shinivas Chimmalgi (TU Delft) 2019.
% Peter J. Prins (TU Delft) 2021.

% This example compares the accuracy-execution tradeoff of the various algorithms for
% the Korteweg-de Vries equation with vanishing boundaries.

clear;
close all;
clc;

%%% Setup parameters %%%

T = [0,16];     % location of the 1st and last sample in the time domain
XI = [0.5,23];  % location of the 1st and last sample in the xi-domain

%%% Setup the signal %%%

A = sym(15);
d = sym(0.5);
exp_t0 = sym(3000);
q_fun = @(t) A*sech((t-log(exp_t0))/d).^2;

%%% Analytic expressions for the nonlinear Fourier transform %%%

delta = sqrt(sym(A)*d^2+1/4);
at = @(xi) 1/2 - 1i*xi*d + delta;
bt = @(xi) 1/2 - 1i*xi*d - delta;
ct = @(xi) 1 - 1i*xi*d;
cgam = @(z) gamma(sym(z));

R_fun = @(xi) exp_t0.^(-2i*xi) .* cgam(ct(xi)-at(xi)-bt(xi)) .* cgam(at(xi)) .* cgam(bt(xi)) ./ ...
                       ( cgam(ct(xi)-at(xi)) .* cgam(ct(xi)-bt(xi)) .* cgam(at(xi)+bt(xi)-ct(xi)) );

%%% Prepare variables to store errors and runtimes %%%

error_FCF2_1 = [];
error_FCF4_2 = [];
error_FCF_RE2_1 = [];
error_FCF_RE4_2 = [];
time_FCF2_1 = [];
time_FCF4_2 = [];
time_FCF_RE2_1 = [];
time_FCF_RE4_2 = [];

%%% Iterate of number of samples, gather errors and runtimes for each %%%

for D = 2.^(8:1:12)
    
    fprintf('Running codes with D=%d...',D);
    t = linspace(T(1),T(2),D);
    q = double(q_fun(t)); % signal samples
    R_exact = double(R_fun(linspace(XI(1),XI(2),D)));
    
    % Compute the continuous part of the nonlinear Fourier transform
    % numerically with different configurations, save correspoding errors
    % and runtimes.
    
    tic
    R_comp = mex_fnft_kdvv(q, T, XI, 'discr_2split4B', 'skip_bs');
    time_FCF2_1 = [time_FCF2_1,toc];
    error_FCF2_1 = [error_FCF2_1, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_kdvv(q, T, XI, 'discr_4split4B', 'skip_bs');
    time_FCF4_2 = [time_FCF4_2,toc];
    error_FCF4_2 = [error_FCF4_2, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_kdvv(q, T, XI, 'discr_2split4B', 'RE', 'skip_bs');
    time_FCF_RE2_1 = [time_FCF_RE2_1,toc];
    error_FCF_RE2_1 = [error_FCF_RE2_1, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_kdvv(q, T, XI, 'discr_4split4B', 'RE', 'skip_bs');
    time_FCF_RE4_2 = [time_FCF_RE4_2,toc];
    error_FCF_RE4_2 = [error_FCF_RE4_2, norm(R_exact-R_comp)/norm(R_exact)];
     
    fprintf('Done.\n');
end

%%% Plot results %%%

lw = 3;
fs = 15;
ms = 6;
alw = 1;
loglog(time_FCF2_1,error_FCF2_1,'-','linewidth',lw,'markersize',ms);                           
hold on
loglog(time_FCF4_2,error_FCF4_2,'--','linewidth',lw,'markersize',ms);
ax = gca;
ax.ColorOrderIndex = 1;
loglog(time_FCF_RE2_1,error_FCF_RE2_1,':','linewidth',lw,'markersize',ms);
loglog(time_FCF_RE4_2,error_FCF_RE4_2,'-.','linewidth',lw,'markersize',ms);
axis tight

set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Times','LineWidth',alw)

xlabel({'Execution Time (s)'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Times');
ylabel({'Relative $L^2$-error in $\rho(\lambda)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Times')
grid on
grid minor
le = legend("FCF$^{[2]}_1$","FCF$^{[4]}_2$","FCF\_RE$^{[2]}_1$","FCF\_RE$^{[4]}_2$");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Times','Location','SouthWest');