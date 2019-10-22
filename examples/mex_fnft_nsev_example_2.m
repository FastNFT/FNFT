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
% Shinivas Chimmalgi(TU Delft) 2019.

% This example compares the accuracy-execution tradeoff of the various algorithms for
% the focusing nonlinear Schroedinger equation with vanishing boundaries.
% It ties to recreate Fig. 13 from doi: 10.1109/ACCESS.2019.2945480.
% NOTE: The results do not exactly match the results in the paper due to slight
% differences in implementations. The overall trends are the same.

clear;
close all;
clc;

%%% Setup parameters %%%

T = [-32,32];  % location of the 1st and last sample in the time domain
XI = [-10,10];  % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation

%%% Setup the signal %%%
qo =  5.4;
lam0 = 3;
q_fun = @(t) qo*sech(t).*exp(-2i*t*lam0); %signal function
R_fun = @(xi) double((-sin(pi*qo)./cosh(pi*(xi-lam0))).*gamma(sym(0.5-1i*(xi-lam0)+qo)).*gamma(sym(0.5-1i*(xi-lam0)-qo))./(gamma(sym(0.5-1i*(xi-lam0)))).^2);
%b_fun = @(xi) double((-sin(pi*qo)./cosh(pi*(xi-lam0))));
%a_fun = @(xi) double(((gamma(sym(0.5-1i*(xi-lam0)))).^2)./(gamma(sym(0.5-1i*(xi-lam0)+qo)).*gamma(sym(0.5-1i*(xi-lam0)-qo))));

error_FCF2_1 = [];
error_FCF4_2 = [];
error_FCF_RE2_1 = [];
error_FCF_RE4_2 = [];
time_FCF2_1 = [];
time_FCF4_2 = [];
time_FCF_RE2_1 = [];
time_FCF_RE4_2 = [];

for D = 2.^(10:1:14)
    t = linspace(T(1),T(2),D);
    q = q_fun(t); %signal samples
    R_exact = R_fun(linspace(XI(1),XI(2),D));
    
    tic
    R_comp = mex_fnft_nsev(q, T, XI, kappa, 'discr_2split4B', 'skip_bs');
    time_FCF2_1 = [time_FCF2_1,toc];
    error_FCF2_1 = [error_FCF2_1, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_nsev(q, T, XI, kappa, 'discr_4split4B', 'skip_bs');
    time_FCF4_2 = [time_FCF4_2,toc];
    error_FCF4_2 = [error_FCF4_2, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_nsev(q, T, XI, kappa, 'discr_2split4B', 'RE', 'skip_bs');
    time_FCF_RE2_1 = [time_FCF_RE2_1,toc];
    error_FCF_RE2_1 = [error_FCF_RE2_1, norm(R_exact-R_comp)/norm(R_exact)];
    
    tic
    R_comp = mex_fnft_nsev(q, T, XI, kappa, 'discr_4split4B', 'RE', 'skip_bs');
    time_FCF_RE4_2 = [time_FCF_RE4_2,toc];
    error_FCF_RE4_2 = [error_FCF_RE4_2, norm(R_exact-R_comp)/norm(R_exact)];
end

%% Plotting results
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

