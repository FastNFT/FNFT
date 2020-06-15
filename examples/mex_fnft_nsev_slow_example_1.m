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
% Shinivas Chimmalgi(TU Delft) 2020.

% This example compares the accuracy the various slow algorithms for
% the defocusing nonlinear Schroedinger equation with vanishing boundaries.
% The example is taken from Section 5.3 of
% https://doi.org/10.1364/OE.377140. Plots similar to those in Figs. 2 are
% shown.



clear;
close all;
clc;

%%% Setup parameters %%%

T = [-30,30];  % location of the 1st and last sample in the time domain
XI = [-20,20];  % location of the 1st and last sample in the xi-domain
kappa = -1;     % defocusing nonlinear Schroedinger equation
M = 1025;
xi = linspace(XI(1),XI(2),M);

%%% Setup the signal %%%
A = 5.2;
C = 4;
q_fun = @(t) A*(sech(t)).^(1+1i*C);

D = sqrt(kappa*A^2-C^2/4);
b = double((1/(A*2^(1i*C)))*gamma(sym(0.5-1i*(xi+C/2))).*gamma(sym(0.5+1i*(xi-C/2)))./(gamma(sym(-1i*C/2-D)).*gamma(sym(-1i*C/2+D))));
a = double(gamma(sym(0.5-1i*(xi+C/2))).*gamma(sym(0.5-1i*(xi-C/2)))./(gamma(sym(0.5-1i*xi-D)).*gamma(sym(0.5-1i*xi+D))));


errorBO_a = [];
errorCF4_2_a= [];
errorCF4_3_a = [];
errorCF5_3_a= [];
errorCF6_4_a= [];

errorBO_b = [];
errorCF4_2_b= [];
errorCF4_3_b = [];
errorCF5_3_b= [];
errorCF6_4_b= [];

timeBO = [];
timeCF4_2= [];
timeCF4_3 = [];
timeCF5_3= [];
timeCF6_4= [];

po = 7:11;
NT = 2.^po;

for N = NT
    N
    t=linspace(T(1),T(2),2*N+1);
    
    q = q_fun(t);
    tic
    [ab]=mex_fnft_nsev_slow(q, T, XI, kappa,'M',M,'BO','cstype_ab');
    timeBO= [timeBO;toc];
    errorBO_a= [errorBO_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorBO_b= [errorBO_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_nsev_slow(q, T, XI, kappa,'M',M,'CF4_2','cstype_ab');
    timeCF4_2= [timeCF4_2;toc];
    errorCF4_2_a= [errorCF4_2_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF4_2_b= [errorCF4_2_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_nsev_slow(q, T, XI, kappa,'M',M,'CF4_3','cstype_ab');
    timeCF4_3= [timeCF4_3;toc];
    errorCF4_3_a= [errorCF4_3_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF4_3_b= [errorCF4_3_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_nsev_slow(q, T, XI, kappa,'M',M,'CF5_3','cstype_ab');
    timeCF5_3= [timeCF5_3;toc];
    errorCF5_3_a= [errorCF5_3_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF5_3_b= [errorCF5_3_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_nsev_slow(q, T, XI, kappa,'M',M,'CF6_4','cstype_ab');
    timeCF6_4= [timeCF6_4;toc];
    errorCF6_4_a= [errorCF6_4_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF6_4_b= [errorCF6_4_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
end
%% Plotting results

lw = 2;
fs = 14;
ms = 6;
alw = 1;
figure
semilogy(po,errorBO_a,'o-','Color',[0.18,0.21,0.59],'MarkerFaceColor',[0.18,0.21,0.59],'linewidth',lw,'markersize',ms);
hold on
semilogy(po,errorCF4_2_a,'d-','Color',[0.93,0.11,0.16],'MarkerFaceColor',[0.93,0.11,0.16],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF4_3_a,'v-','MarkerFaceColor',[0.14,0.12,0.13],'Color',[0.14,0.12,0.13],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF5_3_a,'s-','MarkerFaceColor',[0.22,0.21,0.59],'Color',[0.22,0.21,0.59],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF6_4_a,'>-','MarkerFaceColor',[0.63 0.08 0.18],'Color',[0.63 0.08 0.18],'linewidth',lw,'markersize',ms);

ylim([1e-26,1e2])
set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'NMSE[$a(\xi)$]'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(D)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
grid on
grid minor
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','NorthEast','NumColumns',2);

figure
semilogy(po,errorBO_b,'o-','Color',[0.18,0.21,0.59],'MarkerFaceColor',[0.18,0.21,0.59],'linewidth',lw,'markersize',ms);
hold on
semilogy(po,errorCF4_2_b,'d-','Color',[0.93,0.11,0.16],'MarkerFaceColor',[0.93,0.11,0.16],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF4_3_b,'v-','MarkerFaceColor',[0.14,0.12,0.13],'Color',[0.14,0.12,0.13],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF5_3_b,'s-','MarkerFaceColor',[0.22,0.21,0.59],'Color',[0.22,0.21,0.59],'linewidth',lw,'markersize',ms);
semilogy(po,errorCF6_4_b,'>-','MarkerFaceColor',[0.63 0.08 0.18],'Color',[0.63 0.08 0.18],'linewidth',lw,'markersize',ms);

ylim([1e-26,1e2])
grid on
grid minor
set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'NMSE[$b(\xi)$]'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(D)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','NorthEast','NumColumns',2);

figure
semilogy(po,timeBO,'o-','Color',[0.18,0.21,0.59],'MarkerFaceColor',[0.18,0.21,0.59],'linewidth',lw,'markersize',ms);
hold on
semilogy(po,timeCF4_2,'d-','Color',[0.93,0.11,0.16],'MarkerFaceColor',[0.93,0.11,0.16],'linewidth',lw,'markersize',ms);
semilogy(po,timeCF4_3,'v-','MarkerFaceColor',[0.14,0.12,0.13],'Color',[0.14,0.12,0.13],'linewidth',lw,'markersize',ms);
semilogy(po,timeCF5_3,'s-','MarkerFaceColor',[0.22,0.21,0.59],'Color',[0.22,0.21,0.59],'linewidth',lw,'markersize',ms);
semilogy(po,timeCF6_4,'>-','MarkerFaceColor',[0.63 0.08 0.18],'Color',[0.63 0.08 0.18],'linewidth',lw,'markersize',ms);

ylim([0,1e2])
grid on
grid minor
set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'Execution time'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(D)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','NorthEast','NumColumns',2);

