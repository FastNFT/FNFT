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
% Shinivas Chimmalgi (TU Delft) 2020.
% Sander Wahls (TU Delft) 2020.
% Peter J. Prins (TU Delft) 2021.

% This example compares the accuracy the various slow algorithms for
% the Korteweg-de Vries equation with vanishing boundaries.

clear;
close all;
clc;

%%% Setup parameters %%%

T = [0,16];     % location of the 1st and last sample in the time domain
XI = [0.5,23];  % location of the 1st and last sample in the xi-domain
M = 1025;
xi = linspace(XI(1),XI(2),M);

%%% Setup the signal %%%

A = sym(15);
d = sym(0.5);
exp_t0 = sym(3000);
q_fun = @(t) A*sech((t-log(exp_t0))/d).^2;

%%% Compute the nonlinear Fourier transform analytically %%%

fprintf('Computing exact spectrum symbolically - please wait ...');                
delta = sqrt(sym(A)*d^2+1/4);
at = 1/2 - 1i*xi*d + delta;
bt = 1/2 - 1i*xi*d - delta;
ct = 1 - 1i*xi*d;
cgam = @(z) gamma(sym(z));
a = double(cgam(ct) .* cgam(at+bt-ct) ./ ( cgam(at) .* cgam(bt) ));
b = double(exp_t0.^(-2i*xi) .* cgam(ct) .* cgam(ct-at-bt) ./ ( cgam(ct-at) .* cgam(ct-bt) ));
fprintf('done\n'); 

%%% Prepare variables to store errors and runtimes %%%

errorBO_a = [];
errorCF4_2_a= [];
errorCF4_3_a = [];
errorCF5_3_a= [];
errorCF6_4_a= [];
errorES4_a= [];
errorTES4_a= [];

errorBO_b = [];
errorCF4_2_b= [];
errorCF4_3_b = [];
errorCF5_3_b= [];
errorCF6_4_b= [];
errorES4_b= [];
errorTES4_b= [];

timeBO = [];
timeCF4_2= [];
timeCF4_3 = [];
timeCF5_3= [];
timeCF6_4= [];
timeES4= [];
timeTES4= [];

%%% Iterate of number of samples, gather errors and runtimes for each %%%

po = 7:11;
DT = 2.^po;
for D = DT
    fprintf('Running codes with D=%d...',D);
    t=linspace(T(1),T(2),2*D+1);
    
    % Compute the continuous part of the nonlinear Fourier transform
    % numerically with different configurations, save correspoding errors
    % and runtimes.
    
    q = double(q_fun(t));
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_BO','cstype_ab');
    timeBO= [timeBO;toc];
    errorBO_a= [errorBO_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorBO_b= [errorBO_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_CF4_2','cstype_ab');
    timeCF4_2= [timeCF4_2;toc];
    errorCF4_2_a= [errorCF4_2_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF4_2_b= [errorCF4_2_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_CF4_3','cstype_ab');    
    timeCF4_3= [timeCF4_3;toc];
    errorCF4_3_a= [errorCF4_3_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF4_3_b= [errorCF4_3_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_CF5_3','cstype_ab');
    timeCF5_3= [timeCF5_3;toc];
    errorCF5_3_a= [errorCF5_3_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF5_3_b= [errorCF5_3_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_CF6_4','cstype_ab');
    timeCF6_4= [timeCF6_4;toc];
    errorCF6_4_a= [errorCF6_4_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorCF6_4_b= [errorCF6_4_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_ES4','cstype_ab');   
    timeES4= [timeES4;toc];
    errorES4_a= [errorES4_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorES4_b= [errorES4_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];    
    
    tic
    [ab]=mex_fnft_kdvv(q, T, XI,'M',M,'discr_TES4','cstype_ab');   
    timeTES4= [timeTES4;toc];
    errorTES4_a= [errorTES4_a;sum((abs(ab(1:M)-a).^2)./(max(1,abs(a)).^2))/1025];
    errorTES4_b= [errorTES4_b;sum((abs(ab(M+1:end)-b).^2)./(max(1,abs(b)).^2))/1025];    
    
    fprintf('Done.\n');

end

%%% Plot results %%%

lw = 2;
fs = 14;
ms = 6;
alw = 1;
figure
semilogy(po,errorBO_a,'o-','linewidth',lw,'markersize',ms);
hold on
semilogy(po,errorCF4_2_a,'d-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF4_3_a,'v-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF5_3_a,'s-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF6_4_a,'>-','linewidth',lw,'markersize',ms);
semilogy(po,errorES4_a,'*-','linewidth',lw,'markersize',ms);
semilogy(po,errorTES4_a,'h-','linewidth',lw,'markersize',ms);

set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'NMSE[$a(\xi)$]'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(N)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
grid on
grid minor
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$","ES4","TES4");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','SouthWest','NumColumns',3);

figure
semilogy(po,errorBO_b,'o-','linewidth',lw,'markersize',ms);
hold on
semilogy(po,errorCF4_2_b,'d-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF4_3_b,'v-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF5_3_b,'s-','linewidth',lw,'markersize',ms);
semilogy(po,errorCF6_4_b,'>-','linewidth',lw,'markersize',ms);
semilogy(po,errorES4_b,'*-','linewidth',lw,'markersize',ms);
semilogy(po,errorTES4_b,'h-','linewidth',lw,'markersize',ms);

grid on
grid minor
set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'NMSE[$b(\xi)$]'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(N)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$","ES4","TES4");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','SouthWest','NumColumns',3);

figure
semilogy(po,timeBO,'o-','linewidth',lw,'markersize',ms);
hold on
semilogy(po,timeCF4_2,'d-','linewidth',lw,'markersize',ms);
semilogy(po,timeCF4_3,'v-','linewidth',lw,'markersize',ms);
semilogy(po,timeCF5_3,'s-','linewidth',lw,'markersize',ms);
semilogy(po,timeCF6_4,'>-','linewidth',lw,'markersize',ms);
semilogy(po,timeES4,'*-','linewidth',lw,'markersize',ms);
semilogy(po,timeTES4,'h-','linewidth',lw,'markersize',ms);

grid on
grid minor
set(gca, 'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',fs,'FontName','Arial','LineWidth',alw)
ylabel({'Execution time (s)'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial')
xlabel({'$\log_2(N)$'},'Interpreter','latex','FontUnits','points','Fontsize',fs,'FontName','Arial');
le = legend("BO","CF$^{[4]}_2$","CF$^{[4]}_3$","CF$^{[5]}_3$","CF$^{[6]}_4$","ES4","TES4");
set(le,'interpreter','latex','FontUnits','points','Fontsize',fs-2,'FontName','Arial','Location','NorthWest','NumColumns',3);