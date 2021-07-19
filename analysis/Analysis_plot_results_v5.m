%% This script plots the error against runtime
clear
close all

ext = "sech15jul_v1";       % what follows after discr_method_
var = 'a';

methods = ["discr_4split4B",
           "discr_4split6B",
           "discr_FTES4_suzuki",
           "discr_CF4_2"];

% methods = ["discr_2split4B",
%            "discr_2split6B",
%            "discr_BO"]

% methods = ["discr_2split4B",
%            "discr_2split6B",
%            "discr_4split4B",
%            "discr_4split6B",
%            "discr_FTES4_suzuki",
%            "discr_BO",
%            "discr_CF4_2"];


for ind = 1:length(methods)
    if strcmp(methods(ind),'discr_2split3A')
        names(ind) = "2split3A";
        line_colour(ind,:) = [1 0 0];
    elseif strcmp(methods(ind),'discr_2split3B')
        names(ind) = "2split3B";
        line_colour(ind,:) = [1 0 1];
    elseif strcmp(methods(ind),'discr_2split4A')
        names(ind) = "2split4A";
        line_colour(ind,:) = [0 0 1];
    elseif strcmp(methods(ind),'discr_2split4B')
        names(ind) = "2split4B";
        line_colour(ind,:) = [0 1 0];
    elseif strcmp(methods(ind),'discr_2split6B')
        names(ind) = "2split6B";
        line_colour(ind,:) = [1 1 0];
    elseif strcmp(methods(ind),'discr_4split4A')
        names(ind) = "4split4A";
        line_colour(ind,:) = [0.8500 0.3250 0.0980];
    elseif strcmp(methods(ind),'discr_4split4B')
        names(ind) = "4split4B";
        line_colour(ind,:) = [0 0.4470 0.7410];
    elseif strcmp(methods(ind),'discr_4split6B')
        names(ind) = "4split6B";
        line_colour(ind,:) = [0.6350 0.0780 0.1840];
    elseif strcmp(methods(ind),'discr_FTES4_suzuki')
        names(ind) = "FTES4\_suzuki";
        line_colour(ind,:) = [0.4660 0.6740 0.1880];
    elseif strcmp(methods(ind),'discr_BO')
        names(ind) = "BO";
        line_colour(ind,:) = [0 1 1];
    elseif strcmp(methods(ind),'discr_CF4_2')
        names(ind) = "CF_2^{[4]}";
        line_colour(ind,:) = [0.9290 0.6940 0.1250];
    elseif strcmp(methods(ind),'discr_RK4')
        names(ind) = "RK_4";
        line_colour(ind,:) = [0.3010 0.7450 0.9330];
    end
end
       
for method_index = 1:length(methods)
    load(strcat(methods(method_index),'_',ext,'.mat'))
    D= Test_results.params.D_values;
    runtime = zeros(1,length(D));
    RMSE = zeros(1,length(D));
    for i=1:length(D)
            runtime(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
            if strcmp(var,'a')
            RMSE(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
        elseif strcmp(var,'b1')
            RMSE(i) = norm(Test_results.exact_sol.b1-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).b1)/norm(Test_results.exact_sol.b1);
        elseif strcmp(var,'b2')
            RMSE(i) = norm(Test_results.exact_sol.b2-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).b2)/norm(Test_results.exact_sol.b2);
        elseif strcmp(var,'rho1')
            RMSE(i) = norm(Test_results.exact_sol.rho1-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).rho1)/norm(Test_results.exact_sol.rho1);
        elseif strcmp(var,'rho2')
            RMSE(i) = norm(Test_results.exact_sol.rho2-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).rho2)/norm(Test_results.exact_sol.rho2);
        end
    end
    loglog(runtime, RMSE,'*-','LineWidth',1.5,'MarkerSize',12,'color',line_colour(method_index,:))
    hold on
    clear Test_results
end

% Legend, labeling and title
legend(names(1),'Location','southwest')
for method_index = 2:length(methods)
    old_legend=findobj(gcf, 'Type', 'Legend');
    legend([old_legend.String,names(method_index)])
end

xlabel('runtime');
ylabel('normalized RMSE');

