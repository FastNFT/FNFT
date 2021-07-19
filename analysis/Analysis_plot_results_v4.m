%% This script plots the execution time against the number of samples
clear
close all

ext = "sech_9jul_v2";       % what follows after discr_method_

methods = ["discr_2split3A",
           "discr_2split3B",
           "discr_2split4A",
           "discr_2split4B",
           "discr_2split6B",
           "discr_BO"];
       
methods = ["discr_4split4A",
           "discr_4split4B",
           "discr_4split6B",
           "discr_FTES4_suzuki",
           "discr_CF4_2"];
%        
% methods = ["discr_2split4B",
%            "discr_2split6B",
%            "discr_4split4B",
%            "discr_4split6B",
%            "discr_FTES4_suzuki",
%            "discr_BO",
%            "discr_CF4_2"];
       
% methods = ["discr_2split4B",
%            "discr_2split6B",
%            "discr_BO"];
       
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
    for i=1:length(D)
            runtime(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
    end
    loglog(D, runtime,'*-','LineWidth',1.5,'MarkerSize',12,'color',line_colour(method_index,:))
    hold on
    clear Test_results
end

% % Plotting some reference lines
% loglog(D, 0.01*D.*(log(D)./log(2)).^2,'k','LineWidth',1.5)
% loglog(D, 0.0000001*D.^2,'k--','LineWidth',1.5)

% Legend, labeling and title
legend(names(1),'Location','northwest')
for method_index = 2:length(methods)
    old_legend=findobj(gcf, 'Type', 'Legend');
    legend([old_legend.String,names(method_index)])
end
xlabel('number of samples');
ylabel('runtime (s)');

old_legend=findobj(gcf, 'Type', 'Legend');
legend([old_legend.String,"order D/log^2(D)"])
old_legend=findobj(gcf, 'Type', 'Legend');
legend([old_legend.String,"2nd order slope"])