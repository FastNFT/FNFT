%% This script plots the RMSE against the number of samples on a loglog scale
%% New format: exact solution is in a separate .mat file
clear
close all
figure

ext = "single_soliton_31aug_v1";       % what follows after discr_method_

% methods = ["discr_2split3A",
%            "discr_2split3B",
%            "discr_2split4A",
%            "discr_2split4B",
%            "discr_2split6B",
%            "discr_4split4A",
%            "discr_4split4B",
%            "discr_4split6B",
%            "discr_FTES4_suzuki",
%            "discr_BO",
%            "discr_CF4_2",
%            "discr_RK4"];

% set for which method we plot the bound states
methods = ["discr_2split3A"];

% indices of which D's to include. 1-14 includes all
D_first = 1;
D_last = 14;
       
%% Don't change values here

% load exact solution
load(strcat("exact_sol_",ext,'.mat'))
exact_sol = Test_results;

for method_index = 1:length(methods)
    load(strcat(methods(method_index),'_',ext,'.mat'))
    D_temp= Test_results.params.D_values(D_first:D_last);
    
    % if bound_states is empty, don't plot those and remove from D:
    D_ind = 1;
    for ind = 1:length(D_temp)
       if ~isempty(Test_results.(strcat("D_",num2str(D_temp(ind)))).bound_states)
          D(D_ind) = D_temp(ind);
          D_ind = D_ind+1;
       end
    end
%     
% for ind = 1:length(D)
%     if D(ind) == 4
%         names(ind) = "D = 4";
%         line_colour(ind,:) = [1 0 0];
%         line_marker(ind) = 'p';
%     elseif D(ind) == 8
%         names(ind) = "D = 8";
%         line_colour(ind,:) = [1 0 1];
%         line_marker(ind) = 'x';
%     elseif D(ind) == 16
%         names(ind) = "D = 16";
%         line_colour(ind,:) = [0 0 1];
%         line_marker(ind) = 's';
%     elseif D(ind) == 32
%         names(ind) = "D = 32";
%         line_colour(ind,:) = [0 1 0];
%         line_marker(ind) = 'o';
%     elseif D(ind) == 64
%         names(ind) = "D = 64";
%         line_colour(ind,:) = [1 1 0];
%         line_marker(ind) = '+';
%     elseif D(ind) == 128
%         names(ind) = "D = 128";
%         line_colour(ind,:) = [0.8500 0.3250 0.0980];
%         line_marker(ind) = '*';
%     elseif D(ind) == 256
%         names(ind) = "D = 256";
%         line_colour(ind,:) = [0 0.4470 0.7410];
%         line_marker(ind) = 'd';
%     elseif D(ind) == 512
%         names(ind) = "D = 512";
%         line_colour(ind,:) = [0.6350 0.0780 0.1840];
%         line_marker(ind) = 'h';
%     elseif D(ind) == 1024
%         names(ind) = "D = 1024";
%         line_colour(ind,:) = [0.4660 0.6740 0.1880];
%         line_marker(ind) = '.';
%     elseif D(ind) == 2048
%         names(ind) = "D = 2048";
%         line_colour(ind,:) = [0 1 1];
%         line_marker(ind) = 'v';
%     elseif D(ind) == 4096
%         names(ind) = "D = 4096";
%         line_colour(ind,:) = [0.9290 0.6940 0.1250];
%         line_marker(ind) = '^';
%     elseif D(ind) == 8192
%         names(ind) = "D = 8192";
%         line_colour(ind,:) = [0.3010 0.7450 0.9330];
%         line_marker(ind) = '<';
%     elseif D(ind) == 16384
%         names(ind) = "D = 16384";
%         line_colour(ind,:) = [0 0 0];
%         line_marker(ind) = '|';
%     elseif D(ind) == 32768
%         names(ind) = "D = 32768";
%         line_colour(ind,:) = [0.6 1 0.50];
%         line_marker(ind) = '_';
%     end
% end


for ind = 1:length(D)
    if D(ind) == 4
        names(ind) = "D = 4";
        line_colour(ind,:) = [1 0 0];
        line_marker(ind) = 'p';
    elseif D(ind) == 8
        names(ind) = "D = 8";
        line_colour(ind,:) = [1 0 1];
        line_marker(ind) = 'x';
    elseif D(ind) == 16
        names(ind) = "D = 16";
        line_colour(ind,:) = [0 0 1];
        line_marker(ind) = 's';
    elseif D(ind) == 32
        names(ind) = "D = 32";
        line_colour(ind,:) = [0 1 0];
        line_marker(ind) = 'o';
    elseif D(ind) == 64
        names(ind) = "D = 64";
        line_colour(ind,:) = [1 1 0];
        line_marker(ind) = '+';
    elseif D(ind) == 128
        names(ind) = "D = 128";
        line_colour(ind,:) = [0.8500 0.3250 0.0980];
        line_marker(ind) = '*';
    elseif D(ind) == 256
        names(ind) = "D = 256";
        line_colour(ind,:) = [0 0.4470 0.7410];
        line_marker(ind) = 'd';
    elseif D(ind) == 512
        names(ind) = "D = 512";
        line_colour(ind,:) = [0.6350 0.0780 0.1840];
        line_marker(ind) = 'h';
    elseif D(ind) == 1024
        names(ind) = "D = 1024";
        line_colour(ind,:) = [0.4660 0.6740 0.1880];
        line_marker(ind) = '.';
    elseif D(ind) == 2048
        names(ind) = "D = 2048";
        line_colour(ind,:) = [0 1 1];
        line_marker(ind) = 'v';
    elseif D(ind) == 4096
        names(ind) = "D = 4096";
        line_colour(ind,:) = [0.9290 0.6940 0.1250];
        line_marker(ind) = '^';
    elseif D(ind) == 8192
        names(ind) = "D = 8192";
        line_colour(ind,:) = [0.3010 0.7450 0.9330];
        line_marker(ind) = '<';
    elseif D(ind) == 16384
        names(ind) = "D = 16384";
        line_colour(ind,:) = [0 0 0];
        line_marker(ind) = '|';
    elseif D(ind) == 32768
        names(ind) = "D = 32768";
        line_colour(ind,:) = [0.6 1 0.50];
        line_marker(ind) = '_';
    end
end

line_colour = [linspace(0,1,length(D))',...    % red
               linspace(1,0,length(D))',...    % green
               linspace(0,0,length(D))'];      % blue



    
    for D_ind = 1:length(D)
        plot(Test_results.(strcat("D_",num2str(D(D_ind)))).bound_states,'color',line_colour(D_ind,:),'LineStyle','none','LineWidth',1.5,'MarkerSize',12,'Marker',line_marker(D_ind))
        hold on
    end
end

% Legend, labeling and title
lgnd = legend(names(1),'Location','southwest');
for method_index = 2:length(D)
    old_legend=findobj(gcf, 'Type', 'Legend');
    legend([old_legend.String, names(method_index)])
end
grid on
xlabel('Re(\lambda_d)')
ylabel('Im(\lambda_d)')
% title(methods(1))

set(lgnd,'color','none');

