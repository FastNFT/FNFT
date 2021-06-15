%% This script plots the RMSE of multiple methods
clear
close all

load('discr_BO_rect')     % discr_<method>_<signal>
D_BO_rect = Test_results.params.D_values;
RMSE_a_BO_rect = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_BO_rect(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_BO_sech')     % discr_<method>_<signal>
D_BO_sech = Test_results.params.D_values;
RMSE_a_BO_sech = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_BO_sech(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split3A_sech')     % discr_<method>_<signal>
D_2split3A_sech = Test_results.params.D_values;
RMSE_a_2split3A_sech = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split3A_sech(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split3A_rect')     % discr_<method>_<signal>
D_2split3A_rect = Test_results.params.D_values;
RMSE_a_2split3A_rect = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split3A_rect(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split4A_sech')     % discr_<method>_<signal>
D_2split4A_sech = Test_results.params.D_values;
RMSE_a_2split4A_sech = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split4A_sech(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split4A_rect')     % discr_<method>_<signal>
D_2split4A_rect = Test_results.params.D_values;
RMSE_a_2split4A_rect = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split4A_rect(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split6B_sech')     % discr_<method>_<signal>
D_2split6B_sech = Test_results.params.D_values;
RMSE_a_2split6B_sech = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split6B_sech(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results

load('discr_2split6B_rect')     % discr_<method>_<signal>
D_2split6B_rect = Test_results.params.D_values;
RMSE_a_2split6B_rect = zeros(1,length(Test_results.params.D_values));
for i=1:length(Test_results.params.D_values)
    RMSE_a_2split6B_rect(i) = norm(Test_results.exact_sol.a-Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a)/norm(Test_results.exact_sol.a);
end
clear Test_results


% Plotting the errors

figure;
loglog(D_2split3A_rect, RMSE_a_2split3A_rect,'*-')
hold on
loglog(D_2split4A_rect, RMSE_a_2split4A_rect,'*-')
loglog(D_2split6B_rect, RMSE_a_2split6B_rect,'*-')
loglog(D_BO_rect, RMSE_a_BO_rect,'*-')
grid on
title('RMSE a coefficient rectangle potential');
xlabel('number of samples');
ylabel('RMSE');

figure;
loglog(D_2split3A_sech, RMSE_a_2split3A_sech,'*-')
hold on
loglog(D_2split4A_sech, RMSE_a_2split4A_sech,'*-')
loglog(D_2split6B_sech, RMSE_a_2split6B_sech,'*-')
grid on
title('RMSE a coefficient sech potential');
xlabel('number of samples');
ylabel('RMSE');


%% determining the order of a method
RMSE = RMSE_a_BO_rect; % Change to determine order of a different coefficient
D_sizes = D_BO_rect;

% RMSE = RMSE_a_2split3A_rect; % Change to determine order of a different coefficient
% D_sizes = D_2split3A_rect;
% 
% RMSE = RMSE_a_2split4A_sech; % Change to determine order of a different coefficient
% D_sizes = D_2split4A_sech;
% 
RMSE = RMSE_a_2split6B_sech; % Change to determine order of a different coefficient
D_sizes = D_2split6B_sech;

RMSE = RMSE_a_BO_sech; % Change to determine order of a different coefficient
D_sizes = D_BO_sech;

errors_divided = zeros(length(D_sizes)-1,1);
D_divided = zeros(length(D_sizes)-1,1);

for index_order = 1:length(D_sizes)-1
    errors_divided(index_order) = RMSE(index_order)/RMSE(index_order+1);
    D_divided(index_order) = D_sizes(index_order+1)/D_sizes(index_order);
end

order = log(errors_divided)./log(D_divided)

