clear
close all
load('discr_BO_sech_changing_M')     % discr_<method>_<signal>
D_BO_sech_changing_M = Test_results.params.D_values;
exec_time_BO_sech_changing_M = zeros(1,length(D_BO_sech_changing_M));
for i=1:length(Test_results.params.D_values)
    exec_time_BO_sech_changing_M(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

load('discr_2split3A_sech_changing_M')     % discr_<method>_<signal>
D_2split3A_sech_changing_M = Test_results.params.D_values;
exec_time_2split3A_sech_changing_M = zeros(1,length(D_2split3A_sech_changing_M));
for i=1:length(Test_results.params.D_values)
    exec_time_2split3A_sech_changing_M(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

load('discr_BO_sech')     % discr_<method>_<signal>
D_BO_sech = Test_results.params.D_values;
exec_time_BO_sech = zeros(1,length(D_BO_sech));
for i=1:length(Test_results.params.D_values)
    exec_time_BO_sech(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

load('discr_2split3A_sech')     % discr_<method>_<signal>
D_2split3A_sech = Test_results.params.D_values;
exec_time_2split3A_sech = zeros(1,length(D_2split3A_sech));
for i=1:length(Test_results.params.D_values)
    exec_time_2split3A_sech(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

load('discr_2split3A_sech_D_powers_of_2')     % discr_<method>_<signal>
D_2split3A_sech_D_powers_of_2 = Test_results.params.D_values;
exec_time_2split3A_sech_D_powers_of_2 = zeros(1,length(D_2split3A_sech_D_powers_of_2));
for i=1:length(Test_results.params.D_values)
    exec_time_2split3A_sech_D_powers_of_2(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

load('discr_2split3A_build_polys')     % discr_<method>_<signal>
D_2split3A_build_polys = Test_results.params.D_values;
exec_time_2split3A_build_polys = zeros(1,length(D_2split3A_build_polys));
for i=1:length(Test_results.params.D_values)
    exec_time_2split3A_build_polys(i) = Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).runtime;
end

%%
figure;
loglog(D_2split3A_sech_changing_M, exec_time_2split3A_sech_changing_M,'*-')
hold on
loglog(D_BO_sech_changing_M, exec_time_BO_sech_changing_M,'*-')
grid on
title('execution time, M=D');
xlabel('number of samples');
ylabel('time [s]');
legend('2split3A','BO','Location','northwest')

figure;
loglog(D_2split3A_sech, exec_time_2split3A_sech,'*-')
hold on
loglog(D_BO_sech, exec_time_BO_sech,'*-')
grid on
title('execution time, M=100 (const)');
xlabel('number of samples');
ylabel('time [s]');
legend('2split3A','BO','Location','northwest')

figure;
loglog(D_2split3A_sech, exec_time_2split3A_sech,'b*-')
hold on
loglog(D_2split3A_sech_changing_M, exec_time_2split3A_sech_changing_M,'bo-.')
loglog(D_2split3A_build_polys, exec_time_2split3A_build_polys,'m+-')
loglog(D_2split3A_sech_D_powers_of_2, exec_time_2split3A_sech_D_powers_of_2,'bx-')

loglog(D_BO_sech, exec_time_BO_sech,'r*-')
loglog(D_BO_sech_changing_M, exec_time_BO_sech_changing_M,'ro-.')
loglog(D_2split3A_sech, 0.0001*D_2split3A_sech.*(log(D_2split3A_sech)./log(2)).^2,'k')
loglog(D_2split3A_build_polys, 0.0001*D_2split3A_build_polys,'k--')
grid on
title('execution time with M=D and M const');
xlabel('number of samples');
ylabel('time [s]');
legend('2split3A, M=100', '2split3A with changing M', '2split3A, building polynomials', '2split3A, M=100, D powers of 2','BO, M=100', 'BO with changing M', 'Slope of D*log^2(D)', 'Slope of D', 'Location','northwest')





