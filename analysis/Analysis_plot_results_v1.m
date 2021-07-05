%% This script plots the results of a single method
clear
close all
load('discr_FTES4_suzuki_sech')     % discr_<method>_<signal>

% Plotting q1, q2 (for the last number of samples D in D_values)
t = Test_results.params.t;
q1 = Test_results.params.q1;
q2 = Test_results.params.q2;

figure;
plot(t, real(q1), t, imag(q1));
title('Time-domain');
xlabel('t');
ylabel('q1(t)');
legend('Real part', 'Imaginary part');

figure;
plot(t, real(q2), t, imag(q2));
title('Time-domain');
xlabel('t');
ylabel('q2(t)');
legend('Real part', 'Imaginary part');

%% Plotting the NFT coeffs and rho for multiple D's
xi = linspace(Test_results.params.XI(1),Test_results.params.XI(2),Test_results.params.M);
figure;
hold on
plot(xi, real(Test_results.exact_sol.rho1));
for i=1:length(Test_results.params.D_values)
plot(xi, real(Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).rho1));
end
title('Continuous spectrum, first element real part');
xlabel('\xi');
ylabel('r(\xi)');
legend('Exact solution');
for i=1:length(Test_results.params.D_values)
old_legend=findobj(gcf, 'Type', 'Legend');
legend([old_legend.String,(strcat('D = ',num2str(Test_results.params.D_values(i))))])
end

%% Plotting the RMSE:
RMSE = zeros(1,length(Test_results.params.D_values));

for i = 1:length(RMSE)
    RMSE(i) = norm(Test_results.(strcat('D_',num2str(Test_results.params.D_values(i)))).a - ...
        Test_results.exact_sol.a);
end

D_sizes = Test_results.params.D_values;

figure;
hold on
loglog(D_sizes, RMSE,'*-');
title(strcat('RMSE of the a coefficient, ',Test_results.params.method));
xlabel('number of samples');
ylabel('RMSE');

%% determining the order of a method
errors_divided = zeros(length(D_sizes)-1,1);
D_divided = zeros(length(D_sizes)-1,1);

for index_order = 1:length(D_sizes)-1
    errors_divided(index_order) = RMSE(index_order)/RMSE(index_order+1);
    D_divided(index_order) = D_sizes(index_order+1)/D_sizes(index_order);
end

order = log(errors_divided)./log(D_divided)
