%% Run tests and collect results for nsev

methods = ["discr_2split3B",
           "discr_2split4A",
           "discr_2split4B",
           "discr_2split6B",
           "discr_4split4A",
           "discr_4split4B"];
           
%methods = flip(methods);        

for method_index = 1:7
    
clearvars -except method_index methods;
close all;

%%% Setup parameters %%%
% signal and method
discretization = char(methods(method_index))    % see help mex_fnft_manakov for list of supported discretizations
signal = 'sech';    % potential function: sech or rect

% xi, t, and dt
T   = [-38.5, 38.5];      % location of the 1st and last sample in the time domain
% D_values = [50 60 70 80 90 128 300 512 1000 1500 2000 2500 4000 6500 16000];        % number of samples
D_values = zeros(1,12);
for k=1:16
    D_values(k) = 2^(k+1);
end
% D_values = (250:250:4000);
% D_values = 1000;
XI = [-pi, pi];   % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation
M=100;          % number of lambda for which we determine the spectrum

% signal parameters rectangle and sech
L = [-2,2];     % support of rectangle potential
A = 0.8;     % amplitude of signal for sech and rectangle

XI_vector = linspace(XI(1),XI(2),M);
Test_results.params.signal = signal;
Test_results.params.T = T;
Test_results.params.kappa = kappa;
Test_results.params.method = discretization;
Test_results.params.M = M;
Test_results.params.XI = XI;
Test_results.params.D_values = D_values;
if strcmp(signal,'sech')
    Test_results.params.ampl = [A];
elseif strcmp(signal,'rect')
    Test_results.params.ampl = [A];
    Test_results.params.support = L;
end  

%%% Setup the signal %%%
for i =1:length(D_values)
    D=D_values(i)
    
%    M = D;      % For checking if fast methods are asymptotically faster
%    XI_vector = linspace(XI(1),XI(2),M);

    eps_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1):eps_t:T(2);
if strcmp(signal,'sech')
    q1 = A*sech(t);               % signal samples
elseif strcmp(signal,'rect')
        q1 = rectangle_function(t, A, L);
else
    error('unknown test signal')
end        

%%% Compute the nonlinear Fourier transform %%%

tStart = tic;
[reflection_and_nft_coeffs] = mex_fnft_nsev(complex(q1), T, XI, kappa, 'M', M, 'cstype_ab', 'skip_bs', discretization);
% mex_fnft_nsev has many options => run "help mex_fnft_nsev" to learn more
tEnd = toc(tStart);

%%% Saving the data
Test_results.(strcat('D_',num2str(D))).a = reflection_and_nft_coeffs(1:M);
Test_results.(strcat('D_',num2str(D))).b = reflection_and_nft_coeffs(M+1:2*M);
Test_results.(strcat('D_',num2str(D))).runtime = tEnd;

% runtime_build_pol(i) = tEnd;

end

%%
% loglog(D_values,runtime_build_pol,'*')
% hold on
% loglog(D_values,0.00001*D_values)
% loglog(D_values,0.0000001*D_values.^2)
% %%
% %%% Getting the exact solution
% if strcmp(signal,'sech')
%     [a_exact, b1_exact, b2_exact] = Manakov_sech_exact(A1, A2, XI_vector, kappa);
% elseif strcmp(signal,'rect')
%     [a_exact, b1_exact, b2_exact] = Manakov_rectangle_exact(A1, A2, XI_vector, kappa, L);
% elseif strcmp(signal,'single_soliton')
%     [a_exact] = Manakov_soliton_exact(xi, eta, XI_vector);
%     Test_results.exact_sol.bound_states = xi+1i*eta;
% end
% Test_results.exact_sol.a = a_exact;
% 
% if ~strcmp(signal,'single_soliton')
% Test_results.exact_sol.b1 = b1_exact;
% Test_results.exact_sol.b2 = b2_exact;
% Test_results.exact_sol.rho1 = b1_exact./a_exact;
% Test_results.exact_sol.rho2 = b2_exact./a_exact;
% end

Test_results.params.t = t;
Test_results.params.q1 = q1;

% save(strcat(discretization,'_',signal,'get_pol_coeffs'),'Test_results')
% save(strcat(discretization,'_',signal,'_D_powers_of_2'),'Test_results')
% save(strcat(discretization,'_',signal,'_powers_2_v2'),'Test_results')

save(strcat(discretization,'_',signal,'_27jul_nse_test_without_fast_size'),'Test_results')

end
%% Auxiliary functions
function y = rectangle_function(t, A, L)
y=zeros(1,length(t));
for i=1:length(t)
    if (t(i)>=L(1) && t(i)<=L(2))
        y(i) = A;
    else
        y(i) = 0;
    end
end
end

