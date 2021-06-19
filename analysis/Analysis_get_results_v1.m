%% Run tests and collect results to analyse performance of methods

clear all;
close all;

%%% Setup parameters %%%

signal = 'sech';    % potential function: sech or rect
T = [-15, 15];      % location of the 1st and last sample in the time domain
D_values = [50 60 70 80 90 128 300 512 1000 1500 2000 2500];        % number of samples
XI = [-7/4, 8/4];   % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation
L = [-2,2];     % support of rectangle potential
discretization = 'discr_2split3A';    % see help mex_fnft_manakov for list of supported discretizations
A1 = 0.8; A2 = 5.2;     % amplitudes of signal q1 and q2
M=100;          % number of lambda for which we determine the spectrum
XI_vector = linspace(XI(1),XI(2),M);

Test_results.params.signal = signal;
Test_results.params.T = T;
Test_results.params.kappa = kappa;
Test_results.params.ampl = [A1 A2];
Test_results.params.method = discretization;
Test_results.params.M = M;
Test_results.params.XI = XI;
Test_results.params.D_values = D_values;

%%% Setup the signal %%%
for i =1:length(D_values)
    D=D_values(i)
eps_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1):eps_t:T(2);
if strcmp(signal,'sech')
    q1 = A1*sech(t);               % signal samples
    q2 = A2*sech(t);
elseif strcmp(signal,'rect')
        q1 = rectangle_function(t, A1, L);
        q2 = rectangle_function(t, A2, L);
else
    error('unknown test signal')
end
    
        

%%% Compute the nonlinear Fourier transform %%%

tStart = tic;
[reflection_and_nft_coeffs] = mex_fnft_manakov(complex(q1), complex(q2), T, XI, kappa, 'M', M, 'cstype_both', discretization);
% mex_fnft_manakov has many options => run "help mex_fnft_manakov" to learn more
tEnd = toc(tStart);

%%% Saving the data
Test_results.(strcat('D_',num2str(D))).rho1 = reflection_and_nft_coeffs(1:M);
Test_results.(strcat('D_',num2str(D))).rho2 = reflection_and_nft_coeffs(M+1:2*M);
Test_results.(strcat('D_',num2str(D))).a = reflection_and_nft_coeffs(2*M+1:3*M);
Test_results.(strcat('D_',num2str(D))).b1 = reflection_and_nft_coeffs(3*M+1:4*M);
Test_results.(strcat('D_',num2str(D))).b2 = reflection_and_nft_coeffs(4*M+1:5*M);
Test_results.(strcat('D_',num2str(D))).runtime = tEnd;

end
%%% Getting the exact solution
if strcmp(signal,'sech')
    [a_exact, b1_exact, b2_exact] = Manakov_sech_exact(A1, A2, XI_vector, kappa);
elseif strcmp(signal,'rect')
    [a_exact, b1_exact, b2_exact] = Manakov_rectangle_exact(A1, A2, XI_vector, kappa, L);
end
Test_results.exact_sol.a = a_exact;
Test_results.exact_sol.b1 = b1_exact;
Test_results.exact_sol.b2 = b2_exact;
Test_results.exact_sol.rho1 = b1_exact./a_exact;
Test_results.exact_sol.rho2 = b2_exact./a_exact;

Test_results.params.t = t;
Test_results.params.q1 = q1;
Test_results.params.q2 = q2;

save(strcat(discretization,'_',signal),'Test_results')



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

