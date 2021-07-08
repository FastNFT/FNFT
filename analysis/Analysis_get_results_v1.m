%% Run tests and collect results to analyse performance of methods

clear all;
close all;

%%% Setup parameters %%%
% signal and method
discretization = 'discr_2split4B';    % see help mex_fnft_manakov for list of supported discretizations
signal = 'single_soliton';    % potential function: sech or rect

% xi, t, and dt
T = [-15, 15];      % location of the 1st and last sample in the time domain
% D_values = [50 60 70 80 90 128 300 512 1000 1500 2000 2500 4000 6500 16000];        % number of samples
D_values = zeros(1,8);
for k=1:8
    D_values(k) = 2^k;
end
D_values = 2^14;
XI = [-7/4, 8/4];   % location of the 1st and last sample in the xi-domain
kappa = +1;     % focusing nonlinear Schroedinger equation
M=100;          % number of lambda for which we determine the spectrum

% signal parameters rectangle and sech
L = [-2,2];     % support of rectangle potential
A1 = 0.8; A2 = 5.2;     % amplitudes of signal q1 and q2 for sech and rectangle
% signal parameters soliton
xi = 4.869;     % parameter defining the velocity of the soliton. real number
eta = 0.56;     % parameter defining the amplitude of the soliton. real number
x = 0.1;          % x coordinate at which we get the potential q(x,t) = q(t)
S = [6, 1+5i];  % vector defining the polarization of both components of q for single soliton



XI_vector = linspace(XI(1),XI(2),M);
Test_results.params.signal = signal;
Test_results.params.T = T;
Test_results.params.kappa = kappa;
Test_results.params.method = discretization;
Test_results.params.M = M;
Test_results.params.XI = XI;
Test_results.params.D_values = D_values;
if strcmp(signal,'sech')
    Test_results.params.ampl = [A1 A2];
elseif strcmp(signal,'rect')
    Test_results.params.ampl = [A1 A2];
    Test_results.params.support = L;
elseif strcmp(signal,'single_soliton')
    Test_results.params.xi = xi;
    Test_results.params.eta = eta;
    Test_results.params.x = x;
    Test_results.params.S = S;
end  

runtime_build_pol = zeros(1,length(D_values));

%%% Setup the signal %%%
for i =1:length(D_values)
    D=D_values(i)
    
%    M = D;      % For checking if fast methods are asymptotically faster
%    XI_vector = linspace(XI(1),XI(2),M);

    eps_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1):eps_t:T(2);
if strcmp(signal,'sech')
    q1 = A1*sech(t);               % signal samples
    q2 = A2*sech(t);
elseif strcmp(signal,'rect')
        q1 = rectangle_function(t, A1, L);
        q2 = rectangle_function(t, A2, L);
elseif strcmp(signal,'single_soliton')
    c = S/norm(S);      % normalize vector c for polarization
    x0 = log(norm(S)^2)/(4*eta);
    q1 = conj(-2*eta*sech(2*eta*(t-x0)+8*xi*eta*x).*exp(2*1i*xi*t+4i*(xi^2-eta^2)*x)*c(1));
    q2 = conj(-2*eta*sech(2*eta*(t-x0)+8*xi*eta*x).*exp(2*1i*xi*t+4i*(xi^2-eta^2)*x)*c(2));
else
    error('unknown test signal')
end        

%%% Compute the nonlinear Fourier transform %%%

tStart = tic;
[reflection_and_nft_coeffs, bound_states] = mex_fnft_manakovv(complex(q1), complex(q2), T, XI, kappa, 'M', M, 'cstype_both', discretization);
% mex_fnft_manakov has many options => run "help mex_fnft_manakov" to learn more
tEnd = toc(tStart);

%%% Saving the data
Test_results.(strcat('D_',num2str(D))).rho1 = reflection_and_nft_coeffs(1:M);
Test_results.(strcat('D_',num2str(D))).rho2 = reflection_and_nft_coeffs(M+1:2*M);
Test_results.(strcat('D_',num2str(D))).a = reflection_and_nft_coeffs(2*M+1:3*M);
Test_results.(strcat('D_',num2str(D))).b1 = reflection_and_nft_coeffs(3*M+1:4*M);
Test_results.(strcat('D_',num2str(D))).b2 = reflection_and_nft_coeffs(4*M+1:5*M);
Test_results.(strcat('D_',num2str(D))).bound_states = bound_states;
Test_results.(strcat('D_',num2str(D))).runtime = tEnd;

runtime_build_pol(i) = tEnd;

end

%%
loglog(D_values,runtime_build_pol,'*')
hold on
loglog(D_values,0.00001*D_values)
loglog(D_values,0.0000001*D_values.^2)
%%
%%% Getting the exact solution
if strcmp(signal,'sech')
    [a_exact, b1_exact, b2_exact] = Manakov_sech_exact(A1, A2, XI_vector, kappa);
elseif strcmp(signal,'rect')
    [a_exact, b1_exact, b2_exact] = Manakov_rectangle_exact(A1, A2, XI_vector, kappa, L);
elseif strcmp(signal,'single_soliton')
    [a_exact] = Manakov_soliton_exact(xi, eta, XI_vector);
    Test_results.exact_sol.bound_states = xi+1i*eta;
end
Test_results.exact_sol.a = a_exact;

if ~strcmp(signal,'single_soliton')
Test_results.exact_sol.b1 = b1_exact;
Test_results.exact_sol.b2 = b2_exact;
Test_results.exact_sol.rho1 = b1_exact./a_exact;
Test_results.exact_sol.rho2 = b2_exact./a_exact;
end

Test_results.params.t = t;
Test_results.params.q1 = q1;
Test_results.params.q2 = q2;

% save(strcat(discretization,'_',signal,'get_pol_coeffs'),'Test_results')
% save(strcat(discretization,'_',signal,'_D_powers_of_2'),'Test_results')
% save(strcat(discretization,'_',signal,'_powers_2_v2'),'Test_results')

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
