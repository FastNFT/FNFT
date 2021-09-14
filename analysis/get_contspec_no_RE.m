%% Run tests and collect results to analyse runtime of methods
% Runs all methods at once, without computing bound states. Runs tests with
% M=D. Computes exact result in separate file

methods = ["discr_2split3A",
           "discr_2split3B",
           "discr_2split4A",
           "discr_2split4B",
           "discr_2split6B",
           "discr_4split4A",
           "discr_4split4B",
           "discr_4split6B",
           "discr_FTES4_suzuki",
           "discr_BO",
           "discr_CF4_2",
           "discr_RK4",
           "exact_sol"];
       
methods = ["discr_4split4B",
           "discr_4split6B",
           "discr_FTES4_suzuki"];
              
for method_index = 1:length(methods)

clearvars -except method_index methods;
close all;

%%% Setup parameters %%%
% signal and method
discretization = char(methods(method_index))    % see help mex_fnft_manakov for list of supported discretizations
signal = 'sech';    % potential function: sech or rect or soliton

% xi, t, and dt
% D_values = zeros(1,14);
% for k=1:14
%     D_values(k) = 2^(k+1);
% end
D_values = (100:100:8000);
kappa = 1;     % focusing nonlinear Schroedinger equation

% Setting location of the 1st and last sample in the time domain
if strcmp(signal, 'sech')
    T   = [-38.5, 38.5];
    XI = [-4.6, 4.6];   % location of the 1st and last sample in the xi-domain
end
if strcmp(signal, 'rect')
    T   = [-2, 2];
    XI = [-250, 250];
end
if strcmp(signal, 'single_soliton')
    T   = [-31.5, 32.5];
    XI = [-110, 120];
end
if strcmp(signal, 'two_soliton')        % todo: set T
    T   = [-31.5, 30];
    XI = [-350, 350];
end

avg_flag = 1;       % if == 1, perform fnft 3x and take compute average runtime

%% signal parameters
if strcmp(signal,'sech') || strcmp(signal,'rect')
    A1 = 0.8; A2 = 5.2;     % amplitudes of signal q1 and q2 for sech and rectangle
end
if strcmp(signal,'rect')
    L = [-2,2];     % support of rectangle potential
end
if strcmp(signal,'single_soliton')% signal parameters soliton
    xi = 4.87;     % parameter defining the velocity of the soliton. real number
    eta = 0.56;     % parameter defining the amplitude of the soliton. real number
    x = 0.1;          % x coordinate at which we get the potential q(x,t) = q(t)
    S = [6, 1+5i];  % vector defining the polarization of both components of q for single soliton
end
if strcmp(signal,'two_soliton') % signal parameters two-soliton
x = 0.1;
alpha = [1.5, 3.6];
eta = [0.56, 1.28];
xi = [4.87, 0.358];
u = [1i, 2; 6, 2.2+3i];
for j=1:2
    u(j,:) = u(j,:)/norm(u(j,:));
end
end

Test_results.params.signal = signal;
Test_results.params.T = T;
Test_results.params.kappa = kappa;
Test_results.params.method = discretization;
Test_results.params.M = 'M=D';
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
elseif strcmp(signal,'two_soliton')
    Test_results.params.x = x;
    Test_results.params.alpha = alpha;
    Test_results.params.eta = eta;
    Test_results.params.xi = xi;
    Test_results.params.u = u;
end  

%%% Setup the signal %%%
for i =1:length(D_values)
    D=D_values(i)
    
M = D;      % For checking if fast methods are asymptotically faster
XI_vector = linspace(XI(1),XI(2),M);

eps_t = (T(2) - T(1)) / (D - 1); % time domain step size
t = T(1)-0.5*eps_t:eps_t:T(2)-0.5*eps_t;
T_f = T-0.5*eps_t;

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
elseif strcmp(signal,'two_soliton')
    [q1,q2] = two_soliton_function(t, x, alpha(1), alpha(2), eta(1), eta(2), xi(1), xi(2), u(1,:), u(2,:));
else
    error('unknown test signal')
end        

%%% Compute the nonlinear Fourier transform %%%

if strcmp(discretization,'discr_RK4')
    tStart = tic;
    [a,b1,b2] = RK4_Manakov(q1, q2, t, XI_vector, kappa);
    tEnd = toc(tStart);     % runtime is just calculated here s.t. all .mat
                            % files have the same forms and data in them,
                            % should NOT be used to compare to other
                            % methods as this runs in matlab and the other
                            % methods use C routines
    reflection_and_nft_coeffs = [b1./a, b2./a, a, b1, b2];
elseif strcmp(discretization,'exact_sol')
    tStart = tic;
    if strcmp(signal,'sech')
        [a, b1, b2] = Manakov_sech_exact(A1, A2, XI_vector, kappa);
        reflection_and_nft_coeffs = [b1./a, b2./a, a, b1, b2];
    elseif strcmp(signal,'rect')
        [a, b1, b2] = Manakov_rectangle_exact(A1, A2, XI_vector, kappa, L);
        reflection_and_nft_coeffs = [b1./a, b2./a, a, b1, b2];
    elseif strcmp(signal,'single_soliton')
        [a] = Manakov_soliton_exact(xi, eta, XI_vector);
        reflection_and_nft_coeffs = [zeros(1,M), zeros(1,M), a, zeros(1,M), zeros(1,M)];
    elseif strcmp(signal,'two_soliton')
        [a] = Manakov_two_soliton_exact(xi(1), xi(2), eta(1), eta(2), XI_vector);
        reflection_and_nft_coeffs = [zeros(1,M), zeros(1,M), a, zeros(1,M), zeros(1,M)];
    end
    tEnd = toc(tStart);     % runtime is just calculated here s.t. all .mat
                            % files have the same forms and data in them,
                            % should NOT be used to compare to other
                            % methods
else
if (avg_flag ==1)
    tEnd_for_avg = zeros(1,3);
for ind = 1:3
    tStart = tic;
    [reflection_and_nft_coeffs] = mex_fnft_manakovv(complex(q1), complex(q2), T_f, XI, kappa, 'M', M, 'cstype_both', 'skip_bs', discretization);
    tEnd_for_avg(ind) = toc(tStart);
end
    tEnd = sum(tEnd_for_avg)/3;
else % not averaging over 3 runs
    tStart = tic;
    [reflection_and_nft_coeffs] = mex_fnft_manakovv(complex(q1), complex(q2), T_f, XI, kappa, 'M', M, 'cstype_both', 'skip_bs', discretization);
    % mex_fnft_manakov has many options => run "help mex_fnft_manakov" to learn more
    tEnd = toc(tStart);
end
end

%%% Saving the data
Test_results.(strcat('D_',num2str(D))).rho1 = reflection_and_nft_coeffs(1:M);
Test_results.(strcat('D_',num2str(D))).rho2 = reflection_and_nft_coeffs(M+1:2*M);
Test_results.(strcat('D_',num2str(D))).a = reflection_and_nft_coeffs(2*M+1:3*M);
Test_results.(strcat('D_',num2str(D))).b1 = reflection_and_nft_coeffs(3*M+1:4*M);
Test_results.(strcat('D_',num2str(D))).b2 = reflection_and_nft_coeffs(4*M+1:5*M);
Test_results.(strcat('D_',num2str(D))).runtime = tEnd;
if ~strcmp(discretization,'exact_sol') && ~strcmp(discretization,'discr_RK4')
Test_results.(strcat('D_',num2str(D))).all_runtimes = tEnd_for_avg;
end

end
Test_results.params.t = t;
Test_results.params.q1 = q1;
Test_results.params.q2 = q2;

% save(strcat(discretization,'_',signal,'get_pol_coeffs'),'Test_results')
% save(strcat(discretization,'_',signal,'_D_powers_of_2'),'Test_results')
% save(strcat(discretization,'_',signal,'_powers_2_v2'),'Test_results')

save(strcat(discretization,'_',signal,'_6sept_v1'),'Test_results')

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

function [q1, q2] = two_soliton_function(x, t, alpha1, alpha2, eta1, eta2, xi1, xi2, u1, u2)
% Calculated values for two soliton (values x and t swapped)
u1 = u1./norm(u1);
u2 = u2./norm(u2);
zeta1 = xi1 + 1i*eta1;
zeta2 = xi2 + 1i*eta2;
lambda1 = -[2*eta1*(u1*u1')/((zeta1-conj(zeta1))*(zeta1-conj(zeta1))),...
            2*eta1*(u1*u1')/((zeta1-conj(zeta2))*(zeta1-conj(zeta1)));...
            2*eta1*(u2*u1')/((zeta1-conj(zeta1))*(zeta1-conj(zeta2))),...
            2*eta1*(u2*u1')/((zeta1-conj(zeta2))*(zeta1-conj(zeta2)))];
lambda2 = -[2*eta2*(u1*u2')/((zeta2-conj(zeta1))*(zeta2-conj(zeta1))),...
            2*eta2*(u1*u2')/((zeta2-conj(zeta2))*(zeta2-conj(zeta1)));...
            2*eta2*(u2*u2')/((zeta2-conj(zeta1))*(zeta2-conj(zeta2))),...
            2*eta2*(u2*u2')/((zeta2-conj(zeta2))*(zeta2-conj(zeta2)))];
    tau1 = 2*eta1*(x+4*xi1*t);
    tau2 = 2*eta2*(x+4*xi2*t);
    theta1 = 2*xi1*x+4*(xi1^2-eta1^2)*t;
    theta2 = 2*xi2*x+4*(xi2^2-eta2^2)*t;
    det_U = (exp(tau1+alpha1)/(2*eta1)).*(exp(tau2+alpha2)/(2*eta2))+...
        (exp(tau1+alpha1)/(2*eta1)).*(lambda1(2,2).*exp(-(tau1+alpha1)+1i*(theta1-theta2))+lambda2(2,2)*exp(-(tau2+alpha2)+1i*(theta2-theta2)))+...
        (exp(tau2+alpha2)/(2*eta2)).*(lambda1(1,1).*exp(-(tau1+alpha1)+1i*(theta1-theta1))+lambda2(1,1)*exp(-(tau2+alpha2)+1i*(theta2-theta1)))+...
        (exp(-(tau1+alpha1)).*exp(-(tau2+alpha2))).*(det(lambda1)+det(lambda2)+det([lambda1(1,1), lambda1(1,2); lambda2(2,1), lambda2(2,2)]) + det([lambda2(1,1), lambda2(1,2); lambda1(2,1), lambda1(2,2)]));
    q1 = (2./det_U).*((exp(tau2+alpha2)./(2*eta2)+(lambda1(2,2)-lambda1(2,1)).*exp(-(tau1+alpha1)+1i*(theta1-theta2))+(lambda2(2,2)-lambda2(2,1)).*exp(-(tau2+alpha2)+1i*(theta2-theta2))).*exp(-1i*theta1)*u1(1)+...
                      (exp(tau1+alpha1)./(2*eta1)+(lambda1(1,1)-lambda1(1,2)).*exp(-(tau1+alpha1)+1i*(theta1-theta1))+(lambda2(1,1)-lambda2(1,2)).*exp(-(tau2+alpha2)+1i*(theta2-theta1))).*exp(-1i*theta2)*u2(1));
    q2 = (2./det_U).*((exp(tau2+alpha2)./(2*eta2)+(lambda1(2,2)-lambda1(2,1)).*exp(-(tau1+alpha1)+1i*(theta1-theta2))+(lambda2(2,2)-lambda2(2,1)).*exp(-(tau2+alpha2)+1i*(theta2-theta2))).*exp(-1i*theta1)*u1(2)+...
                      (exp(tau1+alpha1)./(2*eta1)+(lambda1(1,1)-lambda1(1,2)).*exp(-(tau1+alpha1)+1i*(theta1-theta1))+(lambda2(1,1)-lambda2(1,2)).*exp(-(tau2+alpha2)+1i*(theta2-theta1))).*exp(-1i*theta2)*u2(2));
end




%% sech test case defocusing
