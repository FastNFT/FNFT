clear lam lambda result_exact
clear
eps_t = 0.13;
    kappa = -1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact = []; lam = [2, 1+0.5*i]./2;
    for i = 1:length(lam)
        lambda = lam(i);
        S=eye(3);
        for n=1:D
            P = [-1i*lambda,           q1(n),     q2(n);
                -kappa*conj(q1(n)),   1i*lambda, 0;
                -kappa*conj(q2(n)),   0,      1i*lambda];

            TM = expm(eps_t*P);
            TM_reshaped = reshape(TM.',[1,9]);
            S=TM*S;
            S_reshaped = reshape(S.',[1,9]);
        end
        result_exact = [result_exact,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact.'















% clear
% eps_t = 0.13;
%     kappa = -1; D=8;
%     q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
%     q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
%     result_exact = []; lam = [2, 1+0.5*i];
%     
%     %% getting some values needed for the num. method
% a1  = 0.25 + sqrt(3)/6;
% a2  = 0.25 - sqrt(3)/6;
% 
% c1 = 0.5-sqrt(3)/6;
% c2 = 0.5+sqrt(3)/6;
% 
% % getting the non-equidistant samples
% [q1_c1, q2_c1] =  bandlimited_interpolation_CF24(eps_t,...
%                                     [q1; q2],c1*eps_t);
% [q1_c2, q2_c2] = bandlimited_interpolation_CF24(eps_t,...
%                                     [q1; q2],c2*eps_t);
%                                 
%                                 
% % Values from resampling code in C to be used in the remainder of this
% % program. They are slightly different from the values calculated in this
% % matlab file, because the implementation of bandlimited interpolation is
% % slightly different in the fnft library code
% q1_c1 = [1.959445940890e-001+1.241004022547e-001j, -2.592131203174e-001+3.452237856788e-001j, -3.941834020570e-001+3.791444155415e-001j, -1.864333946816e-001+5.044094791949e-001j, 1.758048139303e-001+4.751459873301e-001j, 4.268725816352e-001+5.002284202529e-001j, 1.899710763545e-001+3.984580528229e-001j, -1.566153115089e-002+3.160613770910e-001j];
% q1_c2 = [-5.590625485913e-002+2.138239999513e-001j, -3.856966939342e-001+4.046327835275e-001j, -3.206394446118e-001+4.243784759145e-001j, 4.180433278990e-002+5.216076384561e-001j, 3.299466175909e-001+4.701354532991e-001j, 3.883736304468e-001+4.656136858304e-001j, -5.143984459071e-002+3.498070606383e-001j, 1.866592749703e-001+1.927728225495e-001j];
% q2_c1 = [7.080038365208e-002+1.346306842042e-001j, -1.040163599220e-001+5.069460117029e-001j, -2.390168143246e-001+5.950595276911e-001j, -6.580700396317e-002+8.019351848364e-001j, 6.022699906875e-002+9.021869491973e-001j, 2.561786336031e-001+9.814499044214e-001j, 6.766428684148e-002+1.083299632198e+000j, 2.384822439043e-002+9.045427290125e-001j];
% q2_c2 = [-6.142131204569e-002+2.931676403723e-001j, -1.704202360708e-001+5.979502372586e-001j, -2.004062366658e-001+6.970919258710e-001j, 5.401780295934e-002+8.692997389206e-001j, 1.411514459906e-001+9.661912518895e-001j, 2.359666842292e-001+9.996280110607e-001j, -5.907644665477e-002+1.129339523823e+000j, 1.300666476041e-001+3.573822940680e-001j];
% 
% % Getting the interlaced "samples".
% Neff = length(q1_c1)+length(q1_c2);
% qeff = zeros(2,Neff);
% reff = zeros(2,Neff);
% 
% qeff(1,2:2:end)=a2*q1_c1+a1*q1_c2;
% qeff(2,2:2:end)=a2*q2_c1+a1*q2_c2;
% reff(1,2:2:end)=-kappa*(a2*conj(q1_c1)+a1*conj(q1_c2));
% reff(2,2:2:end)=-kappa*(a2*conj(q2_c1)+a1*conj(q2_c2));
% 
% qeff(1,1:2:end)=a1*q1_c1+a2*q1_c2;
% qeff(2,1:2:end)=a1*q2_c1+a2*q2_c2;
% reff(1,1:2:end)=-kappa*(a1*conj(q1_c1)+a2*conj(q1_c2));
% reff(2,1:2:end)=-kappa*(a1*conj(q2_c1)+a2*conj(q2_c2));
%     
%     for i = 1:length(lam)
%         lambda = lam(i)/2;
%         S=eye(3);
%         for n=1:2*D         % 2*D because we are effectively doing BO with twice the number of samples
%                 P = [-1i*lambda,           qeff(1,n),     qeff(2,n);
%                 reff(1,n),   1i*lambda, 0;
%                 reff(2,n),   0,      1i*lambda];
% 
%             TM = expm(eps_t*P);
%             S=TM*S;
%         end
%         result_exact = [result_exact,...
%             [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
%     end
%     format long g; result_exact.'
%     
%     
%     function [q1s, q2s] = bandlimited_interpolation_CF24(eps_t, qn, ts)
% % implements bandlimited interpolation (interpolation using the time shift
% % property of the Fourier Transform)
% % 
% % inputs:
% % eps_t timestep size
% % qn    values of the samples of q taken at the points in vector tn
% % ts    desired shift of the samples
% % 
% % output: qs    the shifted vector of samples
% 
% Qn = [fft(qn(1,:)); fft(qn(2,:))];
% N = length(qn(1,:));
% Np = floor(N/2);
% Nn = -floor((N-1)/2);
% Qn = [Qn(1,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t));...
%       Qn(2,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t))];
% q1s = ifft(Qn(1,:));
% q2s = ifft(Qn(2,:));
% end
