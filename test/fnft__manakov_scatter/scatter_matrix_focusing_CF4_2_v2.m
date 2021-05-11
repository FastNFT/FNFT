clear
eps_t = 0.13;
    kappa = +1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact = []; lam = [2, 1+0.5*i];
    
    %% getting some values needed for the num. method
a1  = 0.25 + sqrt(3)/6;
a2  = 0.25 - sqrt(3)/6;

c1 = 0.5-sqrt(3)/6;
c2 = 0.5+sqrt(3)/6;

% getting the non-equidistant samples
[q1_c1, q2_c1] =  bandlimited_interpolation_CF24(eps_t,...
                                    [q1; q2],c1*eps_t);
[q1_c2, q2_c2] = bandlimited_interpolation_CF24(eps_t,...
                                    [q1; q2],c2*eps_t);

% Getting the interlaced "samples".
Neff = length(q1_c1)+length(q1_c2);
qeff = zeros(2,Neff);
reff = zeros(2,Neff);

qeff(1,2:2:end)=a2*q1_c1+a1*q1_c2;
qeff(2,2:2:end)=a2*q2_c1+a1*q2_c2;
reff(1,2:2:end)=-kappa*(a2*conj(q1_c1)+a1*conj(q1_c2));
reff(2,2:2:end)=-kappa*(a2*conj(q2_c1)+a1*conj(q2_c2));

qeff(1,1:2:end)=a1*q1_c1+a2*q1_c2;
qeff(2,1:2:end)=a1*q2_c1+a2*q2_c2;
reff(1,1:2:end)=-kappa*(a1*conj(q1_c1)+a2*conj(q1_c2));
reff(2,1:2:end)=-kappa*(a1*conj(q2_c1)+a2*conj(q2_c2));
    
    for i = 1:length(lam)
        lambda = lam(i)/2;
        S=eye(3);
        for n=1:2*D         % 2*D because we are effectively doing BO with twice the number of samples
                P = [-1i*lambda,           qeff(1,n),     qeff(2,n);
                reff(1,n),   1i*lambda, 0;
                reff(2,n),   0,      1i*lambda];

            TM = expm(eps_t*P);
            S=TM*S;
        end
        result_exact = [result_exact,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact.'
    
    
    function [q1s, q2s] = bandlimited_interpolation_CF24(eps_t, qn, ts)
% implements bandlimited interpolation (interpolation using the time shift
% property of the Fourier Transform)
% 
% inputs:
% eps_t timestep size
% qn    values of the samples of q taken at the points in vector tn
% ts    desired shift of the samples
% 
% output: qs    the shifted vector of samples

Qn = [fft(qn(1,:)); fft(qn(2,:))];
N = length(qn(1,:));
Np = floor(N/2);
Nn = -floor((N-1)/2);
Qn = [Qn(1,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t));...
      Qn(2,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t))];
q1s = ifft(Qn(1,:));
q2s = ifft(Qn(2,:));
end
