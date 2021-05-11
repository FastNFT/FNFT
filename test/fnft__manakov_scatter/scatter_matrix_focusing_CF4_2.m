clear lam lambda result_exact_CF
clear
eps_t = 0.13;
    kappa = +1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact_CF = []; lam = [2, 1+0.5*i];
    
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
Cn1_log = [];
Cn2_log = [];
Cn_log = [];


    for i = 1:length(lam)
        lambda = lam(i);
        S=eye(3);
        for n=1:D
            % determining transfer matrix for this timestep
        A1 = [-1i*lambda      q1_c1(n)         q2_c1(n);
              -kappa*conj(q1_c1(n))  1i*lambda 0;
              -kappa*conj(q2_c1(n))  0                   1i*lambda];
        A2 = [-1i*lambda      q1_c2(n)         q2_c2(n);
              -kappa*conj(q1_c2(n))  1i*lambda 0;
              -kappa*conj(q2_c2(n))  0                   1i*lambda];
        % index_t+1 because of matlab indexing: q(T0) = q[1], q(T0+nh) =
        % q[1+n].
        % notation from article: C(T- + (n+c1)h,lambda) = A1
        Cn2 = a2*A1 + a1*A2;
        Cn1 = a1*A1 + a2*A2;
        Cn1_log = [Cn1_log, [Cn1]];
        Cn2_log = [Cn2_log, [Cn2]];
        Cn_log = [Cn_log, Cn1, Cn2];
        TM = expm(eps_t*Cn2)*expm(eps_t*Cn1);
            TM_reshaped = reshape(TM.',[1,9]);
            S=TM*S;
            S_reshaped = reshape(S.',[1,9]);
        end
        result_exact_CF = [result_exact_CF,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact_CF.'
    
    
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

% code from the paper s 21
Qn = [fft(qn(1,:)); fft(qn(2,:))];
N = length(qn(1,:));
Np = floor(N/2);
Nn = -floor((N-1)/2);
Qn = [Qn(1,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t));...
      Qn(2,:).*exp(2i*pi*[0:Np, Nn:-1]*ts/(N*eps_t))];
q1s = ifft(Qn(1,:));
q2s = ifft(Qn(2,:));
end
