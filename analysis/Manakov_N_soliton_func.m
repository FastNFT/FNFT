function [q1,q2] = Manakov_N_soliton_func(xi,eta,alpha,u_vec,x_vec,t_vec)

% % Example use
% % t = 0;
% % x = linspace(-10,10,2^10);
% % N = 4;
% % xi = -5+10*rand(1,N);
% % eta = 5*rand(1,N);
% % alpha = -2+2*rand(1,N);
% % u = -5+5*rand(N,2);
% % for j=1:length(eta)
% %     u(j,:) = u(j,:)/norm(u(j,:));
% % end
% % [q1,q2] = Manakov_N_soliton(xi,eta,alpha,u,x,t);


N = length(xi);
if length(eta)~=N || length(alpha)~=N
    error("Number of xi, eta, alpha should be same.");
end
if size(u_vec,2)~=2 || size(u_vec,1)~=N
    error("u_vec should be a Nx2 matrix.");
end
if sum(imag(xi)~=0)>0
    error("xi should be real.");
end
if sum(imag(eta)~=0)>0 || sum(eta<=0)>0
    error("eta should be real and greater than 0.");
end
if sum(imag(alpha)~=0)>0
    error("alpha should be real.");
end
for i=1:N
    if abs(norm(u_vec(i,:),2)-1)>1e-15
        error("norm of u should be 1.");
    end
end
q1 = zeros(length(t_vec),length(x_vec));
q2 = zeros(length(t_vec),length(x_vec));

zeta = xi+1i*eta;
lam_jkl = @(j,k,l) -2*eta(l)*(u_vec(j,:)*u_vec(l,:)')/((zeta(l)-conj(zeta(k)))*(zeta(l)-conj(zeta(j))));
for i=1:length(t_vec)
    t = t_vec(i);
    for n=1:length(x_vec)
        x = x_vec(n);
        tau = 2*eta.*(x+4*xi*t);
        Theta = 2*xi*x+4*(xi.^2 - eta.^2)*t;
        U = zeros(N);
        for j=1:N
            for k=1:N
                for l=1:N
                    U(j,k) = U(j,k) + lam_jkl(j,k,l)*exp(-(tau(l)+alpha(l))+1i*(Theta(l)-Theta(j)));
                end
                if j==k
                    U(j,k) = U(j,k) + exp(tau(j)+alpha(j))/(2*eta(j));
                end
            end
        end
        UI = inv(U);
        qtmp = 0;
        for j=1:N
            for k=1:N
                qtmp = qtmp + UI(j,k)*exp(-1i*Theta(k))*u_vec(k,:);
            end
        end
        qtmp = 2*qtmp;
        q1(i,n) =  qtmp(1);
        q2(i,n) =  qtmp(2);
        
    end
end
end