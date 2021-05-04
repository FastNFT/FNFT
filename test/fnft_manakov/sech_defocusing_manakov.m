%% File to determine the exact NFT spectrum of focusing sech potential
% potential function: q = [A1*sech(t); A2*sech(t)]

A1 = sym(0.8);
A2 = sym(5.2);
kappa = sym(-1);
q0 = sqrt(kappa)*sqrt(A1*A1+A2*A2);
hlf = sym(1)/sym(2);
im = sym(1j);
syms lam;

a(lam) = (gamma(hlf-im*lam)^2)/(gamma(hlf-im*lam+q0)*gamma(hlf-im*lam-q0));
b1(lam) = -kappa*sin(pi*q0)*A1*sech(lam*pi)/q0;
b2(lam) = -kappa*sin(pi*q0)*A2*sech(lam*pi)/q0;

q1(lam) = b1/a;
q2(lam) = b2/a;

xi = (-sym(7):sym(8))/4;        % xi has 16 elements, M = 16
XI = [xi(1) xi(end)]
digits(40);
contspec = vpa([(b1(xi)./a(xi)) (b2(xi)./a(xi))]).'
ab = vpa([a(xi) b1(xi) b2(xi)]).'
        
        