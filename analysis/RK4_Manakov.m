function [a,b1,b2] = RK4_Manakov(q1, q2, t, xi, kappa)

N=length(t);

h = t(end)-t(end-1);
ts = t+h/2;
q1s = interp1(t,q1,t(1:end-1)+h/2,'spline');
q2s = interp1(t,q2,t(1:end-1)+h/2,'spline');


r1 = -kappa*conj(q1);
r1s = -kappa*conj(q1s);
r2 = -kappa*conj(q2);
r2s = -kappa*conj(q2s);


q1 = h*q1;
r1 = h*r1;
q2 = h*q2;
r2 = h*r2;
q1s = h*q1s;
r1s = h*r1s;
q2s = h*q2s;
r2s = h*r2s;


phi1 = ones(1,length(xi),'like',1+1i);
phi2 = zeros(1,length(xi),'like',1+1i);
phi3 = zeros(1,length(xi),'like',1+1i);



exp_vec = exp(2i*xi*t(1));
q1exp = q1(1).*exp_vec;
r1exp = r1(1)./exp_vec;
q2exp = q2(1).*exp_vec;
r2exp = r2(1)./exp_vec;

for n = 1:1:N-1
    
    sexp_vec = exp(2i*xi*ts(n));
    q1sexp = q1s(n).*sexp_vec;
    r1sexp = r1s(n)./sexp_vec;
    q2sexp = q2s(n).*sexp_vec;
    r2sexp = r2s(n)./sexp_vec;
    
    k11 = q1exp.*phi2 + q2exp.*phi3;
    k12 = r1exp.*phi1;
    k13 = r2exp.*phi1;    
    
    k21 = q1sexp.*(phi2+k12/2) + q2sexp.*(phi3+k13/2);
    k22 = r1sexp.*(phi1+k11/2);
    k23 = r2sexp.*(phi1+k11/2);
    
    k31 = q1sexp.*(phi2+k22/2) + q2sexp.*(phi3+k23/2);
    k32 = r1sexp.*(phi1+k21/2);
    k33 = r2sexp.*(phi1+k21/2);
    
    exp_vec = exp(2i*xi*t(n+1));
    q1exp = q1(n+1).*exp_vec;
    r1exp = r1(n+1)./exp_vec;
    q2exp = q2(n+1).*exp_vec;
    r2exp = r2(n+1)./exp_vec;
    
    k41 = q1exp.*(phi2+k32) + q2exp.*(phi3+k33);
    k42 = r1exp.*(phi1+k31);
    k43 = r2exp.*(phi1+k31);
    
    phi1 = phi1 + (k11+2*k21+2*k31+k41)/6;
    phi2 = phi2 + (k12+2*k22+2*k32+k42)/6;
    phi3 = phi3 + (k13+2*k23+2*k33+k43)/6;
end
a = phi1;
b1 = phi2;
b2 = phi3;
end
