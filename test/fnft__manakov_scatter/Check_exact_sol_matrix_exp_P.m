syms eps_t lam q1 q2 kappa

P = [-1i*lam,           q1,     q2;
     -kappa*conj(q1),   1i*lam, 0;
     -kappa*conj(q2),   0,      1i*lam];

TM = expm(P*eps_t);

%% simplifying the expressions

syms x1     % = (lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2)
syms x2     % = eps_t*x1*i

TM_11 = subs(TM(1,1),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_11 = subs(TM_11, eps_t*x1*1i, x2);

TM_12 = subs(TM(1,2),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_12 = subs(TM_12, eps_t*x1*1i, x2);

TM_13 = subs(TM(1,3),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_13 = subs(TM_13, eps_t*x1*1i, x2);


TM_21 = subs(TM(2,1),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_21 = subs(TM_21, eps_t*x1*1i, x2);

TM_22 = subs(TM(2,2),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_22 = subs(TM_22, eps_t*x1*1i, x2);

TM_23 = subs(TM(2,3),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_23 = subs(TM_23, eps_t*x1*1i, x2);


TM_31 = subs(TM(3,1),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_31 = subs(TM_31, eps_t*x1*1i, x2);

TM_32 = subs(TM(3,2),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_32 = subs(TM_32, eps_t*x1*1i, x2);

TM_33 = subs(TM(3,3),(lam^2 + kappa*q1*conj(q1) + kappa*q2*conj(q2))^(1/2),x1);
TM_33 = subs(TM_33, eps_t*x1*1i, x2);

%% Rewriting using goniometric functions

syms x3     % = q1*conj(q1)+q2*conj(q2) 

TM_simplified_11 = -(lam/x1)*sinh(x2)+cosh(x2);
TM_simplified_12 = -(1i*q1/x1)*sinh(x2);
TM_simplified_13 = -(1i*q2/x1)*sinh(x2);

TM_simplified_21 = (1i*kappa*conj(q1)/x1)*sinh(x2);
TM_simplified_22 = (q1*conj(q1)/x3)*(cosh(x2) + (lam/x1)*sinh(x2)) + q2*conj(q2)*exp(eps_t*lam*1i)/x3;
TM_simplified_23 = (q2*conj(q1)/x3)*(cosh(x2) - exp(eps_t*lam*1i) + (lam/x1)*sinh(x2));

TM_simplified_31 = (1i*kappa*conj(q2)/x1)*sinh(x2);
TM_simplified_32 = (q1*conj(q2)/x3)*(cosh(x2) - exp(eps_t*lam*1i) + (lam/x1)*sinh(x2));
TM_simplified_33 = (q2*conj(q2)/x3)*(cosh(x2) + (lam/x1)*sinh(x2)) + q1*conj(q1)*exp(eps_t*lam*1i)/x3;

%% Check if the simplifications are correct
% Setting values to test
eps_t_val = 0.13;
lam_val = 1+0.5*i;
q1_val = -0.058200013523445 + 0.337731590275575i;
q2_val = -0.030555007099809 + 1.049552283193580i;
kappa_val = 1;

P_val = [-1i*lam_val,           q1_val,     q2_val;
     -kappa_val*conj(q1_val),   1i*lam_val, 0;
     -kappa_val*conj(q2_val),   0,      1i*lam_val];
TM_num = expm(eps_t_val*P_val);

TM_func(eps_t, lam, q1, q2, kappa) = TM;
TM_simplified_func(eps_t, lam, q1, q2, kappa, x1, x2, x3) = ...
    [TM_simplified_11, TM_simplified_12, TM_simplified_13;
     TM_simplified_21, TM_simplified_22, TM_simplified_23;
     TM_simplified_31, TM_simplified_32, TM_simplified_33];

x1_val = (lam_val^2 + kappa_val*q1_val*conj(q1_val) + kappa_val*q2_val*conj(q2_val))^(1/2);
x2_val =  eps_t_val*x1_val*1i;
x3_val = q1_val*conj(q1_val)+q2_val*conj(q2_val); 

digits(6);
a1 = vpa(TM_func(eps_t_val, lam_val, q1_val, q2_val, kappa_val))
a2 = vpa(TM_simplified_func(eps_t_val, lam_val, q1_val, q2_val, kappa_val, x1_val, x2_val, x3_val))
