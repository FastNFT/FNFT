eps_t = 0.13;
    kappa = +1; D=8;
    q1 = 0.4*cos(1:D)+0.5j*sin(0.3*(1:D));
    q2 = 0.21*cos(1:D)+1.05j*sin(0.2*(1:D));
    result_exact = []; lam = [2, 1+0.5*i];
    for i = 1:length(lam)
        lambda = lam(i);
        S=eye(3);
        for n=1:D
            P = [-1i*lambda,           q1(n),     q2(n);
                -kappa*conj(q1(n)),   1i*lambda, 0;
                -kappa*conj(q2(n)),   0,      1i*lambda];

            TM = expm(eps_t*P);
            S=TM*S;
        end
        result_exact = [result_exact,...
            [S(1,1) S(1,2) S(1,3) S(2,1) S(2,2) S(2,3) S(3,1) S(3,2) S(3,3)]];
    end
    format long g; result_exact.'