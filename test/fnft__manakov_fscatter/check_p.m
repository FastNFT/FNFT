% File to check value of p determined by manakov_fscatter

% Setting the values of the BE matrices
eps_t = 0.13;
kappa = 1;
D=4;
qlist = (0.41*cos(1:2*D)+0.59j*sin(0.28*(1:2*D)))*50;


% pre allocate p
% p11 = zeros(56,1); p12 = zeros(56,1); p13 = zeros(56,1);
% p21 = zeros(56,1); p22 = zeros(56,1); p23 = zeros(56,1);
% p31 = zeros(56,1); p32 = zeros(56,1); p33 = zeros(56,1);

p11 = zeros(28,1); p12 = zeros(28,1); p13 = zeros(28,1);
p21 = zeros(28,1); p22 = zeros(28,1); p23 = zeros(28,1);
p31 = zeros(28,1); p32 = zeros(28,1); p33 = zeros(28,1);

i=1;
for j = D:-1:1
q = qlist(2*j-1:2*j);
r = -conj(q);

B = [0      q(1)    q(2);
     r(1)   0       0;
     r(2)   0       0];

e_B = expm(B*eps_t);
e_1_3B = expm(B*eps_t/3);

BE_31_1 = e_1_3B(1,1); BE_31_2 = e_1_3B(1,2); BE_31_3 = e_1_3B(1,3);
BE_32_1 = e_1_3B(2,1); BE_32_2 = e_1_3B(2,2); BE_32_3 = e_1_3B(2,3);
BE_33_1 = e_1_3B(3,1); BE_33_2 = e_1_3B(3,2); BE_33_3 = e_1_3B(3,3);

BE1_1 = e_B(1,1); BE1_2 = e_B(1,2); BE1_3 = e_B(1,3);
BE2_1 = e_B(2,1); BE2_2 = e_B(2,2); BE2_3 = e_B(2,3);
BE3_1 = e_B(3,1); BE3_2 = e_B(3,2); BE3_3 = e_B(3,3);

% coefficients for z^0
p11(7*i) = (9*BE_31_1^3)/8 - BE1_1/8 + (9*BE_31_1*BE_31_2*BE_32_1)/8 + (9*BE_31_1*BE_31_3*BE_33_1)/8;
p12(7*i) = (9*BE_31_1^2*BE_31_2)/8 + (9*BE_32_1*BE_31_2^2)/8 + (9*BE_31_3*BE_33_1*BE_31_2)/8 - BE1_2/8;
p13(7*i) = (9*BE_31_1^2*BE_31_3)/8 + (9*BE_33_1*BE_31_3^2)/8 + (9*BE_31_2*BE_32_1*BE_31_3)/8 - BE1_3/8;

% z^2
p21(5+7*(i-1)) = (9*BE_31_1^2*BE_32_1)/8 + (9*BE_31_1*BE_32_1*BE_32_2)/8 + (9*BE_31_1*BE_32_3*BE_33_1)/8;
p22(5+7*(i-1)) = (9*BE_31_1*BE_31_2*BE_32_1)/8 + (9*BE_31_2*BE_32_1*BE_32_2)/8 + (9*BE_31_2*BE_32_3*BE_33_1)/8;
p23(5+7*(i-1)) = (9*BE_31_1*BE_31_3*BE_32_1)/8 + (9*BE_31_3*BE_32_1*BE_32_2)/8 + (9*BE_31_3*BE_32_3*BE_33_1)/8;
p31(5+7*(i-1)) = (9*BE_31_1^2*BE_33_1)/8 + (9*BE_31_1*BE_32_1*BE_33_2)/8 + (9*BE_31_1*BE_33_1*BE_33_3)/8;
p32(5+7*(i-1)) = (9*BE_31_1*BE_31_2*BE_33_1)/8 + (9*BE_31_2*BE_32_1*BE_33_2)/8 + (9*BE_31_2*BE_33_1*BE_33_3)/8;
p33(5+7*(i-1)) = (9*BE_31_1*BE_31_3*BE_33_1)/8 + (9*BE_31_3*BE_32_1*BE_33_2)/8 + (9*BE_31_3*BE_33_1*BE_33_3)/8;

% z^4
p11(3+7*(i-1)) = (9*BE_31_1*BE_31_2*BE_32_1)/8 + (9*BE_31_1*BE_31_3*BE_33_1)/8 + (9*BE_31_2*BE_32_1*BE_32_2)/8 + (9*BE_31_2*BE_32_3*BE_33_1)/8 + (9*BE_31_3*BE_32_1*BE_33_2)/8 + (9*BE_31_3*BE_33_1*BE_33_3)/8;
p12(3+7*(i-1)) = (9*BE_31_2*BE_32_2^2)/8 + (9*BE_31_1*BE_31_2*BE_32_2)/8 + (9*BE_31_1*BE_31_3*BE_33_2)/8 + (9*BE_31_2*BE_32_3*BE_33_2)/8 + (9*BE_31_3*BE_32_2*BE_33_2)/8 + (9*BE_31_3*BE_33_2*BE_33_3)/8;
p13(3+7*(i-1)) = (9*BE_31_3*BE_33_3^2)/8 + (9*BE_31_1*BE_31_2*BE_32_3)/8 + (9*BE_31_1*BE_31_3*BE_33_3)/8 + (9*BE_31_2*BE_32_2*BE_32_3)/8 + (9*BE_31_2*BE_32_3*BE_33_3)/8 + (9*BE_31_3*BE_32_3*BE_33_2)/8;

% z^6
p21(1+7*(i-1)) = (9*BE_31_2*BE_32_1^2)/8 - BE2_1/8 + (9*BE_32_1*BE_32_2^2)/8 + (9*BE_31_3*BE_32_1*BE_33_1)/8 + (9*BE_32_1*BE_32_3*BE_33_2)/8 + (9*BE_32_2*BE_32_3*BE_33_1)/8 + (9*BE_32_3*BE_33_1*BE_33_3)/8;
p22(1+7*(i-1)) = (9*BE_32_2^3)/8 - BE2_2/8 + (9*BE_31_2*BE_32_1*BE_32_2)/8 + (9*BE_31_3*BE_32_1*BE_33_2)/8 + (9*BE_32_2*BE_32_3*BE_33_2)/4 + (9*BE_32_3*BE_33_2*BE_33_3)/8;
p23(1+7*(i-1)) = (9*BE_32_2^2*BE_32_3)/8 + (9*BE_32_2*BE_32_3*BE_33_3)/8 + (9*BE_33_2*BE_32_3^2)/8 + (9*BE_32_3*BE_33_3^2)/8 + (9*BE_31_2*BE_32_1*BE_32_3)/8 + (9*BE_31_3*BE_32_1*BE_33_3)/8 - BE2_3/8;
p31(1+7*(i-1)) = (9*BE_31_3*BE_33_1^2)/8 - BE3_1/8 + (9*BE_33_1*BE_33_3^2)/8 + (9*BE_31_2*BE_32_1*BE_33_1)/8 + (9*BE_32_1*BE_32_2*BE_33_2)/8 + (9*BE_32_1*BE_33_2*BE_33_3)/8 + (9*BE_32_3*BE_33_1*BE_33_2)/8;
p32(1+7*(i-1)) = (9*BE_32_2^2*BE_33_2)/8 + (9*BE_32_2*BE_33_2*BE_33_3)/8 + (9*BE_31_2*BE_33_1*BE_32_2)/8 + (9*BE_32_3*BE_33_2^2)/8 + (9*BE_33_2*BE_33_3^2)/8 + (9*BE_31_3*BE_33_1*BE_33_2)/8 - BE3_2/8;
p33(1+7*(i-1)) = (9*BE_33_3^3)/8 - BE3_3/8 + (9*BE_31_2*BE_32_3*BE_33_1)/8 + (9*BE_31_3*BE_33_1*BE_33_3)/8 + (9*BE_32_2*BE_32_3*BE_33_2)/8 + (9*BE_32_3*BE_33_2*BE_33_3)/4;

i=i+1;
end

p = [p11; p12; p13; p21; p22; p23; p31; p32; p33];

% First multiply matrices 1 and 2:
    ind1 = 1:7; ind2 = 8:14;
    s11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    s12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    s13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    s21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    s22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    s23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    s31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    s32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    s33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    ind1 = 15:21; ind2 = 22:28;
    t11 = conv(p11(ind1),p11(ind2)) + conv(p12(ind1),p21(ind2)) + conv(p13(ind1),p31(ind2));
    t12 = conv(p11(ind1),p12(ind2)) + conv(p12(ind1),p22(ind2)) + conv(p13(ind1),p32(ind2));
    t13 = conv(p11(ind1),p13(ind2)) + conv(p12(ind1),p23(ind2)) + conv(p13(ind1),p33(ind2));
    t21 = conv(p21(ind1),p11(ind2)) + conv(p22(ind1),p21(ind2)) + conv(p23(ind1),p31(ind2));
    t22 = conv(p21(ind1),p12(ind2)) + conv(p22(ind1),p22(ind2)) + conv(p23(ind1),p32(ind2));
    t23 = conv(p21(ind1),p13(ind2)) + conv(p22(ind1),p23(ind2)) + conv(p23(ind1),p33(ind2));
    t31 = conv(p31(ind1),p11(ind2)) + conv(p32(ind1),p21(ind2)) + conv(p33(ind1),p31(ind2));
    t32 = conv(p31(ind1),p12(ind2)) + conv(p32(ind1),p22(ind2)) + conv(p33(ind1),p32(ind2));
    t33 = conv(p31(ind1),p13(ind2)) + conv(p32(ind1),p23(ind2)) + conv(p33(ind1),p33(ind2));

    r11 = conv(s11,t11) + conv(s12,t21) + conv(s13,t31);
    r12 = conv(s11,t12) + conv(s12,t22) + conv(s13,t32);
    r13 = conv(s11,t13) + conv(s12,t23) + conv(s13,t33);
    r21 = conv(s21,t11) + conv(s22,t21) + conv(s23,t31);
    r22 = conv(s21,t12) + conv(s22,t22) + conv(s23,t32);
    r23 = conv(s21,t13) + conv(s22,t23) + conv(s23,t33);
    r31 = conv(s31,t11) + conv(s32,t21) + conv(s33,t31);
    r32 = conv(s31,t12) + conv(s32,t22) + conv(s33,t32);
    r33 = conv(s31,t13) + conv(s32,t23) + conv(s33,t33);
    
    %% horners method
    % to check the values calculated by poly_eval
    z_cur = zlist(2);       % which of the 5 different z's we are evaluating
    r_cur = r12;            % matrix element
    tmp = r_cur(1);
    tmp_log = zeros(24,1);
    for i = 2:length(r11)
        tmp = r_cur(i)+tmp*z_cur;
        tmp_log(i-1) = tmp;
    end
