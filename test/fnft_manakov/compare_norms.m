% For D = 512
M = 16;

contspec_num = [8.857878018066e-002+1.235928273764e-001j, 1.412952988180e-001+5.618741686497e-002j, 1.458301777789e-001+-4.306909021516e-002j, 7.339667135989e-002+-1.331702516198e-001j, -6.236758747611e-002+-1.386783015314e-001j, -1.520337376168e-001+-2.670189038379e-003j, -2.473800832870e-002+1.500313908059e-001j, 1.520571842539e-001+-1.924579873106e-016j, -2.473800832869e-002+-1.500313908059e-001j, -1.520337376168e-001+2.670189038391e-003j, -6.236758747611e-002+1.386783015314e-001j, 7.339667135990e-002+1.331702516197e-001j, 1.458301777789e-001+4.306909021515e-002j, 1.412952988176e-001+-5.618741686527e-002j, 8.857878018080e-002+-1.235928273761e-001j, 2.038116848616e-002+-1.506850862719e-001j, 5.757620711743e-001+8.033533779466e-001j, 9.184194423173e-001+3.652182096223e-001j, 9.478961555627e-001+-2.799490863986e-001j, 4.770783638393e-001+-8.656066355289e-001j, -4.053893185947e-001+-9.014089599542e-001j, -9.882192945092e-001+-1.735622874946e-002j, -1.607970541365e-001+9.752040402386e-001j, 9.883716976505e-001+7.330542922333e-016j, -1.607970541365e-001+-9.752040402386e-001j, -9.882192945092e-001+1.735622874954e-002j, -4.053893185947e-001+9.014089599541e-001j, 4.770783638394e-001+8.656066355281e-001j, 9.478961555626e-001+2.799490863985e-001j, 9.184194423147e-001+-3.652182096243e-001j, 5.757620711752e-001+-8.033533779449e-001j, 1.324775951600e-001+-9.794530607676e-001j];
contspec_exact = [8.787570611436e-004+2.408088150114e-004j, 1.875609592411e-003+-6.897800621394e-004j, 2.167561273750e-003+-3.809743560787e-003j, -2.707950757808e-003+-9.225347158596e-003j, -1.981270952499e-002+-7.248964492162e-003j, -3.377817869481e-002+3.172732676694e-002j, 4.362974545088e-002+9.078327397334e-002j, 1.631306710314e-001+0.000000000000e+000j, 4.362974545088e-002+-9.078327397334e-002j, -3.377817869481e-002+-3.172732676694e-002j, -1.981270952499e-002+7.248964492162e-003j, -2.707950757808e-003+9.225347158596e-003j, 2.167561273750e-003+3.809743560787e-003j, 1.875609592411e-003+6.897800621394e-004j, 8.787570611436e-004+-2.408088150114e-004j, 2.832750817195e-004+-3.038702430588e-004j, 5.711920897433e-003+1.565257297574e-003j, 1.219146235067e-002+-4.483570403906e-003j, 1.408914827938e-002+-2.476333314512e-002j, -1.760167992575e-002+-5.996475653088e-002j, -1.287826119124e-001+-4.711826919906e-002j, -2.195581615162e-001+2.062276239851e-001j, 2.835933454307e-001+5.900912808267e-001j, 1.060349361704e+000+0.000000000000e+000j, 2.835933454307e-001+-5.900912808267e-001j, -2.195581615162e-001+-2.062276239851e-001j, -1.287826119124e-001+4.711826919906e-002j, -1.760167992575e-002+5.996475653088e-002j, 1.408914827938e-002+2.476333314512e-002j, 1.219146235067e-002+4.483570403906e-003j, 5.711920897433e-003+-1.565257297574e-003j, 1.841288031177e-003+-1.975156579882e-003j];

ab_num = [-3.997128586174e-018+1.107760412527e-018j, -3.887017770492e-018+-1.419598613831e-018j, -2.047866921310e-018+-3.585194636521e-018j, 1.154375242838e-018+-3.950656689592e-018j, 3.830381201243e-018+-1.405063012317e-018j, 2.870750428945e-018+2.693742257059e-018j, -1.484715278655e-018+3.091117434982e-018j, -2.804099716273e-018+6.361374339702e-032j, -1.484715278655e-018+-3.091117434982e-018j, 2.870750428945e-018+-2.693742257059e-018j, 3.830381201243e-018+1.405063012317e-018j, 1.154375242838e-018+3.950656689592e-018j, -2.047866921310e-018+3.585194636521e-018j, -3.887017770493e-018+1.419598613831e-018j, -3.997128586174e-018+-1.107760412528e-018j, -2.824667805602e-018+-3.052166997601e-018j, 1.056048642364e-029+-1.237363375902e-029j, 1.607989895452e-029+-2.576663250210e-030j, 1.318975782487e-029+9.569127370706e-030j, 8.958918300099e-031+1.625764843186e-029j, -1.328983112622e-029+9.211018982083e-030j, -1.303364560077e-029+-8.612860943339e-030j, 4.835234724405e-030+-1.273183656134e-029j, 1.113940051250e-029+9.941192459103e-043j, 4.835234724408e-030+1.273183656134e-029j, -1.303364560077e-029+8.612860943340e-030j, -1.328983112622e-029+-9.211018982082e-030j, 8.958918300098e-031+-1.625764843186e-029j, 1.318975782487e-029+-9.569127370705e-030j, 1.607989895452e-029+2.576663250210e-030j, 1.056048642364e-029+1.237363375903e-029j, 1.101665161205e-030+1.620853236624e-029j, 6.864316175366e-029+-8.042861943366e-029j, 1.045193432044e-028+-1.674831112636e-029j, 8.573342586166e-029+6.219932790959e-029j, 5.823296895065e-030+1.056747148071e-028j, -8.638390232043e-029+5.987162338354e-029j, -8.471869640499e-029+-5.598359613171e-029j, 3.142902570863e-029+-8.275693764873e-029j, 7.240610333126e-029+6.380730273057e-042j, 3.142902570865e-029+8.275693764872e-029j, -8.471869640501e-029+5.598359613171e-029j, -8.638390232043e-029+-5.987162338353e-029j, 5.823296895064e-030+-1.056747148071e-028j, 8.573342586166e-029+-6.219932790958e-029j, 1.045193432044e-028+1.674831112636e-029j, 6.864316175365e-029+8.042861943367e-029j, 7.160823547830e-030+1.053554603806e-028j];
ab_exact = [-9.644260197536e-001+2.642849739162e-001j, -9.384622196167e-001+-3.451318071639e-001j, -4.943100550496e-001+-8.688079880664e-001j, 2.810892903972e-001+-9.576046680321e-001j, 9.302059085550e-001+-3.403385888745e-001j, 6.972263717197e-001+6.548940701021e-001j, -3.611240215511e-001+7.514144455356e-001j, -6.818433382465e-001+0.000000000000e+000j, -3.611240215511e-001+-7.514144455356e-001j, 6.972263717197e-001+-6.548940701021e-001j, 9.302059085550e-001+3.403385888745e-001j, 2.810892903972e-001+9.576046680321e-001j, -4.943100550496e-001+8.688079880664e-001j, -9.384622196167e-001+3.451318071639e-001j, -9.644260197536e-001+-2.642849739162e-001j, -6.818818551989e-001+-7.314572246134e-001j, 9.111383262032e-004+0.000000000000e+000j, 1.998253780620e-003+0.000000000000e+000j, 4.381382970647e-003+0.000000000000e+000j, 9.595411460231e-003+0.000000000000e+000j, 2.089700181070e-002+0.000000000000e+000j, 4.432907513454e-002+0.000000000000e+000j, 8.397161261306e-002+0.000000000000e+000j, 1.112295613064e-001+0.000000000000e+000j, 8.397161261306e-002+0.000000000000e+000j, 4.432907513454e-002+0.000000000000e+000j, 2.089700181070e-002+0.000000000000e+000j, 9.595411460231e-003+0.000000000000e+000j, 4.381382970647e-003+0.000000000000e+000j, 1.998253780620e-003+0.000000000000e+000j, 9.111383262032e-004+0.000000000000e+000j, 4.154282228850e-004+0.000000000000e+000j, 5.922399120321e-003+0.000000000000e+000j, 1.298864957403e-002+0.000000000000e+000j, 2.847898930921e-002+0.000000000000e+000j, 6.237017449150e-002+0.000000000000e+000j, 1.358305117695e-001+0.000000000000e+000j, 2.881389883745e-001+0.000000000000e+000j, 5.458154819849e-001+0.000000000000e+000j, 7.229921484916e-001+0.000000000000e+000j, 5.458154819849e-001+0.000000000000e+000j, 2.881389883745e-001+0.000000000000e+000j, 1.358305117695e-001+0.000000000000e+000j, 6.237017449150e-002+0.000000000000e+000j, 2.847898930921e-002+0.000000000000e+000j, 1.298864957403e-002+0.000000000000e+000j, 5.922399120321e-003+0.000000000000e+000j, 2.700283448752e-003+0.000000000000e+000j];


% result after mult
res11_old = [2.123565488306e-001+2.008804919101e+000j, 3.869491549012e+000+2.896829230226e+000j, 3.515509461454e+000+-6.431076977047e-001j, -1.075928547120e+000+1.378144235338e+000j];
res11_new = [7.749643430445e+000+6.999573424543e+001j, 1.360308837040e+002+1.031460987738e+002j, 1.216927846546e+002+-2.130233673008e+001j, -3.620901713319e+001+4.945932545382e+001j];

a_num = ab_num(1:M);
a_ex = ab_exact(1:M);
b1_num = ab_num(M+1:2*M);
b1_ex = ab_exact(M+1:2*M);

ab_abs_values = [abs(ab_num)', abs(ab_exact)'];
contspec_abs_values = [abs(contspec_num)', abs(contspec_exact)'];
ab_err = (abs(ab_num)-abs(ab_exact))';
contspec_err = (abs(contspec_num)-abs(contspec_exact))';

rel_err_a = zeros(M,1); rel_err_b1 = zeros(M,1); rel_err_b2 = zeros(M,1);
for i = 1:M
    rel_err_a(i) = (abs(ab_num(i))-abs(ab_exact(i)))/abs(ab_exact(i));
    rel_err_b1(i) = (abs(ab_num(i+M))-abs(ab_exact(i+M)))/abs(ab_exact(i+M));
    rel_err_b2(i) = (abs(ab_num(i+2*M))-abs(ab_exact(i+2*M)))/abs(ab_exact(i+2*M));
end

rel_err_contspec = zeros(2*M,1);
for i = 1:2*M
    rel_err_contspec(i) = (abs(contspec_num(i))-abs(contspec_exact(i)))/abs(contspec_exact(i));
end

%% code to try out phase factors
eps_t = 50/511;
t_step_sizes = eps_t;
T = [-25-0.5*eps_t 25+0.5*eps_t];

D = 512;
b1_num = ab_num(M+1:2*M);
b1_exact = ab_exact(M+1:2*M);

phase_factor = +eps_t/3;
phase_factor = -eps_t*D - (T(2)+eps_t*0.5) - (T(1)-eps_t*0.5)% + eps_t/3
%phase_factor = -eps_t*(D-1) + (T1+eps_t*0.5) - (T0-eps_t*0.5)     % also works pretty well for real part, slightly worse for imag part
phase_factor = 0
xi = (-sym(7):sym(8))/4;        % xi has 16 elements, M = 16

b1_C = zeros(M,1);
for i = 1:M
    b1_C(i) = b1_num(i)*exp(1i*xi(i)*phase_factor);
end
b1_side_by_side = [b1_exact', b1_C];
real_err = abs(sum(real(b1_C))+sum(real(b1_exact)))
imag_err = abs(sum(imag(b1_C))-sum(imag(b1_exact)))