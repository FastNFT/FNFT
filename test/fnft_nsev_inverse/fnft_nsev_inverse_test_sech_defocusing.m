% This Matlab script generates the files
% fnft_nsev_inverse_test_sech_defocusing_....inc

clear all;

for D=[2048 4096]
    fnft_nsev_inverse_test_sech_defocusing_aux(D);
end

function fnft_nsev_inverse_test_sech_defocusing_aux(D, filename)
T = 2.0*[-1 1.1];
M = 2*D;
XI = mex_fnft_nsev_inverse_XI(D, T, M);

kappa = -1;
Q = 1;
GAM = 1/25;
F = 1.5;
eps_xi = (XI(2) - XI(1))/(M - 1);
xi = XI(1) + (0:(M-1))*eps_xi;
cgamma = @(z) double(gamma(sym(z)));
d = 0.5 + 1i*(xi*GAM-F);
fp = 0.5 - 1i*(xi*GAM+sqrt(F^2+Q^2));
fm = 0.5 - 1i*(xi*GAM-sqrt(F^2+Q^2));
gp = 1 - 1i*(F+sqrt(F^2+Q^2));
gm = 1 - 1i*(F-sqrt(F^2+Q^2));
contspec_exact = -2^(-2i*F)*Q*cgamma(d).*cgamma(fm).* ...
    cgamma(fp)./(cgamma(conj(d)).*cgamma(gm).*cgamma(gp));
      
eps_t = (T(2) - T(1))/(D - 1);
t = T(1) + (0:D-1)*eps_t;
q_exact = -conj(Q/GAM*sech(t/GAM).^(1-2j*F));

filename = sprintf('fnft_nsev_inverse_test_sech_defocusing_data_%d.inc',D);
fileID = fopen(filename, 'w');
fprintf(fileID, 'const UINT M_%d = %d;\n', D, M);
fprintf(fileID, 'const REAL T_%d[2] = {%.17g, %.17g};\n', D, T(1), T(2));
fprintf(fileID, 'const REAL XI_%d[2] = {%.17g, %.17g};\n', D, XI(1), XI(2));
fprintf(fileID, 'const COMPLEX q_exact_%d[%d] = {\n', D, D);
for n = 1:D
    fprintf(fileID, '    %.17g + %.17g*I', real(q_exact(n)), imag(q_exact(n)));
    if n<D
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, ' };\n');
    end
end
fprintf(fileID, 'COMPLEX contspec_%d[%d] = {\n', D, M);
for n = 1:M
    fprintf(fileID, '    %.17g + %.17g*I', real(contspec_exact(n)), imag(contspec_exact(n)));
    if n<M
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, ' };\n');
    end
end
fclose(fileID);
end