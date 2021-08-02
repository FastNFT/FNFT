%% This function determines the exact solution of the Manakov equation
% with two soliton potential:
function [a] = Manakov_two_soliton_exact(xi1, xi2, eta1, eta2, lam_vector)
zeta1 = xi1 + 1i*eta1;
zeta2 = xi2 + 1i*eta2;
zeta = [zeta1, zeta2];
a = ones(size(lam_vector));
for l=zeta
    a = a.*((lam_vector - l)./(lam_vector - conj(l)));
end
% a = ones(size(lam_vector));
% a = a.*((lam_vector - zeta1).*(lam_vector - zeta2))./((lam_vector - conj(zeta1)).*(lam_vector - conj(zeta2)));
end

