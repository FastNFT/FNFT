/**
 * This file is part of FNFT.  
 *                                                                  
 * FNFT is free software; you can redistribute it and/or
 * modify it under the terms of the version 2 of the GNU General
 * Public License as published by the Free Software Foundation.
 *
 * FNFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                      
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 * Sander Wahls (KIT) 2023.
 **/

/* The genus one solution of the KdV equation is passed to fnft_kdvp and the
resulting main and auxiliary spectra are compared with the exact values.

The data files were generated using the Matlab functions

function [q, mu, L, m] = genus_one(E1, E2, E3, phi0, x, t)
% Computes the genus one solution of the KdV equation u_x+6uu_t+u_xxx=0.
% See p. 144 of the book "Solitons and the Inverse Scattering Transform"
% by Ablowitz and Segur (1981). There is a sign difference because the
% KdV equation in that book is u_x-6uu_t+u_xxx=0.
m = (E3 - E2) / (E3 - E1);
K = ellipke(m);
L = 2/sqrt(E3 - E1)*K;
args = sqrt(E3 - E1)*(x - 2*(E1+E2+E3)*t) + phi0;
[~,cn,~] = ellipj(args, m);
mu = (E3 - E2)*cn.^2 + E2;
q = -(E1+E2+E3) + 2*mu;
end

function test_data(D)
E1 = -1.4001;
E2 = -1.4;
E3 = 0.3;
phi0 = 2.1;
t = 0;

[~, ~, L, ~] = genus_one(E1, E2, E3, phi0, 0, t);
x = linspace(0, L, D+1);
x = x(1:end-1);

ep_x = x(2) - x(1);
[~, mu, ~, ~] = genus_one(E1, E2, E3, phi0, -ep_x/2, t);
aux_spec_exact = mu;
main_spec_exact = [E1, E2, E3];

fprintf("UINT D = %i;\n", D);
fprintf("REAL X[2] = {0, %.60e};\n", L);
fprintf("REAL E[2] = {%.60e, %.60e};\n", E1-(E3-E1), E3+(E3-E1));
fprintf("REAL main_spec_exact[3] = {%.60e, %.60e, %.60e};\n", E1, E2, E3);
fprintf("REAL aux_spec_exact[1] = {%.60e};\n", aux_spec_exact);
fprintf("COMPLEX q[%i] = {\n", D);
for n=1:D
    [qx, ~, ~, ~] = genus_one(E1, E2, E3, phi0, x(n), t);
    fprintf("\t%.60e", qx);
    if n<D
        fprintf(",\n");
    else
        fprintf(" };\n");
end
end

Specifically, test_data(D) was called for D=256, 512, 1024. */

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft_kdvp.h"
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif

static INT test_256(REAL err_bnds[2])
{
#include "fnft_kdvp_test_genus1_data256.inc"
#include "fnft_kdvp_test_genus1_test_fun.inc"
}

static INT test_512(REAL err_bnds[2])
{
#include "fnft_kdvp_test_genus1_data512.inc"
#include "fnft_kdvp_test_genus1_test_fun.inc"
}

static INT test_1024(REAL err_bnds[2])
{
#include "fnft_kdvp_test_genus1_data1024.inc"
#include "fnft_kdvp_test_genus1_test_fun.inc"
}

int main(void)
{
    REAL err_bnds[2] = {0.00013, 0.00016};
    INT ret_code = test_256(err_bnds);
    CHECK_RETCODE(ret_code, failed);

    err_bnds[0] /= 4;
    err_bnds[1] /= 4;
    ret_code = test_512(err_bnds);
    CHECK_RETCODE(ret_code, failed);

    err_bnds[0] /= 4;
    err_bnds[1] /= 4;
    ret_code = test_1024(err_bnds);
    CHECK_RETCODE(ret_code, failed);

    return EXIT_SUCCESS;
failed:
    return EXIT_FAILURE;
}
