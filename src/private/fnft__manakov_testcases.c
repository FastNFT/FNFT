/*
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
* Contributor:
* Lianne de Vries (TU Delft) 2021.
* Original file written by:
* Sander Wahls (TU Delft) 2017-2018.
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__misc.h" // for misc_sech
#include "fnft__manakov_testcases.h"
#include "fnft__manakov_discretization.h" // for nse_discretization_degree
#include "fnft_manakov.h"
#ifdef DEBUG
#include <stdio.h>
#include <string.h> // for memcpy
#endif

INT manakov_testcases(manakov_testcases_t tc, const UINT D,
    COMPLEX ** const q1_ptr, COMPLEX ** const q2_ptr,
    REAL * const T,
    UINT * const M_ptr, COMPLEX ** const contspec_ptr,
    COMPLEX ** const ab_ptr,
    REAL * const XI, UINT * const K_ptr,
    COMPLEX ** const bound_states_ptr,
    COMPLEX ** const normconsts_ptr,
    COMPLEX ** residues_ptr, INT * const kappa_ptr)
{
    UINT i;
    INT ret_code;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q1_ptr == NULL)
        return E_INVALID_ARGUMENT(q1_ptr);
    if (q2_ptr == NULL)
        return E_INVALID_ARGUMENT(q2_ptr);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);
    if (K_ptr == NULL)
        return E_INVALID_ARGUMENT(K_ptr);
    if (contspec_ptr == NULL)
        return E_INVALID_ARGUMENT(contspec_ptr);
    if (ab_ptr == NULL)
        return E_INVALID_ARGUMENT(ab_ptr);
    if (bound_states_ptr == NULL)
        return E_INVALID_ARGUMENT(bound_state_ptr);
    if (normconsts_ptr == NULL)
        return E_INVALID_ARGUMENT(normconst_ptr);
    if (residues_ptr == NULL)
        return E_INVALID_ARGUMENT(residues_ptr);
    if (kappa_ptr == NULL)
        return E_INVALID_ARGUMENT(kappa_ptr);
    
    // Set the number of points in the continuous spectrum *M_ptr and the
    // number of bound states *K_ptr (needed for proper allocation)
    switch (tc) {
        
        case manakov_testcases_SECH_FOCUSING:
            *M_ptr = 16;
            *K_ptr = 3;
            break;
            
        case manakov_testcases_SECH_FOCUSING2:
            *M_ptr = 16;
            *K_ptr = 5;
            break;
            
        case manakov_testcases_SECH_FOCUSING_CONTSPEC:
            *M_ptr = 16;
            *K_ptr = 0;
            break;
            
        case manakov_testcases_SECH_DEFOCUSING:
            *M_ptr = 16;
            *K_ptr = 0;
            break;
            
        case manakov_testcases_TRUNCATED_SOLITON:
            *M_ptr = 16;
            *K_ptr = 0;
            break;
            
        default:
            return E_INVALID_ARGUMENT(tc);
    }
    
    // Allocate memory for results
    *q1_ptr = malloc(D * sizeof(COMPLEX));
    *q2_ptr = malloc(D * sizeof(COMPLEX));
    if (*q1_ptr == NULL || *q2_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_1;
    }
    if (*M_ptr > 0) {
        *contspec_ptr = malloc((*M_ptr) * sizeof(COMPLEX));
        if (*contspec_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_2;
        }
        *ab_ptr = malloc(2*(*M_ptr) * sizeof(COMPLEX));
        if (*ab_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_3;
        }
    } else {
        *contspec_ptr = NULL;
        *ab_ptr = NULL;
    }
    if (*K_ptr > 0) {
        *bound_states_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
        if (*bound_states_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_4;
        }
        *normconsts_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
        if (*normconsts_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_5;
        }
        *residues_ptr = malloc((*K_ptr) * sizeof(COMPLEX));
        if (*residues_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_6;
        }
    } else {
        *bound_states_ptr = NULL;
        *normconsts_ptr = NULL;
        *residues_ptr = NULL;
    }

    // generate test case
    switch (tc) {

    case manakov_testcases_SECH_FOCUSING:

        // This test case can be found in J. Satsuma & N. Yajima, Prog. Theor.
        // Phys. Suppl., No. 55, p. 284ff, 1974.

        // The following MATLAB code has been used to compute the values below.
        /*
        A = sym(3.2);
        hlf = sym(1)/sym(2);
        im = sym(1j);
        syms lam;
        a(lam) = gamma(-im*lam + hlf).^2 ./ ( gamma(-im*lam + A + hlf) ...
            * gamma(-im*lam - A + hlf) );
        da = diff(a, lam);
        b(lam) = im*sin(sym(pi)*A) ./ cosh(sym(pi) * lam);

        bound_states = sym(1j)*(A - (floor(A):-1:1) + hlf)
        normconsts = b(bound_states)
        syms h;
        da_vals = limit(da(bound_states + h), h, 0);
        residues = normconsts ./ da_vals

        xi = (-sym(7):sym(8))/5;
        XI = [xi(1) xi(end)]
        digits(40);
        contspec = vpa(b(xi) ./ a(xi)).'
        ab = vpa([a(xi) b(xi)]).'
        */
        // TODO: write matlab and copy code

        // Time domain
        T[0] = -25.0;
        T[1] = 25.0;

        // Time domain signal. Note that the signal here has an extra factor I
        // when compared to the reference because there a different variant of
        // the Zakharov-Shabat problem is used. By multiplying the signal with
        // I, we get exactly the a() and b() provided there for the signal
        // A*sech(t) without an extra I.
        // Test signal: q1 = 0.8*sech(t), q2 = 5.2*sech(t)
        for (i=0; i<D; i++){
            (*q1_ptr)[i] = I * 0.8 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
            (*q2_ptr)[i] = I * 5.2 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
        }
        // Nonlinear spectral domain
        XI[0] = -7.0 / 5.0;
        XI[1] = 8.0 / 5.0;

        // TODO: fill in right values
        // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
        (*contspec_ptr)[0] = 0.01418788668709013595992021025473684215882
            - 0.002780494135403430187512334448840356894033*I;
        (*contspec_ptr)[1] = 0.02241182921176242546967787829962193792932
            - 0.01523064097337719809817988712025498134135*I;
        (*contspec_ptr)[2] = 0.02465858976907935401910832019686315816969
            - 0.04438144029159159323809699075974745310038*I;
        (*contspec_ptr)[3] = -0.003769955008351912170484991343204667254443
            - 0.09495492126539469847498823350667910603724*I;
        (*contspec_ptr)[4] = -0.1125183095807641031140695325799053429154
            - 0.1368780622436166383108546963189442789103*I;
        (*contspec_ptr)[5] = -0.3233468273512117378446470596682895379774
            - 0.03729229993634457693343890962907380139133*I;
        (*contspec_ptr)[6] = -0.4116152074717953525393539594551519006178
            + 0.3788164254609808107108977070302734836415*I;
        (*contspec_ptr)[7] = 0.7265425280053608858954667574806187496161*I;
        (*contspec_ptr)[8] = 0.4116152074717953525393539594551519006178
            + 0.3788164254609808107108977070302734836415*I;
        (*contspec_ptr)[9] = 0.3233468273512117378446470596682895379774
            - 0.03729229993634457693343890962907380139133*I;
        (*contspec_ptr)[10] = 0.1125183095807641031140695325799053429154
            - 0.1368780622436166383108546963189442789103*I;
        (*contspec_ptr)[11] = 0.003769955008351912170484991343204667254443
            - 0.09495492126539469847498823350667910603724*I;
        (*contspec_ptr)[12] = -0.02465858976907935401910832019686315816969
            - 0.04438144029159159323809699075974745310038*I;
        (*contspec_ptr)[13] = -0.02241182921176242546967787829962193792932
            - 0.01523064097337719809817988712025498134135*I;
        (*contspec_ptr)[14] = -0.01418788668709013595992021025473684215882
            - 0.002780494135403430187512334448840356894033*I;
        (*contspec_ptr)[15] = -0.007616495541673034905010512686300960103961
            + 0.001218250087174144753540506062403435472519*I;

        // a(xi)
        (*ab_ptr)[0] =   0.1922981551089552207414790788104399336051
            - 0.9812300627012406776576449201797836254391*I;
        (*ab_ptr)[1] =   0.561866472778759866094679757854572199885
            - 0.8267843388695356739056977390134312187561*I;
        (*ab_ptr)[2] =   0.8730141099108177240539975329792577203606
            - 0.4850517842024005091448564057792135627085*I;
        (*ab_ptr)[3] =   0.99473134297967765839265519878851030091
            + 0.03949339706100667801810586124473818981564*I;
        (*ab_ptr)[4] =   0.7606487443184108526272926017586175655859
            + 0.6252785106141419778427378472314923857149*I;
        (*ab_ptr)[5] =   0.1089468485775085556086575508888526656452
            + 0.9446351632262268082615533423570763855373*I;
        (*ab_ptr)[6] = - 0.590997041209672540861257403754476421502
            + 0.6421669003309714788364298184082912989208*I;
        (*ab_ptr)[7] = - 0.8090169943749474241022934171828190588602;
        (*ab_ptr)[8] = - 0.590997041209672540861257403754476421502
            - 0.6421669003309714788364298184082912989208*I;
        (*ab_ptr)[9] =   0.1089468485775085556086575508888526656452
            - 0.9446351632262268082615533423570763855373*I;
        (*ab_ptr)[10] =   0.7606487443184108526272926017586175655859
            - 0.6252785106141419778427378472314923857149*I;
        (*ab_ptr)[11] =   0.99473134297967765839265519878851030091
            - 0.03949339706100667801810586124473818981564*I;
        (*ab_ptr)[12] =   0.8730141099108177240539975329792577203606
            + 0.4850517842024005091448564057792135627085*I;
        (*ab_ptr)[13] =   0.561866472778759866094679757854572199885
            + 0.8267843388695356739056977390134312187561*I;
        (*ab_ptr)[14] =   0.1922981551089552207414790788104399336051
            + 0.9812300627012406776576449201797836254391*I;
        (*ab_ptr)[15] = - 0.1579366040628854761179976110287556423348
            + 0.9874191295994487368292392806331584236124*I;

        // b(xi)
        (*ab_ptr)[16] = - 0.01445626483610090114054909784096292396108*I;
        (*ab_ptr)[17] = - 0.02708733591957504791796401942967379252746*I;
        (*ab_ptr)[18] = - 0.05070631655613093722786537802129497025788*I;
        (*ab_ptr)[19] = - 0.09460352468290259430624658322771974638099*I;
        (*ab_ptr)[20] = - 0.1744714072018253986997192102912036857244*I;
        (*ab_ptr)[21] = - 0.3095076615878664435099104772640477540624*I;
        (*ab_ptr)[22] = - 0.4882050485203166754374126644029116788666*I;
        (*ab_ptr)[23] = - 0.5877852522924731291687059546390727685977*I;
        (*ab_ptr)[24] = - 0.4882050485203166754374126644029116788666*I;
        (*ab_ptr)[25] = - 0.3095076615878664435099104772640477540624*I;
        (*ab_ptr)[26] = - 0.1744714072018253986997192102912036857244*I;
        (*ab_ptr)[27] = - 0.09460352468290259430624658322771974638099*I;
        (*ab_ptr)[28] = - 0.05070631655613093722786537802129497025788*I;
        (*ab_ptr)[29] = - 0.02708733591957504791796401942967379252746*I;
        (*ab_ptr)[30] = - 0.01445626483610090114054909784096292396108*I;
        (*ab_ptr)[31] = - 0.007713079680024468575954854247144612396039*I;

        // Discrete spectrum
        (*bound_states_ptr)[0] = 7.0*I / 10.0;
        (*bound_states_ptr)[1] = 17.0*I / 10.0;
        (*bound_states_ptr)[2] = 27*I / 10.0;
        (*normconsts_ptr)[0] = 1.0*I;
        (*normconsts_ptr)[1] = -1.0*I;
        (*normconsts_ptr)[2] = 1.0*I;
        (*residues_ptr)[0] = -1428.0 * GAMMA(2.0/5.0) \
            / ( 25.0 * CPOW(GAMMA(1.0/5.0), 2.0) );
        (*residues_ptr)[1] = -5236.0 * GAMMA(2.0/5.0) \
            / ( 15.0 * CPOW(GAMMA(1.0/5.0), 2.0) );
        (*residues_ptr)[2] = -4284.0 * GAMMA(2.0/5.0) \
            / (11.0 * CPOW(GAMMA(1.0/5.0), 2.0) );

        *kappa_ptr = +1;

        break;

        // TODO: all other cases, line 294 - 705 in fnft__nsev_testcases.c that follow below

    default: // unknown test case

        ret_code = E_INVALID_ARGUMENT(tc);
        goto release_mem_5;
    }

    return SUCCESS;

    // the code below is only executed if an error occurs

release_mem_6:
    free(*residues_ptr);
release_mem_5:
    free(*normconsts_ptr);
release_mem_4:
    free(*bound_states_ptr);
release_mem_3:
    free(*ab_ptr);
release_mem_2:
    free(*contspec_ptr);
release_mem_1:{
    free(*q1_ptr);
    free(*q2_ptr);
}

    return ret_code;
}

// Compares computed with exact nonlinear Fourier spectrum.
static INT manakov_compare_nfs(const UINT M, const UINT K1, const UINT K2,
    COMPLEX const * const contspec_1,
    COMPLEX const * const contspec_2,
    COMPLEX const * const ab_1,
    COMPLEX const * const ab_2,
    COMPLEX const * const bound_states_1,
    COMPLEX const * const bound_states_2,
    COMPLEX const * const normconsts_1,
    COMPLEX const * const normconsts_2,
    COMPLEX const * const residues_1,
    COMPLEX const * const residues_2,
    REAL dists[6])
{
    UINT i, j, min_j = 0;
    REAL dist, min_dist, nrm;

    // Check last argument
    if (dists == NULL)
        return E_INVALID_ARGUMENT(dists);

    if (contspec_1 == NULL || contspec_2 == NULL || M == 0)
        dists[0] = NAN;
    else {
        dists[0] = 0.0;
        nrm = 0.0;
        for (i=0; !(i>=M); i++) {
            dists[0] += CABS(contspec_1[i] - contspec_2[i]);
            nrm += CABS(contspec_2[i]);
        }
        if (nrm > 0)
            dists[0] /= nrm;
    }
    if (ab_1 == NULL || ab_2 == NULL || M == 0) {
        dists[1] = NAN;
        dists[2] = NAN;
    }
    else {
        // error in a
        dists[1] = 0.0;
        nrm = 0.0;
        for (i=0; !(i>=M); i++) {
            dists[1] += CABS(ab_1[i] - ab_2[i]);
            nrm += CABS(ab_2[i]);
        }
        if (nrm > 0)
            dists[1] /= nrm;

        // error in b
        dists[2] = 0.0;
        nrm = 0.0;
        for (i=M; !(i>=2*M); i++) {
            dists[2] += CABS(ab_1[i] - ab_2[i]);
            nrm += CABS(ab_2[i]);
        }
        if (nrm > 0)
            dists[2] /= nrm;

    }
    if (K1 == 0 && K2 == 0) {
        dists[3] = 0.0;
        dists[4] = 0.0;
        dists[5] = 0.0;
    } else if (bound_states_1 == NULL || bound_states_2 == NULL || K1 == 0
    || K2 == 0) {
        dists[3] = NAN;
    } else {
        dists[3] = misc_hausdorff_dist(K1, bound_states_1, K2, bound_states_2);

        if (normconsts_1 == NULL || normconsts_2 == NULL)
            dists[4] = NAN;
        else {
            dists[4] = 0.0;
            nrm = 0.0;
            for (i=0; !(i>=K1); i++) {
                min_dist = INFINITY;
                for (j=0; !(j>=K2); j++) {
                    dist = CABS(bound_states_1[i] - bound_states_2[j]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_j = j;
                    }
                }

                dists[4] += CABS(normconsts_1[i] - normconsts_2[min_j]);
                nrm += CABS(normconsts_2[i]);
           }
            if (nrm > 0)
                dists[4] /= nrm;
        }

        if (residues_1 == NULL || residues_2 == NULL)
            dists[5] = NAN;
        else {
            dists[5] = 0.0;
            nrm = 0.0;
            for (i=0; !(i>=K1); i++) {
                min_dist = INFINITY;
                for (j=0; !(j>=K2); j++) {
                    dist = CABS(bound_states_1[i] - bound_states_2[j]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_j = j;
                    }
                }
                dists[5] += CABS(residues_1[i] - residues_2[min_j]);
                nrm += CABS(residues_2[i]);
            }
            if (nrm > 0)
                dists[5] /= nrm;
        }
    }

    return SUCCESS;
}

INT manakov_testcases_test_fnft(manakov_testcases_t tc, UINT D,
const REAL error_bounds[7], fnft_manakov_opts_t * const opts) {
    COMPLEX * q1 = NULL;
    COMPLEX * q2 = NULL;
    COMPLEX * contspec = NULL;
    COMPLEX * bound_states = NULL;
    COMPLEX * normconsts_and_residues = NULL;
    REAL T[2], XI[2];
    COMPLEX * contspec_exact = NULL;
    COMPLEX * ab_exact = NULL;
    COMPLEX * bound_states_exact = NULL;
    COMPLEX * normconsts_exact = NULL;
    COMPLEX * residues_exact = NULL;
    UINT K, K_exact=0, M;
    INT kappa = 0;
    REAL errs[7] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN };
    INT ret_code;

    // re-defining some variables to make sure everything is the same
       tc = manakov_testcases_SECH_FOCUSING;
       D = 16;
    // TODO: remove this

    // Check inputs
    if (opts == NULL)
        return E_INVALID_ARGUMENT(opts);

    // Load test case
    ret_code = manakov_testcases(tc, D, &q1, &q2, T, &M, &contspec_exact, &ab_exact,
        XI, &K_exact,&bound_states_exact, &normconsts_exact, &residues_exact,
        &kappa);
    CHECK_RETCODE(ret_code, release_mem);
/*
    // Allocate memory
    contspec = malloc(3*M * sizeof(COMPLEX));
    K = manakov_discretization_degree(opts->discretization) * D;
    
    // In every case other than that of slow discretization with no bound states
    // memory has to be alloted for the bound states
    if (!(K == 0 && K_exact == 0)){
        if (K == 0 && K_exact != 0)
            K = K_exact;        
#ifdef DEBUG
        printf("K = %ld, K_exact = %ld\n",K,K_exact);
#endif
        bound_states = malloc(K * sizeof(COMPLEX));
        normconsts_and_residues = malloc(2*D * sizeof(COMPLEX));
        if (bound_states == NULL || normconsts_and_residues == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
        }
    }
    if ( q1 == NULL || q2 == NULL || contspec == NULL) {
            ret_code = E_NOMEM;
            goto release_mem;
    }
    
    if (opts->bound_state_localization == manakov_bsloc_NEWTON){
        memcpy(bound_states,bound_states_exact,K * sizeof(COMPLEX));
    }
    
// Compute the NFT
    opts->contspec_type = fnft_manakov_cstype_BOTH;
    opts->discspec_type = fnft_manakov_dstype_BOTH;
    ret_code = fnft_manakov(D, q1, q2, T, M, contspec, XI, &K, bound_states,
        normconsts_and_residues, kappa, opts);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the errors
    ret_code = manakov_compare_nfs(M, K, K_exact, contspec, contspec_exact,
        contspec+M, ab_exact, bound_states, bound_states_exact,
        normconsts_and_residues, normconsts_exact, normconsts_and_residues+K,
        residues_exact, errs);
    CHECK_RETCODE(ret_code, release_mem);

*/
#ifdef DEBUG
    for (UINT i=0; i<6; i++)
        printf("manakov_testcases_test_fnft: error_bounds[%i] = %2.1e <= %2.1e\n",
            (int)i, errs[i], error_bounds[i]);
    //misc_print_buf(M, contspec, "r_num");
    //misc_print_buf(M, contspec_exact, "r_exact");
    //misc_print_buf(2*M, contspec+M, "ab_num");
     //misc_print_buf(2*M, ab_exact, "ab_exact");
#endif

    // Check if the errors are below the specified bounds. Organized such that
    // the line number tells us which error was too high. The conditions are
    // written in this way to ensure that they fail if an error is NAN.
    if (!(errs[0] <= error_bounds[0])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[1] <= error_bounds[1])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[2] <= error_bounds[2])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[3] <= error_bounds[3])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[4] <= error_bounds[4])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[5] <= error_bounds[5])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    ///// Clean up /////

release_mem:
    free(q1);
    free(q2);
    free(contspec);
    free(ab_exact);
    free(bound_states);
    free(normconsts_and_residues);
    free(bound_states_exact);
    free(contspec_exact);
    free(normconsts_exact);
    free(residues_exact);

    return ret_code;
}
