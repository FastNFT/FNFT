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
* Contributors:
* Sander Wahls (TU Delft) 2018.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h> // for printf
#include "fnft__errwarn.h"
#include "fnft__misc.h" // for misc_sech
#include "fnft__kdvv_testcases.h"

/**
 * Routine to run tests for \link fnft_kdvv \endlink
 */

INT fnft__kdvv_testcases(kdvv_testcases_t tc, const UINT D,
    COMPLEX ** const q_ptr, REAL * const T, 
    UINT * const M_ptr, COMPLEX ** const contspec_ptr,
    COMPLEX ** const ab_ptr, REAL * const XI, UINT * const K_ptr, 
    COMPLEX ** const bound_states_ptr, 
    COMPLEX ** const normconsts_ptr,
    COMPLEX ** const residues_ptr)
{
    UINT i;
    INT ret_code;
    REAL eps_t, t_i;
    
    // Avoid not used warnings. Not yet used, but will be in the future.
    (void) ab_ptr;
    (void) K_ptr;
    (void) bound_states_ptr;
    (void) normconsts_ptr;
    (void) residues_ptr;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q_ptr == NULL)
        return E_INVALID_ARGUMENT(q_ptr);
    if (T == NULL)
        return E_INVALID_ARGUMENT(T);
    if (M_ptr == NULL)
        return E_INVALID_ARGUMENT(M_ptr);
    if (contspec_ptr == NULL)
        return E_INVALID_ARGUMENT(contspec_ptr);

    // Set the number of points in the continuous spectrum *M_ptr and the
    // number of bound states *K_ptr (needed for proper allocation)
    switch (tc) {

    case kdvv_testcases_SECH:
    case kdvv_testcases_RECT:
    case kdvv_testcases_NEGATIVE_RECT:
        *M_ptr = 16;
        break;

    default:
        return E_INVALID_ARGUMENT(tc);
    } 

    // Allocate memory for results
    *q_ptr = malloc(D * sizeof(COMPLEX));
    if (*q_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_1;
    }
    *contspec_ptr = malloc((*M_ptr) * sizeof(COMPLEX));
    if (*contspec_ptr == NULL) {
        ret_code = E_NOMEM;
        goto release_mem_2;
    }

    // generate test case
    switch (tc) {

    case kdvv_testcases_SECH: 

        // This test case can be found in Sec. 5.3.1. of Trogdon et al.,
        // Physica D 241 (2012) 1003-1024. They refer back to the book of
        // Drazin and Johnson, "Solitons: An Introduction", 1996.

        // The following MATLAB code has been used to compute the values below.
        /*
        cgam = @(z) gamma(sym(z));
        at = @(k) 0.5 - 1j*k + (A + 0.25)^0.5;
        bt = @(k) 0.5 - 1j*k - (A + 0.25)^0.5;
        ct = @(k) 1 - 1j*k;
        a = @(k) cgam(at(k)).*cgam(bt(k)) ./ (cgam(ct(k)).*cgam(at(k) + bt(k) - ct(k)));
        rho = @(k) a(k).*cgam(ct(k)).*cgam(ct(k)-at(k)-bt(k)) ./ ( cgam(ct(k)-at(k)).*cgam(ct(k)-bt(k)) );

        xi = (-sym(1)/10 + (-sym(7):sym(8)))/2;
        XI = [xi(1) xi(end)]
        digits(40);
        contspec_exact = vpa(rho(xi)).'
        */

        // Space domain
        T[0] = -16.0;//-10.5;
        T[1] = 15.0;//10.0;
       
        // Space domain signal.
        for (i=0; i<D; i++)
            (*q_ptr)[i] = 3.2 * CPOW(misc_sech(T[0] + i*(T[1] - T[0])/(D - 1)), 2.0);

        // Nonlinear spectral domain
        XI[0] = -71.0 / 20.0;
        XI[1] = 79.0 / 20.0;

        // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
        (*contspec_ptr)[0] =   0.00001969133727937273786256309569057640147214 + 0.00001674031899324291878282874809046115677856*I;
        (*contspec_ptr)[1] =   0.0001042924180868180241029534277241901174759 + 0.00006768141773041039691385238506571730630862*I;
        (*contspec_ptr)[2] =    0.0005497992302347843662438396963070116445916 + 0.0002354192967860323004702520702102621027959*I;
        (*contspec_ptr)[3] =     0.002836080801758817045949972145829401374196 + 0.0004838176592384099640435004607088935813189*I;
        (*contspec_ptr)[4] =       0.01363997503332085032955915136653550260998 - 0.002341609324488054092139191512148319815275*I;
        (*contspec_ptr)[5] =        0.05166842907012522766918606650773513057274 - 0.04189670084801774699210641959191233802322*I;
        (*contspec_ptr)[6] =         0.02522525851207013024776243251692135365253 - 0.3129668484537972516927631727064614422506*I;
        (*contspec_ptr)[7] =        - 0.9662812724225684323774956123872384448904 - 0.1912816707764648616985575230120299971374*I;
        (*contspec_ptr)[8] =       - 0.05828057373883236699003469724963139197661 + 0.4183903740127110822937728254406981597783*I;
        (*contspec_ptr)[9] =        0.06250255474895942181325668026198027247793 + 0.06614727260898870840753897306353062501362*I;
        (*contspec_ptr)[10] =       0.01833190345613142519098859234810492849712 + 0.004789218463862199786383618780644774216483*I;
        (*contspec_ptr)[11] =     0.003914704641391021868806404239021266734245 - 0.0004367271876152405705478883462156928985432*I;
        (*contspec_ptr)[12] =    0.0007653177688926532754047857544255194940605 - 0.0002911773885213392112790377303427788297587*I;
        (*contspec_ptr)[13] =    0.0001455283498004098201178193218482810829457 - 0.0000882963579879275159326990821726513549138*I;
        (*contspec_ptr)[14] =  0.00002748286948200803407581026489337631688606 - 0.00002228910643157474631983979506634092581568*I;
        (*contspec_ptr)[15] =0.000005195158073745829592578906521934403128604 - 0.000005207592023284621408556242469680129799988*I;

        break;
    case kdvv_testcases_RECT:
        
        // This test case can be found in A.R. Osborne, “Non-linear Fourier
        // analysis for the infinite-interval Korteweg-de Vries equation I: An
        // algorithm for the direct scattering transform,” Journal of
        // Computational Physics, Aug. 1991. DOI: 10.1016/0021-9991(91)90223-8,
        // Section 8.

        // The following MATLAB code has been used to compute the values below.
        /*
         ampl = sym(1);
         zeta = sym(0:15).'*sym(pi)/sym(32);
         ell = sym(1)/sym(2);
         
         kappa = sqrt(ampl + zeta.^2);
         gamma = (kappa./zeta - zeta./kappa)/2;
         delta = (kappa./zeta + zeta./kappa)/2;
         T = exp(-2i*zeta*ell) ./ (cos(2*kappa*ell) - 1i*delta.*sin(2*kappa*ell));
         R = 1i*gamma.*sin(2*kappa*ell).*T;
         R(zeta==0) = -1; % Limit value
         
         digits(40)
         T = vpa(T);
         R = vpa(R);
         */
        
        // Space domain
        T[0] = -1.0; // First midpoint, not Lower boundary
        T[1] = 2.0; // Last midpoint, not upper boundary
        
        // Space domain signal.
        eps_t = (T[1] - T[0])/(D-1);
        for (i=0; i<D; i++)
        {
            t_i = T[0] + i*eps_t;
            if (FABS(t_i) == 0.5)
                (*q_ptr)[i] = 0.5;
            else if (FABS(t_i) < 0.5)
                (*q_ptr)[i] = 1.0;
            else
                (*q_ptr)[i] = 0.0;
        }
        
        // Nonlinear spectral domain
        XI[0] = 0;
        XI[1] = 15.0/32.0 * PI;
        
        // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
        (*contspec_ptr)[0]  = - 1.0;
        (*contspec_ptr)[1]  = - 0.950168909294581608551329904244095142701   + 0.2130102017673891500378417840233356288606*I;
        (*contspec_ptr)[2]  = - 0.8259670681548542897571659846169854269094  + 0.3699979087835312790500358270735915565188*I;
        (*contspec_ptr)[3]  = - 0.6767806409539509745025825595042910752823  + 0.4540740713560979106195823030549430788822*I;
        (*contspec_ptr)[4]  = - 0.5382408857133739127874741748685782433198  + 0.4805052452886871588208774257220069438807*I;
        (*contspec_ptr)[5]  = - 0.423964915728128576089461146964133121214   + 0.4718788053853276468294076663651060796482*I;
        (*contspec_ptr)[6]  = - 0.3345086783124015757382713160639633991983  + 0.4453881032356930181388533237109483841467*I;
        (*contspec_ptr)[7]  = - 0.265767781526924377944207428377742244071   + 0.4113655421000193656444506056022786401473*I;
        (*contspec_ptr)[8]  = - 0.2130237512448752827050659061094260154733  + 0.3753282891300210562082932840855622251981*I;
        (*contspec_ptr)[9]  = - 0.1722786988174714080503814995166025434376  + 0.3399954697709277684897246150537867616637*I;
        (*contspec_ptr)[10] = - 0.140472257463270438586495580950958949769   + 0.3065903860853114584210459327329945533711*I;
        (*contspec_ptr)[11] = - 0.1153504456953303488731532430134103922649  + 0.2755726366235698056457955851840100461578*I;
        (*contspec_ptr)[12] = - 0.09527321528843998674080853529682898793209 + 0.2470267901895161686283968721614713648304*I;
        (*contspec_ptr)[13] = - 0.07904675065447655160606575367154906844081 + 0.2208632106193342065286853347418631953582*I;
        (*contspec_ptr)[14] = - 0.06579630203367370114917406544589095133598 + 0.1969199678363887603437522689507922888574*I;
        (*contspec_ptr)[15] = - 0.05487470411329939736278578609268675453902 + 0.1750134338640520324045894423655821274481*I;
       
        break;
        case kdvv_testcases_NEGATIVE_RECT:
            
            // This test case can be found in A.R. Osborne, “Non-linear Fourier
            // analysis for the infinite-interval Korteweg-de Vries equation I: An
            // algorithm for the direct scattering transform,” Journal of
            // Computational Physics, Aug. 1991. DOI: 10.1016/0021-9991(91)90223-8,
            // Section 8.
            
            // The following MATLAB code has been used to compute the values below.
            /*
             ampl = sym(-1);
             zeta = sym(0:15).'*sym(pi)/sym(32);
             ell = sym(1)/sym(2);
             
             kappa = sqrt(ampl + zeta.^2);
             gamma = (kappa./zeta - zeta./kappa)/2;
             delta = (kappa./zeta + zeta./kappa)/2;
             T = exp(-2i*zeta*ell) ./ (cos(2*kappa*ell) - 1i*delta.*sin(2*kappa*ell));
             R = 1i*gamma.*sin(2*kappa*ell).*T;
             R(zeta==0) = -1; % Limit value
             
             digits(40)
             T = vpa(T);
             R = vpa(R);
             */
            
            // Space domain
            T[0] = -1.0; // First midpoint, not Lower boundary
            T[1] = 2.0; // Last midpoint, not upper boundary
            
            // Space domain signal.
            eps_t = (T[1] - T[0])/(D-1);
            for (i=0; i<D; i++)
            {
                t_i = T[0] + i*eps_t;
                if (FABS(t_i) == 0.5)
                    (*q_ptr)[i] = -0.5;
                else if (FABS(t_i) < 0.5)
                    (*q_ptr)[i] = -1.0;
                else
                    (*q_ptr)[i] = 0.0;
            }
            
            // Nonlinear spectral domain
            XI[0] = 0;
            XI[1] = 15.0/32.0 * PI;
            
            // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
            (*contspec_ptr)[0] = -1.0;
            (*contspec_ptr)[1] = - 0.9739467451573985957325199506505369746269 - 0.1555349345978652050355041906667439098412*I;
            (*contspec_ptr)[2] = - 0.9028461287936049160123145171199359789702 - 0.2886658747278307501061098125374803592192*I;
            (*contspec_ptr)[3] = - 0.8035893706069043424417437224776288780565 - 0.3860712396445988976800969426815119972634*I;
            (*contspec_ptr)[4] = - 0.6942959495032771276336443457983276729256 - 0.4458321143826399624140237249662353773942*I;
            (*contspec_ptr)[5] = - 0.5883026354166293688520907843779958639994 - 0.4736712903480039438303182672382621165122*I;
            (*contspec_ptr)[6] = - 0.4928616941529652270973842683008589805061 - 0.4779635646873800727933323613989818639574*I;
            (*contspec_ptr)[7] = - 0.4105790992193202275287820600991870226904 - 0.4665343758579364554211810570573020018263*I;
            (*contspec_ptr)[8] = - 0.3413610425339451031793827543119710816707 - 0.4454554051133676402777299476062467490621*I;
            (*contspec_ptr)[9] = - 0.2838730864160965326338171390016363811241 - 0.4189872470579770003111002353275753860492*I;
            (*contspec_ptr)[10] = - 0.2363914002637033829184948430988984518262 - 0.3899428589360900543040072991771807959595*I;
            (*contspec_ptr)[11] = - 0.1972200838907117762619350871506416854197 - 0.3601034786315031830009142805874726386018*I;
            (*contspec_ptr)[12] = - 0.164859568723129283134171400591204814203 - 0.3305584195242223038998178971998403180682*I;
            (*contspec_ptr)[13] = - 0.138050443728764079251355187952046416493 - 0.3019489387707935514557779796321053408133*I;
            (*contspec_ptr)[14] = - 0.115762402295287875481515944676992871465 - 0.2746325066021757723919274197888779148303*I;
            (*contspec_ptr)[15] = - 0.09716336474331456198864288349325815895773 - 0.2487897763433119534338284333930560885291*I;
            
            break;
    default: // unknown test case

        ret_code = E_INVALID_ARGUMENT(tc);
        goto release_mem_1;
    }

    return SUCCESS;

    // the code below is only executed if an error occurs

release_mem_2:
    free(*contspec_ptr);
release_mem_1:
    free(*q_ptr);
    
    return ret_code;
}

INT kdvv_testcases_test_fnft(kdvv_testcases_t tc, UINT D,
    const REAL eb[6], fnft_kdvv_opts_t * const opts) {
    COMPLEX * q = NULL;
    COMPLEX * contspec = NULL;
    REAL T[2], XI[2];
    COMPLEX * contspec_exact = NULL;
    COMPLEX * ab_exact = NULL;
    COMPLEX * bound_states_exact = NULL;
    COMPLEX * normconsts_exact = NULL;
    COMPLEX * residues_exact = NULL;
    UINT M, K;
    REAL errs[6];
    INT ret_code;

    // Load test case
    ret_code = kdvv_testcases(tc, D, &q, T, &M, &contspec_exact, &ab_exact,
        XI, &K, &bound_states_exact, &normconsts_exact, &residues_exact);
    CHECK_RETCODE(ret_code, release_mem);
 
    // Allocate memory
    contspec = malloc(M * sizeof(COMPLEX));
    if ( q == NULL || contspec == NULL ) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    // Compute the NFT
    ret_code = fnft_kdvv(D, q, T, M, contspec, XI, NULL, NULL, NULL, opts);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the error(s). So far, only the error in the continuous spectrum
    // is computed.
    errs[0] = misc_rel_err(M, contspec, contspec_exact);
    for (UINT i=1; i<6; i++)
        errs[i] = FNFT_INF;
    printf("kdvv_testcases_test_fnft: %2.1e <= %2.1e\n", errs[0], eb[0]);
    //print_buf2(M, contspec, "contspec_test");
    //print_buf2(M, contspec_exact, "contspec_exact");

    // Check if the errors are below the specified bounds. Organized such that
    // the line number tells us which error was too high. The conditions are
    // written in this way to ensure that they fail if an error is NAN.
    if (!(errs[0] <= eb[0])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[1] <= eb[1])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[2] <= eb[2])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[3] <= eb[3])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[4] <= eb[4])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[5] <= eb[5])) {
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }

    ///// Clean up /////

release_mem:
    free(q);
    free(contspec);
    free(contspec_exact);

    return ret_code;
}
