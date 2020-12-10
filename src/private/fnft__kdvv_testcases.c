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
* Sander Wahls (TU Delft) 2017-2018.
* Shrinivas Chimmalgi (TU Delft) 2019-2020.
* Peter J Prins (TU Delft) 2020.
*/
#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__errwarn.h"
#include "fnft__misc.h" // for misc_sech
#include "fnft__kdvv_testcases.h"
#include "fnft__kdv_discretization.h" // for kdv_discretization_degree
#include "fnft_kdvv.h"
#ifdef DEBUG
#include <stdio.h>
#include <string.h> // for memcpy
#endif

INT kdvv_testcases(kdvv_testcases_t tc, const UINT D,
    COMPLEX ** const q_ptr, REAL * const T,
    UINT * const M_ptr, COMPLEX ** const contspec_ptr,
    COMPLEX ** const ab_ptr,
    REAL * const XI, UINT * const K_ptr,
    COMPLEX ** const bound_states_ptr,
    COMPLEX ** const normconsts_ptr,
    COMPLEX ** residues_ptr)
{
    UINT i;
    INT ret_code;
    REAL eps_t;

    // Check inputs
    if (D < 2)
        return E_INVALID_ARGUMENT(D);
    if (q_ptr == NULL)
        return E_INVALID_ARGUMENT(q_ptr);
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

    // Set the number of points in the continuous spectrum *M_ptr and the
    // number of bound states *K_ptr (needed for proper allocation)
    switch (tc) {

        case kdvv_testcases_SECH_SQUARED:
            *M_ptr = 16;
            *K_ptr = 2;
            break;

        case kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH:
            *M_ptr = 16;
            *K_ptr = 2;
            break;

        case kdvv_testcases_SECH:
            *M_ptr = 16;
            *K_ptr = 2;
            break;

        case kdvv_testcases_RECT:
            *M_ptr = 16;
            *K_ptr = 5;
            break;

        case kdvv_testcases_NEGATIVE_RECT:
            *M_ptr = 16;
            *K_ptr = 0;
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

    case kdvv_testcases_SECH_SQUARED:

        // q(t) = 15*sech^2(2*(t-t0)), where t0 = ln(3000).
        //
        // This test case can be calculated using §2.5 in
        // G. L. Lamb, Elements of soliton theory. Wiley, 1980
        // together with the time-shift as in Appendix B in e.g.
        // P. J. Prins and S. Wahls, “Soliton phase shift calculation for the
        // Korteweg–de Vries equation,” IEEE Access, vol. 7, pp. 122914–122930,
        // September 2019, DOI: 10.1109/access.2019.2932256 .
        // The following MATLAB code has been used to compute the values below.
        /*
         A = sym(15);
         d = sym(0.5);
         exp_x0 = sym(3000);
         delta = sqrt(sym(A)*d^2+1/4);
         at = @(xi) 1/2 - 1i*xi*d + delta;
         bt = @(xi) 1/2 - 1i*xi*d - delta;
         ct = @(xi) 1 - 1i*xi*d;
         cgam = @(z) gamma(sym(z));

         rho = @(xi) exp_x0.^(-2i*xi) .* cgam(ct(xi)-at(xi)-bt(xi)) .* cgam(at(xi)) .* cgam(bt(xi)) ./ ...
                   ( cgam(ct(xi)-at(xi)) .* cgam(ct(xi)-bt(xi)) .* cgam(at(xi)+bt(xi)-ct(xi)) );
         a = @(xi) cgam(ct(xi)) .* cgam(at(xi)+bt(xi)-ct(xi)) ./ ( cgam(at(xi)) .* cgam(bt(xi)) );
         b = @(xi) exp_x0.^(-2i*xi) .* cgam(ct(xi)) .* cgam(ct(xi)-at(xi)-bt(xi)) ./ ( cgam(ct(xi)-at(xi)) .* cgam(ct(xi)-bt(xi)) );

         XI = [0.5,23];
         xi = linspace(XI(1),XI(end),16);

         digits(40);
         contspec_exact = vpa(rho(xi)).'
         ab_exact = vpa([a(xi),b(xi)]).'

         bound_states_exact = sym([0.5i; 1.5i])/d
         normconst_exact = [-1.0; 1.0].*exp_x0.^(-2i*bound_states_exact)

         p  = (ceil(delta-3/2):-1:0).';
         z0 = p+0.5-delta;
         Rres = exp_x0.^(-2i.*bound_states_exact) .* -1/(1i*d).*vpa(((-1).^p ./ factorial(p)) .* gamma(z0) .* ...
                            gamma(0.5-z0+delta) .* (cos(pi*delta)/pi) ./ gamma(-z0))

         */
        // Time domain
        T[0] = 0.0;
        T[1] = 16.0;

        // Time domain signal q(t) = 15*sech^2(2*(t-t0)), where t0 = ln(3000)
        for (i=0; i<D; i++)
         (*q_ptr)[i] = 15.0 * CPOW( misc_sech((T[0] + i*(T[1] - T[0])/(D - 1) - CLOG(3000.0)) * 2.0) , 2.0);

        // Nonlinear spectral domain
        XI[0] = 0.5;
        XI[1] = 23.0;

        // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
        (*contspec_ptr)[0] =    0.5318259355148260693807739748884874018125 +
                                   0.5358125960698921725004715585993891673643*I;
        (*contspec_ptr)[1] =    0.07773553930877936261392524113881729723477 +
                                   0.03740502735153968561056788472532282418146*I;
        (*contspec_ptr)[2] =    0.005925749089931136434807130666916924644125 +
                                   0.005655649605657166626330773355419327875735*I;
        (*contspec_ptr)[3] =    0.0001152057194103069939706988198438598137999 +
                                   0.000767811416093258581965862130680664842978*I;
        (*contspec_ptr)[4] =  - 0.00004642350465406691113674333056996134263829 +
                                   0.00005709683094928375907212684004446953418446*I;
        (*contspec_ptr)[5] =  - 0.000006970113555546589397701925994477822350578 -
                                   0.0000002524754626487246021393441167277101330574*I;
        (*contspec_ptr)[6] =  - 0.0000003487663246751790722843411115504658958195 -
                                   0.0000005615739163590524663784802994665295936717*I;
        (*contspec_ptr)[7] =    0.00000002767024065315693664690847131285015797807 -
                                   0.00000005621463172172617330262566591158045042654*I;
        (*contspec_ptr)[8] =    0.000000005920758288043957644193095203095124148974 -
                                   0.0000000004588787869964993694237201441588891286046*I;
        (*contspec_ptr)[9] =    0.0000000003160233962033068750514087762624735731768 +
                                   0.0000000004657611953394703064185926514558556423892*I;
        (*contspec_ptr)[10] = - 0.00000000002403499272223590921567101488246952893851 +
                                   0.00000000004762629684472271684727540672019450998716*I;
        (*contspec_ptr)[11] = - 0.000000000005053801049134480939172276097091203158252 +
                                   0.0000000000001582635686194822401672727116200476709589*I;
        (*contspec_ptr)[12] = - 0.0000000000002395360670397482632572131929638753991074 -
                                   0.0000000000004150769365323891568921366720902470409581*I;
        (*contspec_ptr)[13] =   0.0000000000000241442582749395751659610583728556629004 -
                                   0.00000000000003847355997786204330520518177525481704953*I;
        (*contspec_ptr)[14] =   0.000000000000004292166520809185825505953748638111392475 +
                                   0.0000000000000003335768060497820180155071798409472315344*I;
        (*contspec_ptr)[15] =   0.0000000000000001596518458599739742005323324672207694504 +
                                   0.0000000000000003755093964240090297007710985250023045271*I;

        // a(xi)
        (*ab_ptr)[0] =  - 0.8974787083260529020259013185313583011873 +
                                   1.232783952186963233458995331632565524713*I;
        (*ab_ptr)[1] =  - 0.8746487943185569432414862291030066823274 -
                                   0.4924299286261570020748194765224024097114*I;
        (*ab_ptr)[2] =  - 0.2551420535073444464514296927354953859368 -
                                   0.9669382803000529194212917616330172355992*I;
        (*ab_ptr)[3] =    0.1946202803764687530813509097411114432936 -
                                   0.9808789676985007795823511030326331818302*I;
        (*ab_ptr)[4] =    0.4599314446777565035493357048460048662902 -
                                   0.8879544310446464160752187321583069830986*I;
        (*ab_ptr)[5] =    0.618797471991278640583770720279439515216 -
                                   0.7855505640669153979199478219841801269779*I;
        (*ab_ptr)[6] =    0.7187583726886900957673738857583630023329 -
                                   0.6952599526007111209859036141344864743939*I;
        (*ab_ptr)[7] =    0.7848724411843467271273124184270709601577 -
                                   0.6196573658638523061250661154456663712045*I;
        (*ab_ptr)[8] =    0.8305471621414685024339491717013658365869 -
                                   0.5569483023214571796331531297678237849779*I;
        (*ab_ptr)[9] =    0.8632835298749002160160776783689998775782 -
                                   0.5047192754856230637368469949683900918432*I;
        (*ab_ptr)[10] =   0.8874829691804224040542148510838200241468 -
                                   0.4608405140769433449986007054866931649822*I;
        (*ab_ptr)[11] =   0.9058445052745624266075272372055025667922 -
                                   0.4236103542925777494041737014299986018305*I;
        (*ab_ptr)[12] =   0.9200891746783698178171576064444141792562 -
                                   0.3917089616535167718934216551055520811108*I;
        (*ab_ptr)[13] =   0.9313525976988554180286056893137409736291 -
                                   0.3641185778830763991957105251953477221674*I;
        (*ab_ptr)[14] =   0.9404068639236524338882758193695576497131 -
                                   0.3400513641867667238670562437482301985175*I;
        (*ab_ptr)[15] =   0.9477908531373028791171512608174271436213 -
                                   0.3188926131306016394435730788166750061641*I;
        // b(xi)
        (*ab_ptr)[16] = - 1.137843623474739401358106948560648396313 +
                                   0.1747460820338603763794070223910073422204*I;
        (*ab_ptr)[17] = - 0.04957194078314852104247025513288483673436 -
                                   0.07099556814801487478550869024948262636644*I;
        (*ab_ptr)[18] =   0.003956756312299505472298023383445974853865 -
                                   0.007172827688812989022620138530793675510022*I;
        (*ab_ptr)[19] =   0.0007755514386172862618427037328034335547693 +
                                   0.00003642880594817843953355506521399556742191*I;
        (*ab_ptr)[20] =   0.00002934775447747407190698088049832900472719 +
                                   0.00006748258460722619760253526832839226985268*I;
        (*ab_ptr)[21] = - 0.000004511420889761133273169500530504836609913 +
                                   0.000005319145457143217273183297001742627481681*I;
        (*ab_ptr)[22] = - 0.0000006411185704417375816968217374509248949438 -
                                   0.0000001611526959042578851889803231294513611763*I;
        (*ab_ptr)[23] = - 0.00000001311620128608975552564682272170253207002 -
                                   0.00000006126738366566434524154866781149026585355*I;
        (*ab_ptr)[24] =   0.000000004661897232471458476509763590788974789788 -
                                   0.000000003678676751288641321862253482237959685639*I;
        (*ab_ptr)[25] =   0.00000000050789644605850010875885947393388466041 +
                                   0.0000000002425808692231718010749271679615519788199*I;
        (*ab_ptr)[26] =   0.0000000000006174804161433534109945919861556925174987 +
                                   0.00000000005334382573677348807156385944889332248486*I;
        (*ab_ptr)[27] = - 0.000000000004510915824734781882645143090281769904194 +
                                   0.000000000002284204636967160253785584862235855188206*I;
        (*ab_ptr)[28] = - 0.000000000000382983898043729482010598783782146923735 -
                                   0.0000000000002880793718634050542411231818992040113776*I;
        (*ab_ptr)[29] =   0.000000000000008477879718638687538523755563524632996168 -
                                   0.00000000000004462382301521722991050587824787961640353*I;
        (*ab_ptr)[30] =   0.000000000000004149816105230553591700096706680754101857 -
                                   0.000000000000001145859162662987944883827205223397109093*I;
        (*ab_ptr)[31] =   0.0000000000000002710637318733171478745802860499448999733 +
                                   0.0000000000000003049925768803740636500487501422711915364*I;

        // Discrete spectrum
        (*bound_states_ptr)[0] = 1.0*I;
        (*bound_states_ptr)[1] = 3.0*I;
        (*normconsts_ptr)[0] = -9000000.0;
        (*normconsts_ptr)[1] = 729000000000000000000.0;
        (*residues_ptr)[0] = 22918311.80523292835071926192564206813296*I;
        (*residues_ptr)[1] = 7425533024895468785633.04086390803007508*I;

        break;

    case kdvv_testcases_SECH_SQUARED_LOW_BANDWIDTH:

    // q(t) = 15*sech^2(2*(t-t0)), where t0 = ln(3000).
    //
    // This test case can be calculated using §2.5 in
    // G. L. Lamb, Elements of soliton theory. Wiley, 1980
    // together with the time-shift as in Appendix B in e.g.
    // P. J. Prins and S. Wahls, “Soliton phase shift calculation for the
    // Korteweg–de Vries equation,” IEEE Access, vol. 7, pp. 122914–122930,
    // September 2019, DOI: 10.1109/access.2019.2932256 .
    // The following MATLAB code has been used to compute the values below.
    /*
     A = sym(15);
     d = sym(0.5);
     exp_x0 = sym(3000);
     delta = sqrt(sym(A)*d^2+1/4);
     at = @(xi) 1/2 - 1i*xi*d + delta;
     bt = @(xi) 1/2 - 1i*xi*d - delta;
     ct = @(xi) 1 - 1i*xi*d;
     cgam = @(z) gamma(sym(z));

     rho = @(xi) exp_x0.^(-2i*xi) .* cgam(ct(xi)-at(xi)-bt(xi)) .* cgam(at(xi)) .* cgam(bt(xi)) ./ ...
               ( cgam(ct(xi)-at(xi)) .* cgam(ct(xi)-bt(xi)) .* cgam(at(xi)+bt(xi)-ct(xi)) );
     a = @(xi) cgam(ct(xi)) .* cgam(at(xi)+bt(xi)-ct(xi)) ./ ( cgam(at(xi)) .* cgam(bt(xi)) );
     b = @(xi) exp_x0.^(-2i*xi) .* cgam(ct(xi)) .* cgam(ct(xi)-at(xi)-bt(xi)) ./ ( cgam(ct(xi)-at(xi)) .* cgam(ct(xi)-bt(xi)) );

     XI = [0.125,13.5];
     xi = linspace(XI(1),XI(end),16);

     digits(40);
     contspec_exact = vpa(rho(xi)).'
     ab_exact = vpa([a(xi),b(xi)]).'

     bound_states_exact = sym([0.5i; 1.5i])/d
     normconst_exact = [-1.0; 1.0].*exp_x0.^(-2i*bound_states_exact)

     p  = (ceil(delta-3/2):-1:0).';
     z0 = p+0.5-delta;
     Rres = exp_x0.^(-2i.*bound_states_exact) .* -1/(1i*d).*vpa(((-1).^p ./ factorial(p)) .* gamma(z0) .* ...
                        gamma(0.5-z0+delta) .* (cos(pi*delta)/pi) ./ gamma(-z0))

     */
    // Time domain
    T[0] = 0.0;
    T[1] = 16.0;

    // Time domain signal q(t) = 15*sech^2(2*(t-t0)), where t0 = ln(3000)
    for (i=0; i<D; i++)
     (*q_ptr)[i] = 15.0 * CPOW( misc_sech((T[0] + i*(T[1] - T[0])/(D - 1) - CLOG(3000.0)) * 2.0) , 2.0);

    // Nonlinear spectral domain
    XI[0] = 0.125;
    XI[1] = 13.25;

    // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
    (*contspec_ptr)[0] =    0.5463487523200426134905359052576073032297 +
            0.8148126005042640014940725998115521342038*I;
    (*contspec_ptr)[1] =    0.02591402890211170537565485727704698716401 -
            0.3976934199683518073194559775101118157078*I;
    (*contspec_ptr)[2] =  - 0.08658353191165942665494419253910660327494 +
            0.05920204186878143544152662670198546707553*I;
    (*contspec_ptr)[3] =    0.02275024160399274596833747104071103133584 +
            0.01378994267419491371534774547518038244945*I;
    (*contspec_ptr)[4] =    0.00188354818753865005533004952062200267558 -
            0.006462305433342114296226431070999710676313*I;
    (*contspec_ptr)[5] =  - 0.001688555437343426263439521546978832782004 -
            0.00022037544561397223732405718390865286028*I;
    (*contspec_ptr)[6] =  - 0.00002023711068242244350441617342131483541116 +
            0.0004303172915834161572084237912659988967271*I;
    (*contspec_ptr)[7] =    0.0001089731659069696561658595457552102346889 +
            0.000001370046692437382272197572238954177591503*I;
    (*contspec_ptr)[8] =    0.0000003338953643775709863111373720318524196761 -
            0.00002756813479639915751071964471727677930441*I;
    (*contspec_ptr)[9] =  - 0.000006970113555546589397701925994477822350578 -
            0.0000002524754626487246021393441167277101330574*I;
    (*contspec_ptr)[10] =  - 0.0000001384226387093051507417077552717173844069 +
            0.000001759014252193559897313539725965130709534*I;
    (*contspec_ptr)[11] =    0.0000004423264911123688023761484389977161211358 +
            0.00000005994749298178247431778371818665221756651*I;
    (*contspec_ptr)[12] =    0.00000002261398362324717933238629330708820605174 -
            0.0000001106350146572354380888048821206155377243*I;
    (*contspec_ptr)[13] =  - 0.00000002747792691033381289404818302909082750658 -
            0.000000007812938379155686582857298166787443823349*I;
    (*contspec_ptr)[14] =  - 0.000000002540631619899835897782925007281190829989 +
            0.000000006765581779882701740401254033815618727872*I;
    (*contspec_ptr)[15] =    0.000000001648647408135237716970547024517968767436 +
            0.0000000007902380949367284693755516059666318294607*I;

    // a(xi)
    (*ab_ptr)[0] =  - 0.8211421602571963596790751486473746028336 +
            5.092462056528973394858408181097639796561*I;
    (*ab_ptr)[1] =  - 1.016608586935804951271119396848111307515 +
            0.3941187209295174584677472336597360899871*I;
    (*ab_ptr)[2] =  - 0.9177460252702432554328005704799745414807 -
            0.4109333330690440164030971442772687506617*I;
    (*ab_ptr)[3] =  - 0.5629006270498119948798558227393823717523 -
            0.8269529136075027537457346439236290155814*I;
    (*ab_ptr)[4] =  - 0.2087685683256668103973429561680507655919 -
            0.9779882392327382372884171301030705201432*I;
    (*ab_ptr)[5] =    0.0694892742286887770604657797846631683028 -
            0.9975841521197000639809715327902189293127*I;
    (*ab_ptr)[6] =    0.2747723449112867398798379384353893033374 -
            0.961509409238671673935429644662582704443*I;
    (*ab_ptr)[7] =    0.4248522377700659575427275019695541479854 -
            0.9052627176343867293037919824972118229526*I;
    (*ab_ptr)[8] =    0.535654786091352804579938584328976875112 -
            0.8444370615371761457850807516514664629529*I;
    (*ab_ptr)[9] =    0.618797471991278640583770720279439515216 -
            0.7855505640669153979199478219841801269779*I;
    (*ab_ptr)[10] =    0.682304277940342490985349468533297306344 -
            0.7310683089201864419672492124900350459359*I;
    (*ab_ptr)[11] =    0.7316628535001429819410386562432192890068 -
            0.6816666845374119902989745713558235638208*I;
    (*ab_ptr)[12] =    0.7706525135042528397144070342813094911036 -
            0.6372556029016850833465318339541438113779*I;
    (*ab_ptr)[13] =    0.8019118460262724240317566305889773277387 -
            0.5974423747967135649928215643329180764571*I;
    (*ab_ptr)[14] =    0.8273128475993920381522283875768398655553 -
            0.5617414460381085106558769039895022147831*I;
    (*ab_ptr)[15] =    0.8482057730401525501411755029031228604423 -
            0.5296668448953145723790859553084161702063*I;

    // b(xi)
    (*ab_ptr)[16] =  - 4.59803224598356879442302388203668344958 +
            2.113183311838907665727593112064401272179*I;
    (*ab_ptr)[17] =    0.1303939976960228420187873229649100495649 +
            0.4145117496327245799873257063940069636971*I;
    (*ab_ptr)[18] =    0.1037897846554161561614265498765134159268 -
            0.01875237925560854193941168061445119760949*I;
    (*ab_ptr)[19] =  - 0.001402491991516315352062753119407495278761 -
            0.02657574595798170633665269753866140883135*I;
    (*ab_ptr)[20] =  - 0.006713284370623260367648399680408722414467 -
            0.0004929617220389287430595606758583274641336*I;
    (*ab_ptr)[21] =  - 0.0003371795438967164263461031217814111786488 +
            0.001669162414495812409493850262995881110576*I;
    (*ab_ptr)[22] =    0.0004081935264591172268058087565891677632852 +
            0.0001376974636012027862290614500887269862875*I;
    (*ab_ptr)[23] =    0.00004753774558454658513210540338526266716461 -
            0.00009806727691503475132000646668526600921709*I;
    (*ab_ptr)[24] =  - 0.0000231007020895495202597970112327915235835 -
            0.00001504895696765865148060497800981974717667*I;
    (*ab_ptr)[25] =  - 0.000004511420889761133273169500530504836609913 +
            0.000005319145457143217273183297001742627481681*I;
    (*ab_ptr)[26] =    0.000001191513216162502835779466719364916725139 +
            0.000001301379353627180049189191561584115659266*I;
    (*ab_ptr)[27] =    0.0000003644980714532028249854747303313484259542 -
            0.0000002576578789044047127669131943569054200415*I;
    (*ab_ptr)[28] =  - 0.00000005307525964783388573532654691414202142736 -
            0.00000009967203989501955847693593811132300889078*I;
    (*ab_ptr)[29] =  - 0.00000002670265555302393561817526314585218805089 +
            0.00000001015119006928210603614850433421785633307*I;
    (*ab_ptr)[30] =    0.00000000169861051215999937511823070542421789204 +
            0.000000007024430807993997034775589183628732809955*I;
    (*ab_ptr)[31] =    0.000000001816955167749214224087066872980938747342 -
            0.0000000002029493568102441037678330613882144065464*I;
   
    // Discrete spectrum
    (*bound_states_ptr)[0] = 1.0*I;
    (*bound_states_ptr)[1] = 3.0*I;
    (*normconsts_ptr)[0] = -9000000.0;
    (*normconsts_ptr)[1] = 729000000000000000000.0;
    (*residues_ptr)[0] = 22918311.80523292835071926192564206813296*I;
    (*residues_ptr)[1] = 7425533024895468785633.04086390803007508*I;

    break;
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
         zeta = sym(0.5:15.5).'*sym(pi)/sym(32);

         ell = sym(1)/sym(2);

         kappa = sqrt(ampl + zeta.^2);
         gamma = (kappa./zeta - zeta./kappa)/2;
         delta = (kappa./zeta + zeta./kappa)/2;
         T = exp(-2i*zeta*ell) ./ (cos(2*kappa*ell) - 1i*delta.*sin(2*kappa*ell));
         R = 1i*gamma.*sin(2*kappa*ell).*T;

         digits(40)
         T = vpa(T)
         R = vpa(R)
         */
        
        // Space domain
        T[0] = -1.0; // First midpoint, not Lower boundary
        T[1] = 2.0; // Last midpoint, not upper boundary
        
        // Space domain signal.
        eps_t = (T[1] - T[0])/(D-1);
        for (i=0; i<D; i++)
        {
            if (FABS(T[0] + i*eps_t) == 0.5)
                (*q_ptr)[i] = 0.5;
            else if (FABS(T[0] + i*eps_t) < 0.5)
                (*q_ptr)[i] = 1.0;
            else
                (*q_ptr)[i] = 0.0;
        }
        
        // Nonlinear spectral domain
        XI[0] = 0.5/32.0 * PI;
        XI[1] = 15.5/32.0 * PI;
        
        // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
        (*contspec_ptr)[0]  =  - 0.9870725601151357132059615656352921926505 + 0.1106668024779069100221269716635107304877*I;
        (*contspec_ptr)[1]  =  - 0.8942815149074410593957141404188288844772 + 0.3006083214632606321259647673494588361592*I;
        (*contspec_ptr)[2]  =  - 0.7517029093510107161488576830732266002746 + 0.4206295597197361344635948717938263256352*I;
        (*contspec_ptr)[3]  =   - 0.6048948399952833243051638927785187704323 + 0.473028823684821169748255791803965887681*I;
        (*contspec_ptr)[4]  =  - 0.4778504945251322156446773224766426668616 + 0.4793238390249636904293539788609397958394*I;
        (*contspec_ptr)[5]  =  - 0.3763474370619386743068998268064633149192 + 0.4600799580888626688056742042252798497312*I;
        (*contspec_ptr)[6]  =  - 0.2978533037542501812950159523446268231913 + 0.4288876704419926898290194462483566618216*I;
        (*contspec_ptr)[7]  =  - 0.2376683134714269108688326192379982136498 + 0.3933818871480632195012270772991049487594*I;
        (*contspec_ptr)[8]  =  - 0.1913637415605259149273312618011148544061 + 0.3574730594124660111171758774850972186401*I;
        (*contspec_ptr)[9]  =   - 0.1554156420268148907718456860386422438162 + 0.323011103360301963974093270371986465161*I;
        (*contspec_ptr)[10]  =  - 0.1271905034953081897751387070041458261916 + 0.2907720159569017256255551154712598493094*I;
        (*contspec_ptr)[11]  =   - 0.104764644824235261140851729086816610633 + 0.2609937738957654296620696476801746084408*I;
        (*contspec_ptr)[12]  = - 0.08673956215044220126053393160232166218066 + 0.2336564090694619140641045958452137680397*I;
        (*contspec_ptr)[13]  = - 0.07209443935525737929778535768652561998093 + 0.2086253871411293640983513028978855291533*I;
        (*contspec_ptr)[14]  = - 0.06007786636867152838650839313456544782936 + 0.1857236627095771710417966137341333503595*I;
        (*contspec_ptr)[15]  = - 0.05013091571056271949743226122945714631744 + 0.1647668721023619467965349247926784741708*I;
            // Dummy values to prevent this abolished testcase from causing undefined behavior
            for (UINT j=0; j<(*M_ptr)*2; j++)
                (*ab_ptr)[j]=0.0; // Not correct

            // Discrete spectrum
            for (UINT j=0; j<(*K_ptr); j++) {
                (*bound_states_ptr)[j] = 0; // Not correct
                (*normconsts_ptr)[j] = 0; // Not correct
                (*residues_ptr)[j] = 0; // Not correct
            }

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
             zeta = sym(0.5:15.5).'*sym(pi)/sym(32);

             ell = sym(1)/sym(2);

             kappa = sqrt(ampl + zeta.^2);
             gamma = (kappa./zeta - zeta./kappa)/2;
             delta = (kappa./zeta + zeta./kappa)/2;
             T = exp(-2i*zeta*ell) ./ (cos(2*kappa*ell) - 1i*delta.*sin(2*kappa*ell));
             R = 1i*gamma.*sin(2*kappa*ell).*T;

             digits(40)
             T = vpa(T)
             R = vpa(R)
             */
            
            // Space domain
            T[0] = -1.0; // First midpoint, not Lower boundary
            T[1] = 2.0; // Last midpoint, not upper boundary
            
            // Space domain signal.
            eps_t = (T[1] - T[0])/(D-1);
            for (i=0; i<D; i++)
            {
                if (FABS(T[0] + i*eps_t) == 0.5)
                    (*q_ptr)[i] = -0.5;
                else if (FABS(T[0] + i*eps_t) < 0.5)
                    (*q_ptr)[i] = -1.0;
                else
                    (*q_ptr)[i] = 0.0;
            }
            
            // Nonlinear spectral domain
            XI[0] = 0.5/32.0 * PI;
            XI[1] = 15.5/32.0 * PI;
            
            // Reflection coeffcient b(xi)/a(xi) (correct up to 40 digits)
            (*contspec_ptr)[0] = - 0.9933662141446059100932315762242385330732 - 0.07929705518519969139138112466905044336284*I;
            (*contspec_ptr)[1] = - 0.9431025368383791967365991561582619942712 - 0.2260135123901980111562861182523685355280*I;
            (*contspec_ptr)[2] = - 0.8555387014105028819357912730287807471101 - 0.3421954512691229543848328998205167762542*I;
            (*contspec_ptr)[3] = - 0.7492161643877380978254687870078062166628 - 0.4204150825082221655575943328825895033868*I;
            (*contspec_ptr)[4] = - 0.6403010582986216694290281017953758878661 - 0.4632317581880771713751550430198086335052*I;
            (*contspec_ptr)[5] = - 0.5390154224547847961678452502802902629005 - 0.4782368346748772947825999153625726028036*I;
            (*contspec_ptr)[6] = - 0.4500385219902367771505452893081657054438 - 0.4737894786049764378126476704768496389521*I;
            (*contspec_ptr)[7] = - 0.3744039392634054100382942060902825457698 - 0.4568958244732751736627310527747887949698*I;
            (*contspec_ptr)[8] = - 0.3112559152447967472679146094509503973361 - 0.4326903395570455121530570511284881081973*I;
            (*contspec_ptr)[9] = - 0.2589909242629944188605164103494425958224 - 0.4046560394782979749305520979914802557298*I;
            (*contspec_ptr)[10] = - 0.2158661747246824635078746485207390854714 - 0.3750415455056625172571914298026995862027*I;
            (*contspec_ptr)[11] = - 0.1802728439227340555723550674162770421086 - 0.3452458366296971023135761237543964355674*I;
            (*contspec_ptr)[12] = - 0.1508305280229313467723253327615386963076 - 0.3161092207254913210278840278787325221700*I;
            (*contspec_ptr)[13] = - 0.1263975284135399581284882178878508006411 - 0.2881146080008311951560139177538217304896*I;
            (*contspec_ptr)[14] = - 0.1060469774237907160964022087505637727953 - 0.2615204776144975185370691907861607140949*I;
            (*contspec_ptr)[15] = - 0.08903283712126179155781066343603991582133 - 0.2364465375984681127413182701998077329121*I;
            
            break;
    default: // unknown test case

        ret_code = E_INVALID_ARGUMENT(tc);
        goto release_mem_1;
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
release_mem_1:
    free(*q_ptr);

    return ret_code;
}

// Compares computed with exact nonlinear Fourier spectrum.
static INT kdvv_compare_nfs(const UINT M, const UINT K1, const UINT K2,
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

INT kdvv_testcases_test_fnft(kdvv_testcases_t tc, UINT D,
const REAL error_bounds[6], fnft_kdvv_opts_t * const opts) {
    COMPLEX * q = NULL;
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
    REAL errs[6] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN };
    INT ret_code;

    // Check inputs
    if (opts == NULL)
        return E_INVALID_ARGUMENT(opts);

    // Load test case
    ret_code = kdvv_testcases(tc, D, &q, T, &M, &contspec_exact, &ab_exact,
        XI, &K_exact,&bound_states_exact, &normconsts_exact, &residues_exact);
    CHECK_RETCODE(ret_code, release_mem);

    // Allocate memory
    contspec = malloc(3*M * sizeof(COMPLEX));
    K = kdv_discretization_degree(opts->discretization) * D;
    if (K == 0) //slow discretization
        K = K_exact;
    bound_states = malloc(K * sizeof(COMPLEX));
    normconsts_and_residues = malloc(2*D * sizeof(COMPLEX));
    if ( q == NULL || contspec == NULL || bound_states == NULL
            || normconsts_and_residues == NULL) {
        ret_code = E_NOMEM;
        goto release_mem;
    }

    if (opts->bound_state_localization == kdvv_bsloc_NEWTON){
        memcpy(bound_states,bound_states_exact,K * sizeof(COMPLEX));
    }

// Compute the NFT
    opts->contspec_type = fnft_kdvv_cstype_BOTH;
    opts->discspec_type = fnft_kdvv_dstype_BOTH;
    ret_code = fnft_kdvv(D, q, T, M, contspec, XI, &K, bound_states,
        normconsts_and_residues, opts);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the errors
    ret_code = kdvv_compare_nfs(M, K, K_exact, contspec, contspec_exact,
        contspec+M, ab_exact, bound_states, bound_states_exact,
        normconsts_and_residues, normconsts_exact, normconsts_and_residues+K,
        residues_exact, errs);
    CHECK_RETCODE(ret_code, release_mem);

#ifdef DEBUG
    for (UINT i=0; i<6; i++)
        printf("kdvv_testcases_test_fnft: error_bounds[%i] = %2.1e <= %2.1e\n",
            (int)i, errs[i], error_bounds[i]);
    misc_print_buf(M, contspec, "r_num");
    misc_print_buf(M, contspec_exact, "r_exact");
    misc_print_buf(2*M, contspec+M, "ab_num");
    misc_print_buf(2*M, ab_exact, "ab_exact");
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
    free(q);
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
