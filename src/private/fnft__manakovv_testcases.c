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
#include "fnft__manakovv_testcases.h"
#include "fnft__manakov_discretization.h"
#include "fnft_manakovv.h"
#ifdef DEBUG
#include <stdio.h>
#include <string.h> // for memcpy
#endif

INT manakovv_testcases(manakovv_testcases_t tc, const UINT D,
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
        
        case manakovv_testcases_SECH_FOCUSING:
            *M_ptr = 16;
            *K_ptr = 3;
            break;
            
        case manakovv_testcases_SECH_DEFOCUSING:
            *M_ptr = 16;
            *K_ptr = 0;
            break;

        case manakovv_testcases_RECTANGLE_FOCUSING:
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
        *contspec_ptr = malloc(2*(*M_ptr) * sizeof(COMPLEX));
        if (*contspec_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_2;
        }
        *ab_ptr = malloc(3*(*M_ptr) * sizeof(COMPLEX));
        if (*ab_ptr == NULL) {
            ret_code = E_NOMEM;
            goto release_mem_3;
        }
    } else {
        *contspec_ptr = NULL;
        *ab_ptr = NULL;
    }
    // generate test case
    switch (tc) {

    case manakovv_testcases_SECH_FOCUSING:

        // This test case can be found in J. Satsuma & N. Yajima, Prog. Theor.
        // Phys. Suppl., No. 55, p. 284ff, 1974.

        // The following MATLAB code has been used to compute the values below.
        /*
%% File to determine the exact NFT spectrum of focusing sech potential
% potential function: q = [A1*sech(t); A2*sech(t)]

A1 = sym(0.8);
A2 = sym(5.2);
q0 = sqrt(A1*A1+A2*A2);
hlf = sym(1)/sym(2);
im = sym(1j);
syms lam;

a(lam) = (gamma(hlf-im*lam)^2)/(gamma(hlf-im*lam+q0)*gamma(hlf-im*lam-q0));
b1(lam) = -sin(pi*q0)*A1*sech(lam*pi)/q0;
b2(lam) = -sin(pi*q0)*A2*sech(lam*pi)/q0;

q1(lam) = b1/a;
q2(lam) = b2/a;

xi = (-sym(7):sym(8))/4;        % xi has 16 elements, M = 16
XI = [xi(1) xi(end)]
digits(40);
contspec = vpa([(b1(xi)./a(xi)) (b2(xi)./a(xi))]).'
ab = vpa([a(xi) b1(xi) b2(xi)])'
        */

        // Time domain. sech(25) = sech(-25) = 2.7*10^-11. Outside these boundaries the potential is smaller still, so T domain sufficiently large
        T[0] = -25.0;
        T[1] = 25.0;

        // Test signal: q1 = 0.8*sech(t), q2 = 5.2*sech(t)
        for (i=0; i<D; i++){
            (*q1_ptr)[i] = 0.8 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
            (*q2_ptr)[i] = 5.2 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
        }
        // Nonlinear spectral domain
        XI[0] = -7.0 / 4.0;
        XI[1] = 8.0 / 4.0;

        // Reflection coeffcient [b1(xi)/a(xi) b2(xi)/a(xi)] (correct up to 40 digits)
(*contspec_ptr)[0] =        - 0.0008787570611436085549471989461347153962164 - 0.0002408088150114280992782406033839806080925*I;
(*contspec_ptr)[1] = - 0.001875609592410681149225912268628217023662 + 0.0006897800621393672508546839490768480441995*I;
(*contspec_ptr)[2] =  - 0.002167561273750198024857851801551035092841 + 0.003809743560787260383696824194817150778036*I;
(*contspec_ptr)[3] =    0.002707950757807773719699063337839507741654 + 0.009225347158596330624400431452355098514198*I;
(*contspec_ptr)[4] =     0.01981270952499152857106204097107883830657 + 0.007248964492162440251991514823653416797679*I;
(*contspec_ptr)[5] =      0.03377817869480677821650886139496229783776 - 0.03172732676694056093975135345472212270886*I;
(*contspec_ptr)[6] =    - 0.04362974545087608958418332657822128422822 - 0.09078327397334115453654423637669473804421*I;
(*contspec_ptr)[7] =                                                     -0.1631306710313506786173387647314438294433;
(*contspec_ptr)[8] =    - 0.04362974545087608958418332657822128422822 + 0.09078327397334115453654423637669473804421*I;
(*contspec_ptr)[9] =      0.03377817869480677821650886139496229783776 + 0.03172732676694056093975135345472212270886*I;
(*contspec_ptr)[10] =     0.01981270952499152857106204097107883830657 - 0.007248964492162440251991514823653416797679*I;
(*contspec_ptr)[11] =    0.002707950757807773719699063337839507741654 - 0.009225347158596330624400431452355098514198*I;
(*contspec_ptr)[12] =  - 0.002167561273750198024857851801551035092841 - 0.003809743560787260383696824194817150778036*I;
(*contspec_ptr)[13] = - 0.001875609592410681149225912268628217023662 - 0.0006897800621393672508546839490768480441995*I;
(*contspec_ptr)[14] =- 0.0008787570611436085549471989461347153962164 + 0.0002408088150114280992782406033839806080925*I;
(*contspec_ptr)[15] = - 0.0002832750817195448967966498024149974424002 + 0.000303870243058843773130016649294114360446*I;
(*contspec_ptr)[16] =  - 0.005711920897433455607156793149875650075407 - 0.001565257297574282645308563921995873952601*I;
(*contspec_ptr)[17] =    - 0.0121914623506694274699684297460834106538 + 0.004483570403905887130555445668999512287297*I;
(*contspec_ptr)[18] =    - 0.01408914827937628716157603671008172810347 + 0.02476333314511719249402935726631148005723*I;
(*contspec_ptr)[19] =      0.01760167992575052917804391169595680032075 + 0.05996475653087614905860280444030814034228*I;
(*contspec_ptr)[20] =       0.1287826119124449357119032663120124489927 + 0.04711826919905586163794484635374720918491*I;
(*contspec_ptr)[21] =        0.2195581615162440584073075990672549359455 - 0.2062276239851136461083837974556937976076*I;
(*contspec_ptr)[22] =      - 0.2835933454306945822971916227584383474834 - 0.5900912808267175044875375364485157972873*I;
(*contspec_ptr)[23] =                                                      -1.060349361703779411012701970754384891382;
(*contspec_ptr)[24] =      - 0.2835933454306945822971916227584383474834 + 0.5900912808267175044875375364485157972873*I;
(*contspec_ptr)[25] =        0.2195581615162440584073075990672549359455 + 0.2062276239851136461083837974556937976076*I;
(*contspec_ptr)[26] =       0.1287826119124449357119032663120124489927 - 0.04711826919905586163794484635374720918491*I;
(*contspec_ptr)[27] =      0.01760167992575052917804391169595680032075 - 0.05996475653087614905860280444030814034228*I;
(*contspec_ptr)[28] =    - 0.01408914827937628716157603671008172810347 - 0.02476333314511719249402935726631148005723*I;
(*contspec_ptr)[29] =    - 0.0121914623506694274699684297460834106538 - 0.004483570403905887130555445668999512287297*I;
(*contspec_ptr)[30] =  - 0.005711920897433455607156793149875650075407 + 0.001565257297574282645308563921995873952601*I;
(*contspec_ptr)[31] =  - 0.001841288031177041829178223715697483375602 + 0.001975156579882484525345108220411743342899*I;

        // a(xi)
(*ab_ptr)[0] =    - 0.9644260197535661844055887685714734731134 + 0.264284973916233303521273290494627448225*I;
(*ab_ptr)[1] =- 0.9384622196166738803568342561341871279729 - 0.3451318071638965566800072247374635856464*I;
(*ab_ptr)[2] =- 0.4943100550496050454214875083017573470469 - 0.8688079880664351532328775593907366338957*I;
(*ab_ptr)[3] =  0.2810892903972032141830471188653041109586 - 0.9576046680320670360215962770028529755004*I;
(*ab_ptr)[4] =  0.9302059085550343591590609057736645286585 - 0.3403385888744578025475381557967291611616*I;
(*ab_ptr)[5] =  0.6972263717196523322046195661096254813211 + 0.6548940701020918256290523117187639335591*I;
(*ab_ptr)[6] =- 0.3611240215511071771662267139149164025175 + 0.7514144455355884276706649859308055662765*I;
(*ab_ptr)[7] =                                               -0.6818433382464592089246771552443181892922;
(*ab_ptr)[8] =- 0.3611240215511071771662267139149164025175 - 0.7514144455355884276706649859308055662765*I;
(*ab_ptr)[9] =  0.6972263717196523322046195661096254813211 - 0.6548940701020918256290523117187639335591*I;
(*ab_ptr)[10] =  0.9302059085550343591590609057736645286585 + 0.3403385888744578025475381557967291611616*I;
(*ab_ptr)[11] =  0.2810892903972032141830471188653041109586 + 0.9576046680320670360215962770028529755004*I;
(*ab_ptr)[12] =- 0.4943100550496050454214875083017573470469 + 0.8688079880664351532328775593907366338957*I;
(*ab_ptr)[13] =- 0.9384622196166738803568342561341871279729 + 0.3451318071638965566800072247374635856464*I;
(*ab_ptr)[14] = - 0.9644260197535661844055887685714734731134 - 0.264284973916233303521273290494627448225*I;
(*ab_ptr)[15] =- 0.6818818551989321768088818959213628711837 - 0.7314572246134086248857721886650956234335*I;



    // [b1(xi) b2(xi)]
    (*ab_ptr)[16] =    0.0009111383262031659175316826947484097435264;
    (*ab_ptr)[17] =    0.001998253780619837518165677887285859874399;
    (*ab_ptr)[18] =    0.004381382970647488704597228040424138126556;
    (*ab_ptr)[19] =    0.00959541146023096788944518188346771235739;
    (*ab_ptr)[20] =    0.02089700181069534360371965813183155655875;
    (*ab_ptr)[21] =    0.04432907513453893771132366465340876094225;
    (*ab_ptr)[22] =    0.08397161261305505669722345566748306488696;
    (*ab_ptr)[23] =    0.1112295613064011054907499242659422327261;
    (*ab_ptr)[24] =    0.08397161261305505669722345566748306488696;
    (*ab_ptr)[25] =    0.04432907513453893771132366465340876094225;
    (*ab_ptr)[26] =    0.02089700181069534360371965813183155655875;
    (*ab_ptr)[27] =    0.00959541146023096788944518188346771235739;
    (*ab_ptr)[28] =    0.004381382970647488704597228040424138126556;
    (*ab_ptr)[29] =    0.001998253780619837518165677887285859874399;
    (*ab_ptr)[30] =    0.0009111383262031659175316826947484097435264;
    (*ab_ptr)[31] =    0.0004154282228849761555074645878116921538315;
    (*ab_ptr)[32] =    0.005922399120320578463955937515864663332921;
    (*ab_ptr)[33] =    0.0129886495740289438680769062673580891836;
    (*ab_ptr)[34] =    0.02847898930920867657988198226275689782261;
    (*ab_ptr)[35] =    0.06237017449150129128139368224254013032304;
    (*ab_ptr)[36] =    0.1358305117695197334241777778569051176319;
    (*ab_ptr)[37] =    0.2881389883745030951236038202471569461246;
    (*ab_ptr)[38] =    0.5458154819848578685319524618386399217653;
    (*ab_ptr)[39] =    0.7229921484916071856898745077286245127199;
    (*ab_ptr)[40] =    0.5458154819848578685319524618386399217653;
    (*ab_ptr)[41] =    0.2881389883745030951236038202471569461246;
    (*ab_ptr)[42] =    0.1358305117695197334241777778569051176319;
    (*ab_ptr)[43] =    0.06237017449150129128139368224254013032304;
    (*ab_ptr)[44] =    0.02847898930920867657988198226275689782261;
    (*ab_ptr)[45] =    0.0129886495740289438680769062673580891836;
    (*ab_ptr)[46] =    0.005922399120320578463955937515864663332921;
    (*ab_ptr)[47] =    0.002700283448752345010798519820775998999905;

        *kappa_ptr = +1;

        break;

        case manakovv_testcases_SECH_DEFOCUSING:

        // This test case can be found in J. Satsuma & N. Yajima, Prog. Theor.
        // Phys. Suppl., No. 55, p. 284ff, 1974.

        // The following MATLAB code has been used to compute the values below.
        /*
%% File to determine the exact NFT spectrum of defocusing sech potential
% potential function: q = [A1*sech(t); A2*sech(t)]

A1 = sym(0.8);
A2 = sym(5.2);
kappa = sym(-1);
q0 = sqrt(kappa)*sqrt(A1*A1+A2*A2);
hlf = sym(1)/sym(2);
im = sym(1j);
syms lam;

a(lam) = (gamma(hlf-im*lam)^2)/(gamma(hlf-im*lam+q0)*gamma(hlf-im*lam-q0));
b1(lam) = -kappa*sin(pi*q0)*A1*sech(lam*pi)/q0;
b2(lam) = -kappa*sin(pi*q0)*A2*sech(lam*pi)/q0;

q1(lam) = b1/a;
q2(lam) = b2/a;

xi = (-sym(7):sym(8))/4;        % xi has 16 elements, M = 16
XI = [xi(1) xi(end)]
digits(40);
contspec = vpa([(b1(xi)./a(xi)) (b2(xi)./a(xi))]).'
ab = vpa([a(xi) b1(xi) b2(xi)]).'        
        */

        T[0] = -25.0;
        T[1] = 25.0;

        // Test signal: q1 = 0.8*sech(t), q2 = 5.2*sech(t)
        for (i=0; i<D; i++){
            (*q1_ptr)[i] = 0.8 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
            (*q2_ptr)[i] = 5.2 * misc_sech(T[0] + i*(T[1] - T[0])/(D - 1));
        }
        // Nonlinear spectral domain
        XI[0] = -7.0 / 4.0;
        XI[1] = 8.0 / 4.0;

        // Reflection coeffcient [b1(xi)/a(xi) b2(xi)/a(xi)] (correct up to 40 digits)

(*contspec_ptr)[0] =         0.08870047970469106375367392898889863286297 + 0.1235055147648429924158910678580640755161*I;
(*contspec_ptr)[1] =       0.1413477130935491832574730598175048629 + 0.05605543047020132455539097332292723879078*I;
(*contspec_ptr)[2] =    0.1457940027369356466391083250354616905268 - 0.04319138859676859532183126486832372904435*I;
(*contspec_ptr)[3] =    0.07330222302821077938980498086769721286531 - 0.1332222630880617735621072252436300722421*I;
(*contspec_ptr)[4] =   - 0.0624443181861419862532098918435537801528 - 0.1386437680153468700412003279141211073202*I;
(*contspec_ptr)[5] = - 0.1520347385343029112855061604914521340402 - 0.002612577662114768427800776456992842191349*I;
(*contspec_ptr)[6] =  - 0.02470913753366775732004032844169351861581 + 0.1500361483295894640417645987438715100841*I;
(*contspec_ptr)[7] =                                                   0.1520571842539397649653291997585641621625;
(*contspec_ptr)[8] =  - 0.02470913753366775732004032844169351861581 - 0.1500361483295894640417645987438715100841*I;
(*contspec_ptr)[9] = - 0.1520347385343029112855061604914521340402 + 0.002612577662114768427800776456992842191349*I;
(*contspec_ptr)[10] =   - 0.0624443181861419862532098918435537801528 + 0.1386437680153468700412003279141211073202*I;
(*contspec_ptr)[11] =    0.07330222302821077938980498086769721286531 + 0.1332222630880617735621072252436300722421*I;
(*contspec_ptr)[12] =    0.1457940027369356466391083250354616905268 + 0.04319138859676859532183126486832372904435*I;
(*contspec_ptr)[13] =       0.1413477130935491832574730598175048629 - 0.05605543047020132455539097332292723879078*I;
(*contspec_ptr)[14] =    0.08870047970469106375367392898889863286297 - 0.1235055147648429924158910678580640755161*I;
(*contspec_ptr)[15] =    0.02052893061506637005106835697814793007902 - 0.1506650266712882851191910167229956615927*I;
(*contspec_ptr)[16] =     0.5765531180804919143988805384278411136093 + 0.8027858459714794507032919410774164908549*I;
(*contspec_ptr)[17] =       0.91876013510806969117357488881378160885 + 0.3643602980563086096100413265990270521401*I;
(*contspec_ptr)[18] =     0.9476610177900817031542041127305009884242 - 0.2807440258789958695919032216441042387883*I;
(*contspec_ptr)[19] =     0.4764644496833700660337323756400318836245 - 0.8659447100724015281536969640835954695738*I;
(*contspec_ptr)[20] =   - 0.4058880682099229106458642969830995709932 - 0.9011844920997546552678021314417871975814*I;
(*contspec_ptr)[21] =  - 0.9882258004729689233557900431944388712611 - 0.01698175480374599478070504697045347424377*I;
(*contspec_ptr)[22] =   - 0.1606093939688404225802621348710078710027 + 0.9752349641423315162714698918351648155468*I;
(*contspec_ptr)[23] =                                                   0.9883716976506084722746397984306670540563;
(*contspec_ptr)[24] =   - 0.1606093939688404225802621348710078710027 - 0.9752349641423315162714698918351648155468*I;
(*contspec_ptr)[25] =  - 0.9882258004729689233557900431944388712611 + 0.01698175480374599478070504697045347424377*I;
(*contspec_ptr)[26] =   - 0.4058880682099229106458642969830995709932 + 0.9011844920997546552678021314417871975814*I;
(*contspec_ptr)[27] =     0.4764644496833700660337323756400318836245 + 0.8659447100724015281536969640835954695738*I;
(*contspec_ptr)[28] =     0.9476610177900817031542041127305009884242 + 0.2807440258789958695919032216441042387883*I;
(*contspec_ptr)[29] =       0.91876013510806969117357488881378160885 - 0.3643602980563086096100413265990270521401*I;
(*contspec_ptr)[30] =     0.5765531180804919143988805384278411136093 - 0.8027858459714794507032919410774164908549*I;
(*contspec_ptr)[31] =     0.1334380489979314053319443203579615455136 - 0.9793226733633738532747416086994718003526*I;


        // a(x*I;)
(*ab_ptr)[0] =         36014.77957890819387998363006316493457907 - 50146.55958845024284838306267800767775926*I;
(*ab_ptr)[1] =    125866.4278391712653335265082033866277088 - 49915.8892624013416961381812163906375796*I;
(*ab_ptr)[2] =   284656.6622388098211977172116830015852989 + 84329.37078762829164101515719576396535501*I;
(*ab_ptr)[3] =   313437.6952783077122955790197178919206757 + 569653.6527413646658944278771749083691513*I;
(*ab_ptr)[4] = - 581496.8755704972766884734952966497999795 + 1291084.925900222603502967865629026470843*I;
(*ab_ptr)[5] =  - 3003322.565660918256107034153266867516118 + 51609.3461455888194198302760800078263715*I;
(*ab_ptr)[6] = - 924614.1969422518563918640514757912506358 - 5614342.168400237966910808921658715176766*I;
(*ab_ptr)[7] =                                                7536984.908257509552988682238955806577404;
(*ab_ptr)[8] = - 924614.1969422518563918640514757912506358 + 5614342.168400237966910808921658715176766*I;
(*ab_ptr)[9] =  - 3003322.565660918256107034153266867516118 - 51609.3461455888194198302760800078263715*I;
(*ab_ptr)[10] = - 581496.8755704972766884734952966497999795 - 1291084.925900222603502967865629026470843*I;
(*ab_ptr)[11] =   313437.6952783077122955790197178919206757 - 569653.6527413646658944278771749083691513*I;
(*ab_ptr)[12] =   284656.6622388098211977172116830015852989 - 84329.37078762829164101515719576396535501*I;
(*ab_ptr)[13] =    125866.4278391712653335265082033866277088 + 49915.8892624013416961381812163906375796*I;
(*ab_ptr)[14] =   36014.77957890819387998363006316493457907 + 50146.55958845024284838306267800767775926*I;
(*ab_ptr)[15] =   3800.430566485191297940367385872825733879 + 27891.95323411732437427863237685178500824*I;

    // [b1(xi) b2(xi)]
(*ab_ptr)[16] =                                                 9387.90488076528882790077727523625179291;
(*ab_ptr)[17] =                                                 20588.9883902278989829955794266304810822;
(*ab_ptr)[18] =                                                45143.53681734144429975823944397546038068;
(*ab_ptr)[19] =                                                98866.22863932434613848755029231450495667;
(*ab_ptr)[20] =                                                 215312.054876993313140761968588644965231;
(*ab_ptr)[21] =                                                456744.1944293258085417043507294635441973;
(*ab_ptr)[22] =                                                865200.6937089949201452946795645976315595;
(*ab_ptr)[23] =                                                1146052.702914075425630703042320723335725;
(*ab_ptr)[24] =                                                865200.6937089949201452946795645976315595;
(*ab_ptr)[25] =                                                456744.1944293258085417043507294635441973;
(*ab_ptr)[26] =                                                 215312.054876993313140761968588644965231;
(*ab_ptr)[27] =                                                98866.22863932434613848755029231450495667;
(*ab_ptr)[28] =                                                45143.53681734144429975823944397546038068;
(*ab_ptr)[29] =                                                 20588.9883902278989829955794266304810822;
(*ab_ptr)[30] =                                                 9387.90488076528882790077727523625179291;
(*ab_ptr)[31] =                                                4280.360653339364090407615272124275346637;
(*ab_ptr)[32] =                                                61021.38172497437738135505228903563665391;
(*ab_ptr)[33] =                                                133828.4245364813433894712662730981270343;
(*ab_ptr)[34] =                                                293432.9893127193879484285563858404924744;
(*ab_ptr)[35] =                                                642630.4861556082499001690769000442822183;
(*ab_ptr)[36] =                                                1399528.356700456535414952795826192274001;
(*ab_ptr)[37] =                                                2968837.263790617755521078279741513037282;
(*ab_ptr)[38] =                                                5623804.509108466980944415417169884605137;
(*ab_ptr)[39] =                                                7449342.568941490266599569775084701682215;
(*ab_ptr)[40] =                                                5623804.509108466980944415417169884605137;
(*ab_ptr)[41] =                                                2968837.263790617755521078279741513037282;
(*ab_ptr)[42] =                                                1399528.356700456535414952795826192274001;
(*ab_ptr)[43] =                                                642630.4861556082499001690769000442822183;
(*ab_ptr)[44] =                                                293432.9893127193879484285563858404924744;
(*ab_ptr)[45] =                                                133828.4245364813433894712662730981270343;
(*ab_ptr)[46] =                                                61021.38172497437738135505228903563665391;
(*ab_ptr)[47] =                                                27822.34424670586658764949926880778975314;

        *kappa_ptr = -1;

        break;

    case manakovv_testcases_RECTANGLE_FOCUSING:

    // For this tescase the MZS matrix P is piecewise constant and the solution is easily found:
        // The following MATLAB code has been used to compute the values below.
        /*
%% File to determine the exact NFT spectrum of rectangle potential focusing manakov equation
% potential function: q = [A1; A2] for L(1)<t<L(2), q = [0; 0] otherwise

A1 = sym(0.8);
A2 = sym(5.2);
L = [-2 3];
kappa = sym(1);
XI_vector = (-sym(7):sym(8))/4;        % xi has 16 elements, M = 16

a=zeros(1,length(XI_vector));
b1=zeros(1,length(XI_vector));
b2=zeros(1,length(XI_vector));
for i=1:length(XI_vector)
    L1 = L(1); L2 = L(2);
    lam = XI_vector(i);
    P0 = [-1j*lam, 0, 0;
          0, 1j*lam, 0;
          0, 0, 1j*lam];
    PL = [-1j*lam, A1, A2;
          -kappa*conj(A1), 1j*lam, 0;
          -kappa*conj(A2), 0, 1j*lam];
    
      % For the first part -inf<t<L1, we have v = [exp(j*lam*t); 0; 0]
      % From this we get the BC at L1:
      v_L1 = [exp(-1j*lam*L1); 0; 0];
      
      % L1<t<L2, v=exp(PL*t)*c2 where we get the constant c2 from v_L1
      c2 = expm(PL*L1)\v_L1;
      v_L2 = expm(PL*L2)*c2;
      
      % L2<t, v = exp(P0*t)*c2 where we get c3 from v_L2
      c3 = expm(P0*L2)\v_L2;
%      v_fun = @(t)expm(PL*t)*c3;
%       v_end = expm(P0*t)*c3;
%       ab(1) = v_end(1)*exp(1j*lam*t);
%       ab(2) = v_end(2)*exp(-1j*lam*t);
%       ab(3) = v_end(3)*exp(-1j*lam*t);
% expm(P0*t)*[epx(j*lam*t); epx(-j*lam*t); epx(-j*lam*t)] cancel against
% each other (becomes [1; 1; 1], so we do not need ot specify at which t
% and we just get c3 as the answer
    a(i) = c3(1);
    b1(i) = c3(2);
    b2(i) = c3(3);
end
contspec = vpa([(b1./a) (b2./a)]).'
ab = vpa([a b1 b2]).'
        */

        // Time domain. q1,2(t) = 0 for L(1)<t<L(2), so choose T[0]<L(1) and T[1]>L(2)
        T[0] = -5.0;
        T[1] = 5.0;
        REAL L[2] = {-2.0, 3.0};
        REAL eps_t = (T[1]-T[0])/(D-1);

        // Test signal: q = [A1; A2] for L(1)<t<L(2), q = [0; 0] otherwise
        for (i=0; i<D; i++){
            if (T[0]+i*eps_t>=L[0] && T[0]+i*eps_t<=L[1]){
            (*q1_ptr)[i] = 0.8;
            (*q2_ptr)[i] = 5.2;
            }
            else{
            (*q1_ptr)[i] = 0.0;
            (*q2_ptr)[i] = 0.0;
            }
        }
        // Nonlinear spectral domain
        XI[0] = -7.0 / 4.0;
        XI[1] = 8.0 / 4.0;

        // Reflection coeffcient [b1(xi)/a(xi) b2(xi)/a(xi)] (correct up to 40 digits)
(*contspec_ptr)[0] =- 0.026062734927461455453068239762615 - 0.083118540185744507753362597668456*I;
(*contspec_ptr)[1] = - 0.18017636533205713678107429132069 + 0.014271311686288888631590765498913*I;
(*contspec_ptr)[2] = - 0.081717097538348423824139388216281 + 0.34476662876851466021932424155239*I;
(*contspec_ptr)[3] =     0.4552421404279109462542862729606 + 0.58906635331663270704893875517882*I;
(*contspec_ptr)[4] =     0.66651666775858309232916099063004 + 0.4007723270173131346894024318317*I;
(*contspec_ptr)[5] =    0.44932281393344897235664348045248 - 0.22286308263152837527698579833668*I;
(*contspec_ptr)[6] = - 0.075113760364707085614988102406642 - 0.38374043607692320811963782034582*I;
(*contspec_ptr)[7] =                                         -0.36202399287992614862474738401943;
(*contspec_ptr)[8] = - 0.075113760364707085614988102406642 + 0.38374043607692320811963782034582*I;
(*contspec_ptr)[9] =    0.44932281393344897235664348045248 + 0.22286308263152837527698579833668*I;
(*contspec_ptr)[10] =     0.66651666775858309232916099063004 - 0.4007723270173131346894024318317*I;
(*contspec_ptr)[11] =     0.4552421404279109462542862729606 - 0.58906635331663270704893875517882*I;
(*contspec_ptr)[12] = - 0.081717097538348423824139388216281 - 0.34476662876851466021932424155239*I;
(*contspec_ptr)[13] = - 0.18017636533205713678107429132069 - 0.014271311686288888631590765498913*I;
(*contspec_ptr)[14] =- 0.026062734927461455453068239762615 + 0.083118540185744507753362597668456*I;
(*contspec_ptr)[15] = 0.016345346761064319385248211347061 + 0.0093432109285461941577954902982128*I;
(*contspec_ptr)[16] =  - 0.16940777702849943442409141880489 - 0.54027051120733926570238736530882*I;
(*contspec_ptr)[17] =  - 1.1711463746583714584659219326568 + 0.092763525960877751819211312067637*I;
(*contspec_ptr)[18] =   - 0.53116113399926467852907308042631 + 2.2409830869953455412257881107507*I;
(*contspec_ptr)[19] =       2.9590739127814211784084363898728 + 3.828931296558113039907311758725*I;
(*contspec_ptr)[20] =      4.3323583404307894895168828952592 + 2.6050201256125355975257207319373*I;
(*contspec_ptr)[21] =       2.920598290567417709695519079105 - 1.4486100371049344115448320735595*I;
(*contspec_ptr)[22] =    - 0.4882394423705959662918019148492 - 2.4943128345000005197107384447008*I;
(*contspec_ptr)[23] =                                          -2.3531559537195194664604969148058;
(*contspec_ptr)[24] =    - 0.4882394423705959662918019148492 + 2.4943128345000005197107384447008*I;
(*contspec_ptr)[25] =       2.920598290567417709695519079105 + 1.4486100371049344115448320735595*I;
(*contspec_ptr)[26] =      4.3323583404307894895168828952592 - 2.6050201256125355975257207319373*I;
(*contspec_ptr)[27] =       2.9590739127814211784084363898728 - 3.828931296558113039907311758725*I;
(*contspec_ptr)[28] =   - 0.53116113399926467852907308042631 - 2.2409830869953455412257881107507*I;
(*contspec_ptr)[29] =  - 1.1711463746583714584659219326568 - 0.092763525960877751819211312067637*I;
(*contspec_ptr)[30] =  - 0.16940777702849943442409141880489 + 0.54027051120733926570238736530882*I;
(*contspec_ptr)[31] =   0.10624475394691809682079508547758 + 0.060730871035550262893032424926787*I;

        // a(xi)
(*ab_ptr)[0] =        0.76842151850790341161001606451464 + 0.40303701190105989393686058974708*I;
(*ab_ptr)[1] =- 0.0053085791801664762579182621493601 + 0.64375338202718679347924535250058*I;
(*ab_ptr)[2] =   - 0.33548292877868279848740939996787 + 0.2073160938127976138378727455347*I;
(*ab_ptr)[3] = - 0.19935631425978672659482526796637 - 0.017418064471645532120813726351116*I;
(*ab_ptr)[4] =  - 0.1877204314188311018440913358063 - 0.039742357177227530251606424371857*I;
(*ab_ptr)[5] =  - 0.16629080878327068271360644757806 - 0.23774563516814151142853006604128*I;
(*ab_ptr)[6] =    0.15545410698654801850970841314847 - 0.32739956195514313685990259727987*I;
(*ab_ptr)[7] =                                          0.38724787039322466908330966361973;
(*ab_ptr)[8] =    0.15545410698654801850970841314847 + 0.32739956195514313685990259727987*I;
(*ab_ptr)[9] =  - 0.16629080878327068271360644757806 + 0.23774563516814151142853006604128*I;
(*ab_ptr)[10] =  - 0.1877204314188311018440913358063 + 0.039742357177227530251606424371857*I;
(*ab_ptr)[11] = - 0.19935631425978672659482526796637 + 0.017418064471645532120813726351116*I;
(*ab_ptr)[12] =   - 0.33548292877868279848740939996787 - 0.2073160938127976138378727455347*I;
(*ab_ptr)[13] =- 0.0053085791801664762579182621493601 - 0.64375338202718679347924535250058*I;
(*ab_ptr)[14] =    0.76842151850790341161001606451464 - 0.40303701190105989393686058974708*I;
(*ab_ptr)[15] =    0.80637785454689814290674121366465 + 0.57849423805943522847172744150157*I;

        // [b1(x*I;) b2(x*I;)]
(*ab_ptr)[16] =  0.013472681720611728514724880767517 - 0.074374321672823442974653573855903*I;
(*ab_ptr)[17] = - 0.0082307246622527578672201542531184 - 0.1160649049319692482917787401675*I;
(*ab_ptr)[18] =  - 0.04406097953983694609503274364215 - 0.13260458782378395370393775465345*I;
(*ab_ptr)[19] = - 0.080494999491298332983113539285114 - 0.12536353400384045442272906711878*I;
(*ab_ptr)[20] =  - 0.10919115945241222864758867672208 - 0.10172209760305597758822671039525*I;
(*ab_ptr)[21] = - 0.12770297926953097511493240290292 - 0.069764455535414729436816116958653*I;
(*ab_ptr)[22] = - 0.13731319321595736315799740623333 - 0.035061814564757910750447678083219*I;
(*ab_ptr)[23] =                                         -0.14019302027400332155693263302965;
(*ab_ptr)[24] = - 0.13731319321595736315799740623333 + 0.035061814564757910750447678083219*I;
(*ab_ptr)[25] = - 0.12770297926953097511493240290292 + 0.069764455535414729436816116958653*I;
(*ab_ptr)[26] =  - 0.10919115945241222864758867672208 + 0.10172209760305597758822671039525*I;
(*ab_ptr)[27] = - 0.080494999491298332983113539285114 + 0.12536353400384045442272906711878*I;
(*ab_ptr)[28] =  - 0.04406097953983694609503274364215 + 0.13260458782378395370393775465345*I;
(*ab_ptr)[29] = - 0.0082307246622527578672201542531184 + 0.1160649049319692482917787401675*I;
(*ab_ptr)[30] =  0.013472681720611728514724880767517 + 0.074374321672823442974653573855903*I;
(*ab_ptr)[31] = 0.0077755319658742196031742111017593 + 0.016989847303499373831758845199147*I;
(*ab_ptr)[32] =   0.087572431183976229274179559070035 - 0.48343309087335234464077871052723*I;
(*ab_ptr)[33] = - 0.053499710304642920932760574714848 - 0.75442188205780014165213742671767*I;
(*ab_ptr)[34] =  - 0.28639636700894016696494759344205 - 0.86192982085459568519780759743298*I;
(*ab_ptr)[35] =  - 0.52321749669343919908470752488938 - 0.81486297102496296762552674408653*I;
(*ab_ptr)[36] =  - 0.70974253644067952784268982213689 - 0.66119363441986389595683704101248*I;
(*ab_ptr)[37] =  - 0.83006936525195129661369719542563 - 0.45346896098019573440041085632402*I;
(*ab_ptr)[38] =  - 0.89253575590372280501583190925885 - 0.22790179467092644416403857121622*I;
(*ab_ptr)[39] =                                         -0.91125463178102150685333526780596;
(*ab_ptr)[40] =  - 0.89253575590372280501583190925885 + 0.22790179467092644416403857121622*I;
(*ab_ptr)[41] =  - 0.83006936525195129661369719542563 + 0.45346896098019573440041085632402*I;
(*ab_ptr)[42] =  - 0.70974253644067952784268982213689 + 0.66119363441986389595683704101248*I;
(*ab_ptr)[43] =  - 0.52321749669343919908470752488938 + 0.81486297102496296762552674408653*I;
(*ab_ptr)[44] =  - 0.28639636700894016696494759344205 + 0.86192982085459568519780759743298*I;
(*ab_ptr)[45] = - 0.053499710304642920932760574714848 + 0.75442188205780014165213742671767*I;
(*ab_ptr)[46] =   0.087572431183976229274179559070035 + 0.48343309087335234464077871052723*I;
(*ab_ptr)[47] =    0.050540957778182425685908896184628 + 0.1104340074727459403147733496553*I;

        *kappa_ptr = +1;

        break;

    default: // unknown test case
        printf("unknown test case");
        ret_code = E_INVALID_ARGUMENT(tc);
        goto release_mem_5;
    }

    return SUCCESS;

    // the code below is only executed if an error occurs

release_mem_5:
    free(*normconsts_ptr);
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
static INT manakovv_compare_nfs(const UINT M, const UINT K1, const UINT K2,
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
    REAL dists[5])
{
    UINT i, j, min_j = 0;
    REAL dist, min_dist, nrm;

    // Check last argument
    if (dists == NULL)
        return E_INVALID_ARGUMENT(dists);

    if (contspec_1 == NULL || contspec_2 == NULL || M == 0){
        dists[0] = NAN;
        dists[1] = NAN;
    }
    else {
        // dist. refelction coef 1
        dists[0] = 0.0;
        nrm = 0.0;
        for (i=0; !(i>=M); i++) {
            dists[0] += CABS(contspec_1[i] - contspec_2[i]);
            nrm += CABS(contspec_2[i]);
        }
        if (nrm > 0)
            dists[0] /= nrm;
        // dist. refelction coef 2
        dists[1] = 0.0;
        nrm = 0.0;
        for (i=M; !(i>=2*M); i++) {
            dists[1] += CABS(contspec_1[i] - contspec_2[i]);
            nrm += CABS(contspec_2[i]);
        }
        if (nrm > 0)
            dists[1] /= nrm;
    }

    if (ab_1 == NULL || ab_2 == NULL || M == 0) {
        dists[2] = NAN; // error a
        dists[3] = NAN; // error b1
        dists[4] = NAN; // error b2
    }
    else {
        // error in a
        dists[2] = 0.0;
        nrm = 0.0;
        for (i=0; !(i>=M); i++) {
            dists[2] += CABS(ab_1[i] - ab_2[i]);
            nrm += CABS(ab_2[i]);
        }
        if (nrm > 0)
            dists[2] /= nrm;

        // error in b1
        dists[3] = 0.0;
        nrm = 0.0;
        for (i=M; !(i>=2*M); i++) {
            dists[3] += CABS(ab_1[i] - ab_2[i]);
            nrm += CABS(ab_2[i]);
        }
        if (nrm > 0)
            dists[3] /= nrm;

        // error in b2
        dists[4] = 0.0;
        nrm = 0.0;
        for (i=2*M; !(i>=3*M); i++) {
            dists[4] += CABS(ab_1[i] - ab_2[i]);
            nrm += CABS(ab_2[i]);
        }
        if (nrm > 0)
            dists[4] /= nrm;
    }

    return SUCCESS;
}

INT manakovv_testcases_test_fnft(manakovv_testcases_t tc, UINT D,
const REAL error_bounds[5], fnft_manakovv_opts_t * const opts) {
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
    REAL errs[5] = {
        FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN, FNFT_NAN};
    INT ret_code;

    // Check inputs
    if (opts == NULL)
        return E_INVALID_ARGUMENT(opts);

    // Load test case
    ret_code = manakovv_testcases(tc, D, &q1, &q2, T, &M, &contspec_exact, &ab_exact,
        XI, &K_exact,&bound_states_exact, &normconsts_exact, &residues_exact,
        &kappa);
    CHECK_RETCODE(ret_code, release_mem);

    // Allocate memory
    contspec = malloc(5*M * sizeof(COMPLEX));
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
    
    if (opts->bound_state_localization == manakovv_bsloc_NEWTON){
        memcpy(bound_states,bound_states_exact,K * sizeof(COMPLEX));
    }
    
// Compute the NFT
    opts->contspec_type = fnft_manakovv_cstype_BOTH;
    opts->discspec_type = fnft_manakovv_dstype_BOTH;

    ret_code = fnft_manakovv(D, q1, q2, T, M, contspec, XI, &K, bound_states,
        normconsts_and_residues, kappa, opts);
    CHECK_RETCODE(ret_code, release_mem);

    // Compute the errors
    ret_code = manakovv_compare_nfs(M, K, K_exact, contspec, contspec_exact,
        contspec+2*M, ab_exact, bound_states, bound_states_exact,
        normconsts_and_residues, normconsts_exact, normconsts_and_residues+K,
        residues_exact, errs);
    CHECK_RETCODE(ret_code, release_mem);
//    return E_INVALID_ARGUMENT(kappa);       // here to stop executing

#ifdef DEBUG
    for (UINT i=0; i<5; i++)
        printf("manakovv_testcases_test_fnft: error_bounds[%i] = %2.1e <= %2.1e\n",
            (int)i, errs[i], error_bounds[i]);
    misc_print_buf(3*M, contspec+2*M, "ab_num");
    misc_print_buf(3*M, ab_exact, "ab_exact");

printf("errors\n");
for (UINT i = 0; i<5; i++){
    printf("%f \n", errs[i]);
}
#endif


    // Check if the errors are below the specified bounds. Organized such that
    // the line number tells us which error was too high. The conditions are
    // written in this way to ensure that they fail if an error is NAN.
    if (!(errs[0] <= error_bounds[0])) {    // reflection coef 1
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[1] <= error_bounds[1])) {    // reflection coef 2
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[2] <= error_bounds[2])) {    // a coef
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[3] <= error_bounds[3])) {    // b1
        ret_code = E_TEST_FAILED;
        goto release_mem;
    }
    if (!(errs[4] <= error_bounds[4])) {    // b2
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
