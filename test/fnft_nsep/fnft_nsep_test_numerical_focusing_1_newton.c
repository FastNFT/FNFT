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
 * Sander Wahls (TU Delft) 2020.
 * Shrinivas Chimmalgi (TU Delft) 2020.
 */

#define FNFT_ENABLE_SHORT_NAMES
#include "fnft.h"
#include "fnft_nsep.h"
#include "fnft_nse_discretization_t.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"
#ifdef DEBUG
#include <stdio.h>
#endif

// This test signal was constructed numerically using a squared
// eigenfunctions approach (currently undocumented).

INT fnft_nsep_test_numerical_focusing_newton(fnft_nsep_opts_t opts, REAL error_bounds[2])
{
    INT ret_code = SUCCESS;
    COMPLEX q[257] = {
      4.915570541397106e+00 + 0.000000000000000e+00*I,
      4.842611768461482e+00 - 5.669162201596430e-03*I,
      4.770298742512298e+00 - 1.128814751205249e-02*I,
      4.698669583139001e+00 - 1.685399390052625e-02*I,
      4.627760071499383e+00 - 2.236392104082632e-02*I,
      4.557603719555073e+00 - 2.781532493156223e-02*I,
      4.488231841785164e+00 - 3.320577232373459e-02*I,
      4.419673628866724e+00 - 3.853299499549146e-02*I,
      4.351956222844923e+00 - 4.379488391113549e-02*I,
      4.285104793356807e+00 - 4.898948329825951e-02*I,
      4.219142614495446e+00 - 5.411498467512241e-02*I,
      4.154091141955718e+00 - 5.916972085614094e-02*I,
      4.089970090122955e+00 - 6.415215996182128e-02*I,
      4.026797508798991e+00 - 6.906089945686594e-02*I,
      3.964589859296719e+00 - 7.389466023734968e-02*I,
      3.903362089661035e+00 - 7.865228078577839e-02*I,
      3.843127708800615e+00 - 8.333271141078030e-02*I,
      3.783898859344990e+00 - 8.793500858584583e-02*I,
      3.725686389062726e+00 - 9.245832939987546e-02*I,
      3.668499920700066e+00 - 9.690192613046247e-02*I,
      3.612347920124617e+00 - 1.012651409488811e-01*I,
      3.557237762672274e+00 - 1.055474007646886e-01*I,
      3.503175797620481e+00 - 1.097482122159185e-01*I,
      3.450167410720109e+00 - 1.138671568101258e-01*I,
      3.398217084745376e+00 - 1.179038862194380e-01*I,
      3.347328458016572e+00 - 1.218581177331265e-01*I,
      3.297504380888518e+00 - 1.257296298682472e-01*I,
      3.248746970178738e+00 - 1.295182581403742e-01*I,
      3.201057661549003e+00 - 1.332238909933624e-01*I,
      3.154437259841594e+00 - 1.368464658880374e-01*I,
      3.108885987388902e+00 - 1.403859655483631e-01*I,
      3.064403530321719e+00 - 1.438424143631193e-01*I,
      3.020989082904536e+00 - 1.472158749408878e-01*I,
      2.978641389934364e+00 - 1.505064448155092e-01*I,
      2.937358787242486e+00 - 1.537142532989499e-01*I,
      2.897139240342526e+00 - 1.568394584782060e-01*I,
      2.857980381271579e+00 - 1.598822443526140e-01*I,
      2.819879543674012e+00 - 1.628428181077110e-01*I,
      2.782833796178502e+00 - 1.657214075217180e-01*I,
      2.746839974121205e+00 - 1.685182585005340e-01*I,
      2.711894709668877e+00 - 1.712336327370598e-01*I,
      2.677994460396121e+00 - 1.738678054906419e-01*I,
      2.645135536371524e+00 - 1.764210634823812e-01*I,
      2.613314125807249e+00 - 1.788937029020661e-01*I,
      2.582526319326461e+00 - 1.812860275225056e-01*I,
      2.552768132902448e+00 - 1.835983469170757e-01*I,
      2.524035529522417e+00 - 1.858309747763635e-01*I,
      2.496324439628115e+00 - 1.879842273198573e-01*I,
      2.469630780384245e+00 - 1.900584217987197e-01*I,
      2.443950473824592e+00 - 1.920538750857680e-01*I,
      2.419279463923933e+00 - 1.939709023489240e-01*I,
      2.395613732643221e+00 - 1.958098158044468e-01*I,
      2.372949314992693e+00 - 1.975709235464717e-01*I,
      2.351282313157554e+00 - 1.992545284493950e-01*I,
      2.330608909728007e+00 - 2.008609271398525e-01*I,
      2.310925380071861e+00 - 2.023904090353209e-01*I,
      2.292228103894117e+00 - 2.038432554458987e-01*I,
      2.274513576012904e+00 - 2.052197387369745e-01*I,
      2.257778416394099e+00 - 2.065201215495025e-01*I,
      2.242019379471315e+00 - 2.077446560758034e-01*I,
      2.227233362788474e+00 - 2.088935833880060e-01*I,
      2.213417414988847e+00 - 2.099671328172702e-01*I,
      2.200568743185927e+00 - 2.109655213810425e-01*I,
      2.188684719734817e+00 - 2.118889532568962e-01*I,
      2.177762888432003e+00 - 2.127376193007834e-01*I,
      2.167800970167249e+00 - 2.135116966078650e-01*I,
      2.158796868044836e+00 - 2.142113481145701e-01*I,
      2.150748671996310e+00 - 2.148367222401679e-01*I,
      2.143654662902116e+00 - 2.153879525665060e-01*I,
      2.137513316234346e+00 - 2.158651575549516e-01*I,
      2.132323305239142e+00 - 2.162684402991099e-01*I,
      2.128083503664997e+00 - 2.165978883128208e-01*I,
      2.124792988056224e+00 - 2.168535733519528e-01*I,
      2.122451039611252e+00 - 2.170355512700030e-01*I,
      2.121057145621991e+00 - 2.171438619062567e-01*I,
      2.120611000490146e+00 - 2.171785290068138e-01*I,
      2.121112506333200e+00 - 2.171395601775020e-01*I,
      2.122561773177864e+00 - 2.170269468688455e-01*I,
      2.124959118741475e+00 - 2.168406643930505e-01*I,
      2.128305067804240e+00 - 2.165806719727854e-01*I,
      2.132600351165705e+00 - 2.162469128222654e-01*I,
      2.137845904184281e+00 - 2.158393142607382e-01*I,
      2.144042864890861e+00 - 2.153577878590608e-01*I,
      2.151192571672322e+00 - 2.148022296197012e-01*I,
      2.159296560517358e+00 - 2.141725201907476e-01*I,
      2.168356561804800e+00 - 2.134685251154666e-01*I,
      2.178374496635032e+00 - 2.126900951173669e-01*I,
      2.189352472679839e+00 - 2.118370664226799e-01*I,
      2.201292779540991e+00 - 2.109092611210144e-01*I,
      2.214197883596753e+00 - 2.099064875658026e-01*I,
      2.228070422322403e+00 - 2.088285408156148e-01*I,
      2.242913198053554e+00 - 2.076752031187750e-01*I,
      2.258729171181832e+00 - 2.064462444420792e-01*I,
      2.275521452749891e+00 - 2.051414230461900e-01*I,
      2.293293296423714e+00 - 2.037604861094187e-01*I,
      2.312048089811572e+00 - 2.023031704022727e-01*I,
      2.331789345102341e+00 - 2.007692030148935e-01*I,
      2.352520688990876e+00 - 1.991583021398898e-01*I,
      2.374245851853927e+00 - 1.974701779134116e-01*I,
      2.396968656146474e+00 - 1.957045333167941e-01*I,
      2.420693003980793e+00 - 1.938610651417124e-01*I,
      2.445422863842587e+00 - 1.919394650223830e-01*I,
      2.471162256413916e+00 - 1.899394205371757e-01*I,
      2.497915239453110e+00 - 1.878606163834985e-01*I,
      2.525685891690335e+00 - 1.857027356291728e-01*I,
      2.554478295694233e+00 - 1.834654610437565e-01*I,
      2.584296519663441e+00 - 1.811484765134093e-01*I,
      2.615144598089018e+00 - 1.787514685434914e-01*I,
      2.647026511244834e+00 - 1.762741278522362e-01*I,
      2.679946163450488e+00 - 1.737161510597992e-01*I,
      2.713907360055871e+00 - 1.710772424766443e-01*I,
      2.748913783094387e+00 - 1.683571159953747e-01*I,
      2.784968965551152e+00 - 1.655554970901898e-01*I,
      2.822076264191516e+00 - 1.626721249282033e-01*I,
      2.860238830894610e+00 - 1.597067545969355e-01*I,
      2.899459582438455e+00 - 1.566591594521163e-01*I,
      2.939741168681145e+00 - 1.535291335901282e-01*I,
      2.981085939086402e+00 - 1.503164944490885e-01*I,
      3.023495907539100e+00 - 1.470210855428135e-01*I,
      3.066972715396189e+00 - 1.436427793318942e-01*I,
      3.111517592736602e+00 - 1.401814802347203e-01*I,
      3.157131317746666e+00 - 1.366371277833816e-01*I,
      3.203814174211385e+00 - 1.330096999267499e-01*I,
      3.251565907059530e+00 - 1.292992164847888e-01*I,
      3.300385675947456e+00 - 1.255057427552579e-01*I,
      3.350272006828809e+00 - 1.216293932769238e-01*I,
      3.401222741488856e+00 - 1.176703357509263e-01*I,
      3.453234985062203e+00 - 1.136287951188429e-01*I,
      3.506305051473824e+00 - 1.095050578021197e-01*I,
      3.560428406821714e+00 - 1.052994761014441e-01*I,
      3.615599610755853e+00 - 1.010124727518136e-01*I,
      3.671812255804280e+00 - 9.664454563711990e-02*I,
      3.729058904703785e+00 - 9.219627265978364e-02*I,
      3.787331025853057e+00 - 8.766831675627872e-02*I,
      3.846618926833729e+00 - 8.306143106278964e-02*I,
      3.906911686179549e+00 - 7.837646421699505e-02*I,
      3.968197083408249e+00 - 7.361436579484484e-02*I,
      4.030461527564761e+00 - 6.877619186301426e-02*I,
      4.093689984263164e+00 - 6.386311064800879e-02*I,
      4.157865901451972e+00 - 5.887640830447290e-02*I,
      4.222971134198401e+00 - 5.381749475972517e-02*I,
      4.288985868520386e+00 - 4.868790963228817e-02*I,
      4.355888544595315e+00 - 4.348932819884900e-02*I,
      4.423655779728128e+00 - 3.822356737991860e-02*I,
      4.492262291178770e+00 - 3.289259173641981e-02*I,
      4.561680819293107e+00 - 2.749851944269435e-02*I,
      4.631882051390888e+00 - 2.204362820068575e-02*I,
      4.702834546630877e+00 - 1.653036107819166e-02*I,
      4.774504662444198e+00 - 1.096133222526082e-02*I,
      4.846856482937568e+00 - 5.339332437522906e-03*I,
      4.919851749836012e+00 + 3.326654778102039e-04*I,
      4.993449796474768e+00 + 6.051501572873321e-03*I,
      5.067607485423983e+00 + 1.181382391236223e-02*I,
      5.142279150410953e+00 + 1.761608406067786e-02*I,
      5.217416543215617e+00 + 2.345453295551546e-02*I,
      5.292968786092312e+00 + 2.932521721085940e-02*I,
      5.368882330629127e+00 + 3.522397592020718e-02*I,
      5.445100923755999e+00 + 4.114643801527263e-02*I,
      5.521565581445249e+00 + 4.708802022242415e-02*I,
      5.598214571256547e+00 + 5.304392570636310e-02*I,
      5.674983404357590e+00 + 5.900914345010334e-02*I,
      5.751804837677618e+00 + 6.497844842230510e-02*I,
      5.828608887332887e+00 + 7.094640262048496e-02*I,
      5.905322853980824e+00 + 7.690735704112391e-02*I,
      5.981871360715303e+00 + 8.285545462427060e-02*I,
      6.058176404717978e+00 + 8.878463426703712e-02*I,
      6.134157422897789e+00 + 9.468863592402618e-02*I,
      6.209731372621654e+00 + 1.005610068804010e-01*I,
      6.284812828005621e+00 + 1.063951092340576e-01*I,
      6.359314092235232e+00 + 1.121841286233270e-01*I,
      6.433145326644383e+00 + 1.179210842568719e-01*I,
      6.506214696806211e+00 + 1.235988402654802e-01*I,
      6.578428535899397e+00 + 1.292101183962201e-01*I,
      6.649691525653995e+00 + 1.347475120725910e-01*I,
      6.719906894834377e+00 + 1.402035018173714e-01*I,
      6.788976635197969e+00 + 1.455704720333982e-01*I,
      6.856801734657984e+00 + 1.508407291211585e-01*I,
      6.923282427334311e+00 + 1.560065209086537e-01*I,
      6.988318459742434e+00 + 1.610600573352496e-01*I,
      7.051809372615120e+00 + 1.659935323502457e-01*I,
      7.113654797118524e+00 + 1.707991469299448e-01*I,
      7.173754764668431e+00 + 1.754691331515027e-01*I,
      7.232010028761925e+00 + 1.799957792004192e-01*I,
      7.288322397428798e+00 + 1.843714552032252e-01*I,
      7.342595074782548e+00 + 1.885886397672365e-01*I,
      7.394733009513303e+00 + 1.926399470597223e-01*I,
      7.444643248640414e+00 + 1.965181542957692e-01*I,
      7.492235294109871e+00 + 2.002162294471999e-01*I,
      7.537421460150543e+00 + 2.037273590104518e-01*I,
      7.580117228843396e+00 + 2.070449756355999e-01*I,
      7.620241601511607e+00 + 2.101627854306478e-01*I,
      7.657717443332311e+00 + 2.130747947391159e-01*I,
      7.692471818592464e+00 + 2.157753361906398e-01*I,
      7.724436313989750e+00 + 2.182590938226288e-01*I,
      7.753547347415643e+00 + 2.205211270738334e-01*I,
      7.779746459638188e+00 + 2.225568934491598e-01*I,
      7.802980586553981e+00 + 2.243622696746394e-01*I,
      7.823202309540536e+00 + 2.259335711507168e-01*I,
      7.840370081894862e+00 + 2.272675695473537e-01*I,
      7.854448429181239e+00 + 2.283615083717810e-01*I,
      7.865408121897904e+00 + 2.292131163853272e-01*I,
      7.873226318686212e+00 + 2.298206187312898e-01*I,
      7.877886678971584e+00 + 2.301827456875493e-01*I,
      7.879379443820632e+00 + 2.302987389494636e-01*I,
      7.877701484435846e+00 + 2.301683553980793e-01*I,
      7.872856317684927e+00 + 2.297918683068153e-01*I,
      7.864854088627705e+00 + 2.291700659837446e-01*I,
      7.853711520146298e+00 + 2.283042478576744e-01*I,
      7.839451830169684e+00 + 2.271962180461916e-01*I,
      7.822104617258528e+00 + 2.258482764651882e-01*I,
      7.801705715594834e+00 + 2.242632075610280e-01*I,
      7.778297020682324e+00 + 2.224442667668374e-01*I,
      7.751926287401861e+00 + 2.203951648106698e-01*I,
      7.722646902152412e+00 + 2.181200500100323e-01*I,
      7.690517631100150e+00 + 2.156234887099260e-01*I,
      7.655602346795380e+00 + 2.129104440399835e-01*I,
      7.617969735372232e+00 + 2.099862531628186e-01*I,
      7.577692986867995e+00 + 2.068566032107153e-01*I,
      7.534849471137697e+00 + 2.035275061030107e-01*I,
      7.489520402007791e+00 + 2.000052724496091e-01*I,
      7.441790492175563e+00 + 1.962964847354119e-01*I,
      7.391747601506901e+00 + 1.924079699917801e-01*I,
      7.339482381284963e+00 + 1.883467721533560e-01*I,
      7.285087916880837e+00 + 1.841201242922688e-01*I,
      7.228659371197661e+00 + 1.797354209124503e-01*I,
      7.170293631209538e+00 + 1.752001904844269e-01*I,
      7.110088959779438e+00 + 1.705220683902949e-01*I,
      7.048144654667493e+00 + 1.657087704274210e-01*I,
      6.984560716588303e+00 + 1.607680670152978e-01*I,
      6.919437528015124e+00 + 1.557077582374598e-01*I,
      6.852875544236857e+00 + 1.505356498354747e-01*I,
      6.784974997847509e+00 + 1.452595302467035e-01*I,
      6.715835617864997e+00 + 1.398871487788209e-01*I,
      6.645556364427713e+00 + 1.344261949947644e-01*I,
      6.574235179785983e+00 + 1.288842793638543e-01*I,
      6.501968756058727e+00 + 1.232689152156579e-01*I,
      6.428852320323296e+00 + 1.175875020406799e-01*I,
      6.354979437174563e+00 + 1.118473101484821e-01*I,
      6.280441828779566e+00 + 1.060554666852808e-01*I,
      6.205329212516692e+00 + 1.002189430179189e-01*I,
      6.129729155918959e+00 + 9.434454346242564e-02*I,
      6.053726948610743e+00 + 8.843889533303768e-02*I,
      5.977405490841972e+00 + 8.250844028090516e-02*I,
      5.900845198197560e+00 + 7.655942688966112e-02*I,
      5.824123921802808e+00 + 7.059790447508869e-02*I,
      5.747316883389937e+00 + 6.462971803956485e-02*I,
      5.670496624616373e+00 + 5.866050433389933e-02*I,
      5.593732969821312e+00 + 5.269568896338130e-02*I,
      5.517093001421384e+00 + 4.674048447593922e-02*I,
      5.440641047207875e+00 + 4.079988937508964e-02*I,
      5.364438678735176e+00 + 3.487868799470606e-02*I,
      5.288544719912481e+00 + 2.898145116663033e-02*I,
      5.213015265039893e+00 + 2.311253762214577e-02*I,
      5.137903705486382e+00 + 1.727609606494724e-02*I,
      5.063260764184878e+00 + 1.147606785153594e-02*I,
      4.989134537175214e+00 + 5.716190219260737e-03*I,
      4.915570541397106e+00 + 0.000000000000000e+00*I
    };
    const REAL T[2] = {0, 6.628591515456010e-01};

    const COMPLEX mainspec_exact[6] = {-5*I, -2*I, -1*I, 1*I, 2*I, 5*I};
    const COMPLEX auxspec_exact[2] = {
        -8.363052072415789e-02 + 2.059857349398656e+00*I,
        -1.400213346927766e-01 + 8.184100833800795e-01*I
    };
    const UINT K_exact = sizeof(mainspec_exact) / sizeof(mainspec_exact[0]);
    const UINT M_exact = sizeof(auxspec_exact) / sizeof(auxspec_exact[0]);

    const UINT D = sizeof(q) / sizeof(q[0]);
    COMPLEX mainspec[12] = {-5*I+0.1, -2*I-0.1*I, -1*I-0.1, 1*I+0.03*I, 2*I+0.001-0.7*I, 5*I-0.4+0.1*I, -5*I+0.3*I-0.34, -2*I+0.4, -1*I+0.1, 1*I+0.1*I, 2*I-0.213, 5*I+0.023*I};
    COMPLEX auxspec[2] = {
        -8e-02 + 2e+00*I,
        -1e-01 + 8e-01*I
    };
    UINT K = sizeof(mainspec) / sizeof(mainspec[0]) / opts.points_per_spine;
    UINT M = sizeof(auxspec) / sizeof(auxspec[0]);

    // Compute main and auxiliary spectrum
    REAL phase_shift = CARG(q[D-1]/q[0]);
    ret_code = fnft_nsep(D-1, q, T, phase_shift, &K, mainspec, &M, auxspec, NULL, +1, &opts);
    CHECK_RETCODE(ret_code, leave_fun);

    // Check that the main and auxiliary spectrum are correct
    REAL dist = misc_hausdorff_dist(K_exact, mainspec_exact, K, mainspec);
#ifdef DEBUG
    printf("dist(mainspec) = %g\n", dist);
#endif
    if (dist > error_bounds[0]) {
#ifdef DEBUG
        misc_print_buf(K, mainspec, "mainspec");
#endif
        ret_code = E_TEST_FAILED;
        goto leave_fun;
    }
    dist = misc_hausdorff_dist(M_exact, auxspec_exact, M, auxspec);
#ifdef DEBUG
    printf("dist(auxspec) = %g\n", dist);
#endif
   if (dist > error_bounds[1]) {
#ifdef DEBUG
        misc_print_buf(M, auxspec, "auxspec");
#endif
       ret_code = E_TEST_FAILED;
        goto leave_fun;
    }

leave_fun:
    return ret_code;
}

int main()
{
    fnft_nsep_opts_t opts = fnft_nsep_default_opts();
    opts.localization = fnft_nsep_loc_NEWTON;
    opts.filtering = fnft_nsep_filt_MANUAL;
    opts.bounding_box[0] = -1;
    opts.bounding_box[1] = 1;
    opts.bounding_box[2] = -10;
    opts.bounding_box[3] = 10;

    REAL error_bounds[2] = { 1.5e-4 /* mainspec */, 1.3e-2 /* auxspec */ };
    opts.discretization = fnft_nse_discretization_BO;
    if (fnft_nsep_test_numerical_focusing_newton(opts, error_bounds) != SUCCESS)
        return EXIT_FAILURE;
 
    error_bounds[0] = 1e-8;
    opts.discretization = fnft_nse_discretization_CF4_2;
    if (fnft_nsep_test_numerical_focusing_newton(opts, error_bounds) != SUCCESS)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
