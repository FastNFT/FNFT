UINT D = 512;
REAL X[2] = {0, 9.597341792231519619349455751944333314895629882812500000000000e+00};
REAL E[2] = {-3.100200000000000066791017161449417471885681152343750000000000e+00, 2.000099999999999766941982670687139034271240234375000000000000e+00};
REAL main_spec_exact[3] = {-1.400099999999999900168745625705923885107040405273437500000000e+00, -1.399999999999999911182158029987476766109466552734375000000000e+00, 2.999999999999999888977697537484345957636833190917968750000000e-01};
REAL aux_spec_exact[1] = {-1.298686971821685354200326401041820645332336425781250000000000e+00};
COMPLEX q[512] = {
	-1.020239673335621155558783357264474034309387207031250000000000e-01,
	-1.112025061469483766529720014659687876701354980468750000000000e-01,
	-1.199674560936290212964649981586262583732604980468750000000000e-01,
	-1.283363532671883966429504653206095099449157714843750000000000e-01,
	-1.363260926517630977627959509845823049545288085937500000000000e-01,
	-1.439529416833704722478159965248778462409973144531250000000000e-01,
	-1.512325545947650340394829981960356235504150390625000000000000e-01,
	-1.581799873932983047097877715714275836944580078125000000000000e-01,
	-1.648097133364148092482537322212010622024536132812500000000000e-01,
	-1.711356387834279679793780815089121460914611816406250000000000e-01,
	-1.771711193151590535421746608335524797439575195312500000000000e-01,
	-1.829289760248729024283420585561543703079223632812500000000000e-01,
	-1.884215118948984191149520484032109379768371582031250000000000e-01,
	-1.936605281832664005037258903030306100845336914062500000000000e-01,
	-1.986573407538752000789372686995193362236022949218750000000000e-01,
	-2.034227962919756826920547609915956854820251464843750000000000e-01,
	-2.079672883543888239898933534277603030204772949218750000000000e-01,
	-2.123007732107078382455256360117346048355102539062500000000000e-01,
	-2.164327854380059257266566419275477528572082519531250000000000e-01,
	-2.203724532371738042968445370206609368324279785156250000000000e-01,
	-2.241285134441204363042743352707475423812866210937500000000000e-01,
	-2.277093262136262730166436085710301995277404785156250000000000e-01,
	-2.311228893577523813007701392052695155143737792968750000000000e-01,
	-2.343768523243716650483747798716649413108825683593750000000000e-01,
	-2.374785298046862003218393510906025767326354980468750000000000e-01,
	-2.404349149614981584477391152177006006240844726562500000000000e-01,
	-2.432526922725988249851525324629619717597961425781250000000000e-01,
	-2.459382499859463777625023794826120138168334960937500000000000e-01,
	-2.484976921852792841605150897521525621414184570312500000000000e-01,
	-2.509368504666418253634674329077824950218200683593750000000000e-01,
	-2.532612952277855100646775099448859691619873046875000000000000e-01,
	-2.554763465737970307145587867125868797302246093750000000000000e-01,
	-2.575870848434300697249454970005899667739868164062500000000000e-01,
	-2.595983607616321187094854394672438502311706542968750000000000e-01,
	-2.615148052245945820004635606892406940460205078125000000000000e-01,
	-2.633408387243512116526744648581370711326599121093750000000000e-01,
	-2.650806804205632083437649271218106150627136230468750000000000e-01,
	-2.667383568675965044292297534411773085594177246093750000000000e-01,
	-2.683177104053751094170365831814706325531005859375000000000000e-01,
	-2.698224072228243564097738271811977028846740722656250000000000e-01,
	-2.712559451029101786900810111546888947486877441406250000000000e-01,
	-2.726216608584621781119494698941707611083984375000000000000000e-01,
	-2.739227374680592852485006005736067891120910644531250000000000e-01,
	-2.751622109213034406138831400312483310699462890625000000000000e-01,
	-2.763429767828151639719180820975452661514282226562500000000000e-01,
	-2.774677964842395816447151446482166647911071777343750000000000e-01,
	-2.785393033534893092451056872960180044174194335937500000000000e-01,
	-2.795600083903373445082252146676182746887207031250000000000000e-01,
	-2.805323057973647671303751849336549639701843261718750000000000e-01,
	-2.814584782750984004451311193406581878662109375000000000000000e-01,
	-2.823407020900243757921543874545022845268249511718750000000000e-01,
	-2.831810519239810197689166670897975564002990722656250000000000e-01,
	-2.839815055132284271621756488457322120666503906250000000000000e-01,
	-2.847439480853029003526444284943863749504089355468750000000000e-01,
	-2.854701766015526054331985505996271967887878417968750000000000e-01,
	-2.861619038130229775163115846225991845130920410156250000000000e-01,
	-2.868207621371525739562002854654565453529357910156250000000000e-01,
	-2.874483073625100360004580579698085784912109375000000000000000e-01,
	-2.880460221885789984241910133277997374534606933593750000000000e-01,
	-2.886153196073752980055360239930450916290283203125000000000000e-01,
	-2.891575461334534580259969516191631555557250976562500000000000e-01,
	-2.896739848886400459093692916212603449821472167968750000000000e-01,
	-2.901658585476090124188885965850204229354858398437500000000000e-01,
	-2.906343321501982934762509103165939450263977050781250000000000e-01,
	-2.910805157861560132914746645838022232055664062500000000000000e-01,
	-2.915054671577870237797469599172472953796386718750000000000000e-01,
	-2.919101940257751159890631242888048291206359863281250000000000e-01,
	-2.922956565432519582259374146815389394760131835937500000000000e-01,
	-2.926627694829835313328203483251854777336120605468750000000000e-01,
	-2.930124043623618668163999245734885334968566894531250000000000e-01,
	-2.933453914706980469873087713494896888732910156250000000000000e-01,
	-2.936625218031401196583374257897958159446716308593750000000000e-01,
	-2.939645489053552829261661827331408858299255371093750000000000e-01,
	-2.942521906329558234460819221567362546920776367187500000000000e-01,
	-2.945261308294777613525639026192948222160339355468750000000000e-01,
	-2.947870209265648355767552857287228107452392578125000000000000e-01,
	-2.950354814698634697833767859265208244323730468750000000000000e-01,
	-2.952721035739753752125125174643471837043762207031250000000000e-01,
	-2.954974503096812199487430916633456945419311523437500000000000e-01,
	-2.957120580265102383066277980105951428413391113281250000000000e-01,
	-2.959164376135934304556940332986414432525634765625000000000000e-01,
	-2.961110757016163219645932258572429418563842773437500000000000e-01,
	-2.962964358085655725005835847696289420127868652343750000000000e-01,
	-2.964729594318384897633222863078117370605468750000000000000000e-01,
	-2.966410670891792555892152449814602732658386230468750000000000e-01,
	-2.968011593107990897522086015669628977775573730468750000000000e-01,
	-2.969536175849176729002465435769408941268920898437500000000000e-01,
	-2.970988052588885430793652631109580397605895996093750000000000e-01,
	-2.972370683979486116754742397461086511611938476562500000000000e-01,
	-2.973687366035617785087197262328118085861206054687500000000000e-01,
	-2.974941237932235971186401002341881394386291503906250000000000e-01,
	-2.976135289435188902018580847652629017829895019531250000000000e-01,
	-2.977272367981389500357636279659345746040344238281250000000000e-01,
	-2.978355185424867990207076218212023377418518066406250000000000e-01,
	-2.979386324464297075564900296740233898162841796875000000000000e-01,
	-2.980368244766813390356219315435737371444702148437500000000000e-01,
	-2.981303288802297224435733369318768382072448730468750000000000e-01,
	-2.982193687401721859941972070373594760894775390625000000000000e-01,
	-2.983041565052322319218092161463573575019836425781250000000000e-01,
	-2.983848944942057990203920780913904309272766113281250000000000e-01,
	-2.984617753764982062136823515174910426139831542968750000000000e-01,
	-2.985349826298802078383687330642715096473693847656250000000000e-01,
	-2.986046909765267542979927384294569492340087890625000000000000e-01,
	-2.986710667983589750917872152058407664299011230468750000000000e-01,
	-2.987342685326619395880243246210739016532897949218750000000000e-01,
	-2.987944470489001247415217221714556217193603515625000000000000e-01,
	-2.988517460076201004426366125699132680892944335937500000000000e-01,
	-2.989063022022775406583150470396503806114196777343750000000000e-01,
	-2.989582458847901413889758259756490588188171386718750000000000e-01,
	-2.990077010755882724879484158009290695190429687500000000000000e-01,
	-2.990547858588867846663106320193037390708923339843750000000000e-01,
	-2.990996126638698626720724860206246376037597656250000000000000e-01,
	-2.991422885324634961534684407524764537811279296875000000000000e-01,
	-2.991829153743119640296299621695652604103088378906250000000000e-01,
	-2.992215902095707313890216028084978461265563964843750000000000e-01,
	-2.992584054000837490150388475740328431129455566406250000000000e-01,
	-2.992934488694904970884635986294597387313842773437500000000000e-01,
	-2.993268043127859101559806731529533863067626953125000000000000e-01,
	-2.993585513958274546553184336517006158828735351562500000000000e-01,
	-2.993887659452587612918250670190900564193725585937500000000000e-01,
	-2.994175201293036714389472763286903500556945800781250000000000e-01,
	-2.994448826298556909364378952886909246444702148437500000000000e-01,
	-2.994709188062691929133052326506003737449645996093750000000000e-01,
	-2.994956908512458326754313020501285791397094726562500000000000e-01,
	-2.995192579391807718991458386881276965141296386718750000000000e-01,
	-2.995416763673208748741672025062143802642822265625000000000000e-01,
	-2.995629996900737168630257656332105398178100585937500000000000e-01,
	-2.995832788467813756483337783720344305038452148437500000000000e-01,
	-2.996025622832649837334884068695828318595886230468750000000000e-01,
	-2.996208960674278110047907830448821187019348144531250000000000e-01,
	-2.996383239991935454327176557853817939758300781250000000000000e-01,
	-2.996548877150360112864291295409202575683593750000000000000000e-01,
	-2.996706267873547879787565761944279074668884277343750000000000e-01,
	-2.996855788189285441092124528950080275535583496093750000000000e-01,
	-2.996997795326739044696751079754903912544250488281250000000000e-01,
	-2.997132628569181278521682543214410543441772460937500000000000e-01,
	-2.997260610063978703010434401221573352813720703125000000000000e-01,
	-2.997382045591652222071843425510451197624206542968750000000000e-01,
	-2.997497225295928657828881114255636930465698242187500000000000e-01,
	-2.997606424376479949955864867661148309707641601562500000000000e-01,
	-2.997709903746028636817300139227882027626037597656250000000000e-01,
	-2.997807910653329521721843775594606995582580566406250000000000e-01,
	-2.997900679273572954741666762856766581535339355468750000000000e-01,
	-2.997988431267599729324047075351700186729431152343750000000000e-01,
	-2.998071376311211011511659307871013879776000976562500000000000e-01,
	-2.998149712595918892077406781027093529701232910156250000000000e-01,
	-2.998223627302243343706322775688022375106811523437500000000000e-01,
	-2.998293297046772387659530068049207329750061035156250000000000e-01,
	-2.998358888303949143505633401218801736831665039062500000000000e-01,
	-2.998420557803695984944170049857348203659057617187500000000000e-01,
	-2.998478452905746216572424600599333643913269042968750000000000e-01,
	-2.998532711951611418044194579124450683593750000000000000000000e-01,
	-2.998583464595019343335025041596964001655578613281250000000000e-01,
	-2.998630832111595090339051239425316452980041503906250000000000e-01,
	-2.998674927688549374238391465041786432266235351562500000000000e-01,
	-2.998715856695040038459865172626450657844543457031250000000000e-01,
	-2.998753716933886259710106969578191637992858886718750000000000e-01,
	-2.998788598875159472356699552619829773902893066406250000000000e-01,
	-2.998820585872334909538494684966281056404113769531250000000000e-01,
	-2.998849754361381236833494767779484391212463378906250000000000e-01,
	-2.998876174043365594457100087311118841171264648437500000000000e-01,
	-2.998899908050951523819094290956854820251464843750000000000000e-01,
	-2.998921013099229426757119654212146997451782226562500000000000e-01,
	-2.998939539621190419893537182360887527465820312500000000000000e-01,
	-2.998955531888207737267748598242178559303283691406250000000000e-01,
	-2.998969028115809898338284256169572472572326660156250000000000e-01,
	-2.998980060554958804175385012058541178703308105468750000000000e-01,
	-2.998988655569121419830480590462684631347656250000000000000000e-01,
	-2.998994833697219419832435960415750741958618164062500000000000e-01,
	-2.998998609702723250336475757649168372154235839843750000000000e-01,
	-2.998999992608912812386279256315901875495910644531250000000000e-01,
	-2.998998985720433552160102408379316329956054687500000000000000e-01,
	-2.998995586631179044445616455050185322761535644531250000000000e-01,
	-2.998989787218566682724940619664266705513000488281250000000000e-01,
	-2.998981573624126539812095870729535818099975585937500000000000e-01,
	-2.998970926220372312798190250759944319725036621093750000000000e-01,
	-2.998957819563936588735941768391057848930358886718750000000000e-01,
	-2.998942222334775031811204826226457953453063964843750000000000e-01,
	-2.998924097261346233267431671265512704849243164062500000000000e-01,
	-2.998903401031576265722833340987563133239746093750000000000000e-01,
	-2.998880084189377015491118072532117366790771484375000000000000e-01,
	-2.998854091016509570977177645545452833175659179687500000000000e-01,
	-2.998825359399468482024531112983822822570800781250000000000000e-01,
	-2.998793820681115995796517381677404046058654785156250000000000e-01,
	-2.998759399496653266226076084421947598457336425781250000000000e-01,
	-2.998722013593586588342532195383682847023010253906250000000000e-01,
	-2.998681573635253450049731327453628182411193847656250000000000e-01,
	-2.998637982987402139656296640168875455856323242187500000000000e-01,
	-2.998591137487371938163960294332355260848999023437500000000000e-01,
	-2.998540925195265494096474867546930909156799316406250000000000e-01,
	-2.998487226126558269356792152393609285354614257812500000000000e-01,
	-2.998429911965518890326620748965069651603698730468750000000000e-01,
	-2.998368845758694334335814346559345722198486328125000000000000e-01,
	-2.998303881587798258578914101235568523406982421875000000000000e-01,
	-2.998234864221189788224819494644179940223693847656250000000000e-01,
	-2.998161628743103435112971055787056684494018554687500000000000e-01,
	-2.998084000159768613968935824232175946235656738281250000000000e-01,
	-2.998001792981477287014513422036543488502502441406250000000000e-01,
	-2.997914810779578331789707590360194444656372070312500000000000e-01,
	-2.997822845717363904327612544875591993331909179687500000000000e-01,
	-2.997725678053728692873391992179676890373229980468750000000000e-01,
	-2.997623075618438548417543643154203891754150390625000000000000e-01,
	-2.997514793257689547090194537304341793060302734375000000000000e-01,
	-2.997400572248714034628847002750262618064880371093750000000000e-01,
	-2.997280139682025890124350553378462791442871093750000000000000e-01,
	-2.997153207809768460379018506500869989395141601562500000000000e-01,
	-2.997019473358704111376482615014538168907165527343750000000000e-01,
	-2.996878616806131212513264472363516688346862792968750000000000e-01,
	-2.996730301617085423515618458623066544532775878906250000000000e-01,
	-2.996574173440929023115586460335180163383483886718750000000000e-01,
	-2.996409859265449782128598599229007959365844726562500000000000e-01,
	-2.996236966526497624840885691810399293899536132812500000000000e-01,
	-2.996055082170960837117945629870519042015075683593750000000000e-01,
	-2.995863771670923547674192377598956227302551269531250000000000e-01,
	-2.995662577986624164338991249678656458854675292968750000000000e-01,
	-2.995451020475763392880708124721422791481018066406250000000000e-01,
	-2.995228593746541712050657224608585238456726074218750000000000e-01,
	-2.994994766451788414940438087796792387962341308593750000000000e-01,
	-2.994748980021171291809878312051296234130859375000000000000000e-01,
	-2.994490647328631460766246163984760642051696777343750000000000e-01,
	-2.994219151291779290602335095172747969627380371093750000000000e-01,
	-2.993933843399960714748431200860068202018737792968750000000000e-01,
	-2.993634042167516717825037630973383784294128417968750000000000e-01,
	-2.993319031508541172570403432473540306091308593750000000000000e-01,
	-2.992988059029308978153949283296242356300354003906250000000000e-01,
	-2.992640334234315524497560545569285750389099121093750000000000e-01,
	-2.992275026641677548866482538869604468345642089843750000000000e-01,
	-2.991891263803467815307612909236922860145568847656250000000000e-01,
	-2.991488129226231862389795423950999975204467773437500000000000e-01,
	-2.991064660186868451319242012687027454376220703125000000000000e-01,
	-2.990619845438602375509162811795249581336975097656250000000000e-01,
	-2.990152622801725001977501960936933755874633789062500000000000e-01,
	-2.989661876633356030197319341823458671569824218750000000000000e-01,
	-2.989146435170302318340418423758819699287414550781250000000000e-01,
	-2.988605067738703269242250826209783554077148437500000000000000e-01,
	-2.988036481823939105595400178572162985801696777343750000000000e-01,
	-2.987439319993843156453294795937836170196533203125000000000000e-01,
	-2.986812156668072759657661663368344306945800781250000000000000e-01,
	-2.986153494725938273290921642910689115524291992187500000000000e-01,
	-2.985461761944798730894490290665999054908752441406250000000000e-01,
	-2.984735307260670822415704606100916862487792968750000000000000e-01,
	-2.983972396842209384715260966913774609565734863281250000000000e-01,
	-2.983171209968973336401631968328729271888732910156250000000000e-01,
	-2.982329834704238180620450293645262718200683593750000000000000e-01,
	-2.981446263352309777872051199665293097496032714843750000000000e-01,
	-2.980518387689703452281264617340639233589172363281250000000000e-01,
	-2.979543993959090641965303802862763404846191406250000000000000e-01,
	-2.978520757614320224604398390511050820350646972656250000000000e-01,
	-2.977446237804368678325772634707391262054443359375000000000000e-01,
	-2.976317871583300522786430519772693514823913574218750000000000e-01,
	-2.975132967832943009511836862657219171524047851562500000000000e-01,
	-2.973888700884090852127883408684283494949340820312500000000000e-01,
	-2.972582103821590493453186354599893093109130859375000000000000e-01,
	-2.971210061457774109783258609240874648094177246093750000000000e-01,
	-2.969769302958131795833196520106866955757141113281250000000000e-01,
	-2.968256394102217754493722168263047933578491210937500000000000e-01,
	-2.966667729162084654603859235066920518875122070312500000000000e-01,
	-2.964999522379616614387032313970848917961120605468750000000000e-01,
	-2.963247799023265294238171918550506234169006347656250000000000e-01,
	-2.961408386003925308216366829583421349525451660156250000000000e-01,
	-2.959476902028455036486320750555023550987243652343750000000000e-01,
	-2.957448747268656141784504143288359045982360839843750000000000e-01,
	-2.955319092522263879629917937563732266426086425781250000000000e-01,
	-2.953082867841554381982405175222083926200866699218750000000000e-01,
	-2.950734750603980494076949980808421969413757324218750000000000e-01,
	-2.948269152998137521137778094271197915077209472656250000000000e-01,
	-2.945680208897059060291212517768144607543945312500000000000000e-01,
	-2.942961760089701783726923167705535888671875000000000000000000e-01,
	-2.940107341840052512793590722139924764633178710937500000000000e-01,
	-2.937110167741976418653848668327555060386657714843750000000000e-01,
	-2.933963113836512981436044356087222695350646972656250000000000e-01,
	-2.930658701956856404535756155382841825485229492187500000000000e-01,
	-2.927189082264707309377627098001539707183837890625000000000000e-01,
	-2.923546014940154869066191167803481221199035644531250000000000e-01,
	-2.919720850985600968385824671713635325431823730468750000000000e-01,
	-2.915704512102528234152032382553443312644958496093750000000000e-01,
	-2.911487469598190713782059901859611272811889648437500000000000e-01,
	-2.907059722277596236494900949764996767044067382812500000000000e-01,
	-2.902410773274133326538049004739150404930114746093750000000000e-01,
	-2.897529605770454708135730470530688762664794921875000000000000e-01,
	-2.892404657559239922193228267133235931396484375000000000000000e-01,
	-2.887023794391381237289806449552997946739196777343750000000000e-01,
	-2.881374282057227453890391188906505703926086425781250000000000e-01,
	-2.875442757144308636441110138548538088798522949218750000000000e-01,
	-2.869215196412904234080087917391210794448852539062500000000000e-01,
	-2.862676884728680981595516641391441226005554199218750000000000e-01,
	-2.855812381489291063019209104822948575019836425781250000000000e-01,
	-2.848605485479800414339024428045377135276794433593750000000000e-01,
	-2.841039198089361228483085142215713858604431152343750000000000e-01,
	-2.833095684819544324284379399614408612251281738281250000000000e-01,
	-2.824756235012229055314492143224924802780151367187500000000000e-01,
	-2.816001219722910065001997281797230243682861328125000000000000e-01,
	-2.806810047663028662157103099161759018898010253906250000000000e-01,
	-2.797161119132685058730203309096395969390869140625000000000000e-01,
	-2.787031777863040460374577378388494253158569335937500000000000e-01,
	-2.776398260685692953586567455204203724861145019531250000000000e-01,
	-2.765235644944290527291741454973816871643066406250000000000000e-01,
	-2.753517793561837123661462101154029369354248046875000000000000e-01,
	-2.741217297675557773573018494062125682830810546875000000000000e-01,
	-2.728305416749501333129046543035656213760375976562500000000000e-01,
	-2.714752016073931351058945438126102089881896972656250000000000e-01,
	-2.700525501559383201310993172228336334228515625000000000000000e-01,
	-2.685592751732608363113286031875759363174438476562500000000000e-01,
	-2.669919046841049414808821893529966473579406738281250000000000e-01,
	-2.653467994972644738993494684109464287757873535156250000000000e-01,
	-2.636201455097961776630199892679229378700256347656250000000000e-01,
	-2.618079456942714600131694169249385595321655273437500000000000e-01,
	-2.599060117600164865336864750133827328681945800781250000000000e-01,
	-2.579099554794961335346670239232480525970458984375000000000000e-01,
	-2.558151796713103998115457216044887900352478027343750000000000e-01,
	-2.536168688316107200364513118984177708625793457031250000000000e-01,
	-2.513099794062418901319233555113896727561950683593750000000000e-01,
	-2.488892296964562156347255950095131993293762207031250000000000e-01,
	-2.463490893917796853429535985924303531646728515625000000000000e-01,
	-2.436837687243738059805764351040124893188476562500000000000000e-01,
	-2.408872072402554742609481763793155550956726074218750000000000e-01,
	-2.379530621838470416662403295049443840980529785156250000000000e-01,
	-2.348746964936356818043350358493626117706298828125000000000000e-01,
	-2.316451664082448402837144385557621717453002929687500000000000e-01,
	-2.282572086839440572703097132034599781036376953125000000000000e-01,
	-2.247032274266187457101295876782387495040893554687500000000000e-01,
	-2.209752805434308520204922388074919581413269042968750000000000e-01,
	-2.170650658219983597518876194953918457031250000000000000000000e-01,
	-2.129639066477291287071693659527227282524108886718750000000000e-01,
	-2.086627373732108381432226451579481363296508789062500000000000e-01,
	-2.041520883571430466929541580611839890480041503906250000000000e-01,
	-1.994220706943528043098012858536094427108764648437500000000000e-01,
	-1.944623606628930190254322951659560203552246093750000000000000e-01,
	-1.892621839192965005338464834494516253471374511718750000000000e-01,
	-1.838102994784969190789070125902071595191955566406250000000000e-01,
	-1.780949835211740328588803095044568181037902832031250000000000e-01,
	-1.721040130779774024460948567138984799385070800781250000000000e-01,
	-1.658246496475865860986687039257958531379699707031250000000000e-01,
	-1.592436228137343867672370834043249487876892089843750000000000e-01,
	-1.523471139353924641568482911679893732070922851562500000000000e-01,
	-1.451207399940765974122314219130203127861022949218750000000000e-01,
	-1.375495376930500057710560213308781385421752929687500000000000e-01,
	-1.296179479149475177734984754351899027824401855468750000000000e-01,
	-1.213098006570172415763408935163170099258422851562500000000000e-01,
	-1.126083005771167933062315569259226322174072265625000000000000e-01,
	-1.034960132984981129311563563533127307891845703125000000000000e-01,
	-9.395485263761838012896987493149936199188232421875000000000000e-02,
	-8.396606893663793869109213119372725486755371093750000000000000e-02,
	-7.351023870095962919890553166624158620834350585937500000000000e-02,
	-6.256725576219457352067365718539804220199584960937500000000000e-02,
	-5.111632420834633094841592537704855203628540039062500000000000e-02,
	-3.913595334577468776160458219237625598907470703125000000000000e-02,
	-2.660395498159351390654592250939458608627319335937500000000000e-02,
	-1.349744334077707108576760219875723123550415039062500000000000e-02,
	2.071620410424301894636300858110189437866210937500000000000000e-04,
	1.453413007985560057022667024284601211547851562500000000000000e-02,
	2.950840816220035023320633627008646726608276367187500000000000e-02,
	4.515560799189843876888517115730792284011840820312500000000000e-02,
	6.150198680445795318405544094275683164596557617187500000000000e-02,
	7.857442345709753084292970015667378902435302734375000000000000e-02,
	9.640038886993096411970327608287334442138671875000000000000000e-02,
	1.150079102612719950116115796845406293869018554687500000000000e-01,
	1.344255285877982153408538579242303967475891113281250000000000e-01,
	1.546822485687728843117838550824671983718872070312500000000000e-01,
	1.758074806431841885512312728678807616233825683593750000000000e-01,
	1.978309741802588561654374643694609403610229492187500000000000e-01,
	2.207827412379996800950721080880612134933471679687500000000000e-01,
	2.446929701420947900203373137628659605979919433593750000000000e-01,
	2.695919281393388544643130444455891847610473632812500000000000e-01,
	2.955098523675250277165105217136442661285400390625000000000000e-01,
	3.224768283775039101612946979003027081489562988281250000000000e-01,
	3.505226554454541343375240103341639041900634765625000000000000e-01,
	3.796766979247676587760906841140240430831909179687500000000000e-01,
	4.099677219094584046388263232074677944183349609375000000000000e-01,
	4.414237165163585885352404147852212190628051757812500000000000e-01,
	4.740716991429096971444323571631684899330139160156250000000000e-01,
	5.079375041238203269955420182668603956699371337890625000000000e-01,
	5.430455542945262781273640939616598188877105712890625000000000e-01,
	5.794186150747857588783062965376302599906921386718750000000000e-01,
	6.170775308136313697104924358427524566650390625000000000000000e-01,
	6.560409432899105830472308298340067267417907714843750000000000e-01,
	6.963249924420034098204723704839125275611877441406250000000000e-01,
	7.379429996085067333666529520996846258640289306640625000000000e-01,
	7.809051338002497288925951579585671424865722656250000000000000e-01,
	8.252180617938886797446684795431792736053466796875000000000000e-01,
	8.708845831399727455135462150792591273784637451171875000000000e-01,
	9.179032515141345438536291112541221082210540771484375000000000e-01,
	9.662679842080730274744837515754625201225280761718750000000000e-01,
	1.015967661957445811182765282865148037672042846679687500000000e+00,
	1.066985721733452407988806953653693199157714843750000000000000e+00,
	1.119299745582001737176369715598411858081817626953125000000000e+00,
	1.172881049074021353462171646242495626211166381835937500000000e+00,
	1.227694273427362592343570213415659964084625244140625000000000e+00,
	1.283696985868436613031917659100145101547241210937500000000000e+00,
	1.340839293310869040354305070650298148393630981445312500000000e+00,
	1.399063474930185790867653850000351667404174804687500000000000e+00,
	1.458303639694966680195875596837140619754791259765625000000000e+00,
	1.518485415363354240625426427868660539388656616210937500000000e+00,
	1.579525675854905886197343534149695187807083129882812500000000e+00,
	1.641332314243649559415416661067865788936614990234375000000000e+00,
	1.703804068873963295516205107560381293296813964843750000000000e+00,
	1.766830410255569194788449749466963112354278564453125000000000e+00,
	1.830291496431678766043660289142280817031860351562500000000000e+00,
	1.894058204414277746252537326654419302940368652343750000000000e+00,
	1.957992245027075028929175459779798984527587890625000000000000e+00,
	2.021946368071851019010409800102934241294860839843750000000000e+00,
	2.085764664125197231925312735256738960742950439453125000000000e+00,
	2.149282968467778776044951882795430719852447509765625000000000e+00,
	2.212329371641800879189077022601850330829620361328125000000000e+00,
	2.274724839920290619943443743977695703506469726562500000000000e+00,
	2.336283947559594231790924823144450783729553222656250000000000e+00,
	2.396815721102572904754879346000961959362030029296875000000000e+00,
	2.456124594221426704621080716606229543685913085937500000000000e+00,
	2.514011469660566877593055323814041912555694580078125000000000e+00,
	2.570274882792689119526130525628104805946350097656250000000000e+00,
	2.624712259173392236988320291857235133647918701171875000000000e+00,
	2.677121256318836461929322467767633497714996337890625000000000e+00,
	2.727301177786958241711090522585436701774597167968750000000000e+00,
	2.775054445573628569832180801313370466232299804687500000000000e+00,
	2.820188114900094689119214308448135852813720703125000000000000e+00,
	2.862515413728829116735141724348068237304687500000000000000000e+00,
	2.901857287862253276244928201776929199695587158203125000000000e+00,
	2.938043931309806566787301562726497650146484375000000000000000e+00,
	2.970916280807200937630341286421753466129302978515625000000000e+00,
	3.000327452979780318287339468952268362045288085937500000000000e+00,
	3.026144102696783910744215972954407334327697753906250000000000e+00,
	3.048247681684608245689105388009920716285705566406250000000000e+00,
	3.066535577466045214833911813911981880664825439453125000000000e+00,
	3.080922114160854263786859519314020872116088867187500000000000e+00,
	3.091339398601030019619884114945307374000549316406250000000000e+00,
	3.097737997544327104293415686697699129581451416015625000000000e+00,
	3.100087434461361812765289869275875389575958251953125000000000e+00,
	3.098376497361218895321144373156130313873291015625000000000000e+00,
	3.092613352333047682662936495034955441951751708984375000000000e+00,
	3.082825460833978237218389040208421647548675537109375000000000e+00,
	3.069059302158565927243216719944030046463012695312500000000000e+00,
	3.051379905891901334058502470725215971469879150390625000000000e+00,
	3.029870202388833977380500073195435106754302978515625000000000e+00,
	3.004630202351624035372879006899893283843994140625000000000000e+00,
	2.975776019321503262915484810946509242057800292968750000000000e+00,
	2.943438751290420984929596670554019510746002197265625000000000e+00,
	2.907763239624849838094178267056122422218322753906250000000000e+00,
	2.868906725034569671350936914677731692790985107421875000000000e+00,
	2.827037421392662785990523843793198466300964355468750000000000e+00,
	2.782333028809559127836337211192585527896881103515625000000000e+00,
	2.734979207490494168553141207667067646980285644531250000000000e+00,
	2.685168033582744584464307990856468677520751953125000000000000e+00,
	2.633096457478643959149167130817659199237823486328125000000000e+00,
	2.578964783926000148284174429136328399181365966796875000000000e+00,
	2.522975191860799881737875693943351507186889648437500000000000e+00,
	2.465330310173668326712004272849299013614654541015625000000000e+00,
	2.406231863713893304179691767785698175430297851562500000000000e+00,
	2.345879401782093864881062472704797983169555664062500000000000e+00,
	2.284469119227324540588597301393747329711914062500000000000000e+00,
	2.222192778101689913938798781600780785083770751953125000000000e+00,
	2.159236735690232222850681864656507968902587890625000000000000e+00,
	2.095781082670708173054663348011672496795654296875000000000000e+00,
	2.031998893208212031424864107975736260414123535156250000000000e+00,
	1.968055586986250382608432119013741612434387207031250000000000e+00,
	1.904108401542659478877794754225760698318481445312500000000000e+00,
	1.840305971835716647433400794398039579391479492187500000000000e+00,
	1.776788012722152476641213070251978933811187744140625000000000e+00,
	1.713685098989309985029194649541750550270080566406250000000000e+00,
	1.651118536748350740239743572601582854986190795898437500000000e+00,
	1.589200319355628110429279331583529710769653320312500000000000e+00,
	1.528033160577054738737956540717277675867080688476562500000000e+00,
	1.467710597430426755494181634276174008846282958984375000000000e+00,
	1.408317155018317778925052152771968394517898559570312500000000e+00,
	1.349928565680445169405743399693164974451065063476562500000000e+00,
	1.292612034933002140846269867324735969305038452148437500000000e+00,
	1.236426546902576140141150062845554202795028686523437500000000e+00,
	1.181423202287495977103048971912357956171035766601562500000000e+00,
	1.127645582270778268352273698837962001562118530273437500000000e+00,
	1.075130132251301162682466383557766675949096679687500000000000e+00,
	1.023906559736719179554143011046107858419418334960937500000000e+00,
	9.739982412412215850849861453752964735031127929687500000000000e-01,
	9.254226335396360525464842794463038444519042968750000000000000e-01,
	8.781916851381847965996030325186438858509063720703125000000000e-01,
	8.323122443210928533829928710474632680416107177734375000000000e-01,
	7.877864606150899451364466585800983011722564697265625000000000e-01,
	7.446121769748099694652410107664763927459716796875000000000000e-01,
	7.027833104263951913281971428659744560718536376953125000000000e-01,
	6.622902193114570401633045548805966973304748535156250000000000e-01,
	6.231200556473037721616492490284144878387451171875000000000000e-01,
	5.852571014597653942246324731968343257904052734375000000000000e-01,
	5.486830882528161978939351683948189020156860351562500000000000e-01,
	5.133774990545860728019533780752681195735931396484375000000000e-01,
	4.793178527222465845625265501439571380615234375000000000000000e-01,
	4.464799704014152936792925174813717603683471679687500000000000e-01,
	4.148382242195611802060284389881417155265808105468750000000000e-01,
	3.843657684497205906382077955640852451324462890625000000000000e-01,
	3.550347535126094911106520157773047685623168945312500000000000e-01,
	3.268165232941226960861058614682406187057495117187500000000000e-01,
	2.996817963428495268374263105215504765510559082031250000000000e-01,
	2.736008315814841829194392630597576498985290527343750000000000e-01,
	2.485435792181007741419307421892881393432617187500000000000000e-01,
	2.244798175807338758147579937940463423728942871093750000000000e-01,
	2.013792766230131370264189172303304076194763183593750000000000e-01,
	1.792117488618298892788516241125762462615966796875000000000000e-01,
	1.579471885113248319498779892455786466598510742187500000000000e-01,
	1.375557995729286275832237151917070150375366210937500000000000e-01,
	1.180081136294308663536867243237793445587158203125000000000000e-01,
	9.927505807402203785727579088415950536727905273437500000000000e-02,
	8.132801548335466534922488790471106767654418945312500000000000e-02,
	6.413887481848901117587047338020056486129760742187500000000000e-02,
	4.768007510939398230220831464976072311401367187500000000000000e-02,
	3.192464224872626488149762735702097415924072265625000000000000e-02,
	1.684621948916342759616782132070511579513549804687500000000000e-02,
	2.419092206420891955076513113453984260559082031250000000000000e-03,
	-1.138179254248683491823612712323665618896484375000000000000000e-02,
	-2.458081116852817160634003812447190284729003906250000000000000e-02,
	-3.720165294247967935348242463078349828720092773437500000000000e-02,
	-4.926731319550503584991929528769105672836303710937500000000000e-02,
	-6.080008992535157474890183948446065187454223632812500000000000e-02,
	-7.182158343446287673828010156285017728805541992187500000000000e-02,
	-8.235269865457839344458079722244292497634887695312500000000000e-02,
	-9.241364983953070932898299361113458871841430664062500000000000e-02 };