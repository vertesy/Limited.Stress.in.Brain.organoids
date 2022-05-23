# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0001666","response to hypoxia",1.48244180147681,7.82615259180876,0.939142145437665,0,"response to hypoxia"),
c("GO:0033555","multicellular organismal response to stress",0.388929598106082,0.740019298352358,0.891457911373118,0.28087355,"response to hypoxia"),
c("GO:0009612","response to mechanical stimulus",1.13296882926554,3.31265032258439,0.942521421097263,0.61008128,"response to hypoxia"),
c("GO:0070482","response to oxygen levels",1.6797249309509,6.54149359214967,0.951483854386118,0.63903314,"response to hypoxia"),
c("GO:0006816","calcium ion transport",1.34716194126599,1.59908055162017,0.985758964594704,0,"calcium ion transport"),
c("GO:0099003","vesicle-mediated transport in synapse",0.68767262273829,0.845478648711363,0.983817006253751,0.26070135,"calcium ion transport"),
c("GO:0006836","neurotransmitter transport",0.744039231159461,1.06364713712358,0.986489645395138,0.26293183,"calcium ion transport"),
c("GO:0030900","forebrain development",2.11938447663604,19.1719520261484,0.753794723093031,0,"forebrain development"),
c("GO:0007611","learning or memory",1.48244180147681,9.49647760045875,0.870747301972189,0.15551755,"forebrain development"),
c("GO:0001503","ossification",1.53317174905586,2.53514695389433,0.901140994697146,0.15616907,"forebrain development"),
c("GO:0045165","cell fate commitment",1.36407192379235,3.72782873349808,0.860918022139749,0.19603783,"forebrain development"),
c("GO:0021700","developmental maturation",1.38661856716081,1.29259079815714,0.882827818510304,0.19643302,"forebrain development"),
c("GO:0060560","developmental growth involved in morphogenesis",0.653852657685587,11.845515162968,0.837761387730946,0.20518773,"forebrain development"),
c("GO:0031099","regeneration",0.834225804633335,5.15914645804107,0.874485088642258,0.2109731,"forebrain development"),
c("GO:0009791","post-embryonic development",0.49038949326419,0.43530790252655,0.854537087562007,0.22405679,"forebrain development"),
c("GO:0014706","striated muscle tissue development",1.60081167916126,6.4571521744761,0.835928774386304,0.22818543,"forebrain development"),
c("GO:0048736","appendage development",0.997688969054732,1.25723638861055,0.844108032061167,0.24286361,"forebrain development"),
c("GO:0060173","limb development",0.997688969054732,1.25723638861055,0.833942617159077,0.24286361,"forebrain development"),
c("GO:0055123","digestive system development",0.783495857054281,1.45189324166974,0.839533975190043,0.25483643,"forebrain development"),
c("GO:0003002","regionalization",1.87137139958289,5.30738394360304,0.833549020405915,0.26236655,"forebrain development"),
c("GO:0007389","pattern specification process",2.49704075305789,4.83234009182832,0.828218838120948,0.27239764,"forebrain development"),
c("GO:0045444","fat cell differentiation",0.60875937094865,2.06530010304368,0.871685073387092,0.28285616,"forebrain development"),
c("GO:0001655","urogenital system development",1.89391804295136,10.0497861037957,0.82431193311472,0.28359753,"forebrain development"),
c("GO:0048880","sensory system development",2.11374781579392,5.33610410979385,0.822215815663474,0.28763647,"forebrain development"),
c("GO:0001649","osteoblast differentiation",0.749675892001578,1.24058298558833,0.8552060590877,0.28926624,"forebrain development"),
c("GO:0048608","reproductive structure development",2.37867087537343,5.67579988321462,0.811442158665272,0.29210972,"forebrain development"),
c("GO:0061458","reproductive system development",2.39558085789978,4.88245694093677,0.819765747065717,0.29689989,"forebrain development"),
c("GO:0007517","muscle organ development",1.59517501831915,2.9725688422331,0.807084897041699,0.36172652,"forebrain development"),
c("GO:0001822","kidney development",1.62899498337185,9.49647760045875,0.758864889491191,0.36267856,"forebrain development"),
c("GO:0003206","cardiac chamber morphogenesis",0.698945944422524,0.735276599098,0.790257331002782,0.37043778,"forebrain development"),
c("GO:0001654","eye development",2.05738120737275,5.48782809610141,0.770276692185289,0.37362159,"forebrain development"),
c("GO:0048732","gland development",2.26593765853109,1.06020677902028,0.799744688235144,0.37834269,"forebrain development"),
c("GO:0048568","embryonic organ development",2.53086071811059,0.806046579207435,0.797317468047776,0.38733688,"forebrain development"),
c("GO:0021675","nerve development",0.456569528211488,0.711049135192077,0.817540155520918,0.39947972,"forebrain development"),
c("GO:0003015","heart process",0.597486049264416,2.03991088428749,0.891623082407429,0.40848228,"forebrain development"),
c("GO:0006936","muscle contraction",1.27388535031847,2.50132286340554,0.882722589508245,0.44454939,"forebrain development"),
c("GO:0003012","muscle system process",1.54444507074009,5.15761670220189,0.883198323955291,0.45476382,"forebrain development"),
c("GO:0061564","axon development",2.12502113747816,18.1762923839648,0.735962128346539,0.47920506,"forebrain development"),
c("GO:0051216","cartilage development",0.884955752212389,3.27119348894803,0.797557527510529,0.49083393,"forebrain development"),
c("GO:0042698","ovulation cycle",0.383292937263965,0.910130198490448,0.905693525913019,0.49306989,"forebrain development"),
c("GO:0048762","mesenchymal cell differentiation",0.952595682317795,1.65476550127481,0.788571845938908,0.49494746,"forebrain development"),
c("GO:0061448","connective tissue development",1.17806211600248,1.98644235997517,0.843774367473725,0.50720783,"forebrain development"),
c("GO:0060485","mesenchyme development",1.29643199368694,1.87655018563006,0.792790422621213,0.51293442,"forebrain development"),
c("GO:0050954","sensory perception of mechanical stimulus",0.992052308212615,0.660769240324493,0.885928279511463,0.52241673,"forebrain development"),
c("GO:0060537","muscle tissue development",1.69663491347726,7.03113522699136,0.837768346147134,0.52974018,"forebrain development"),
c("GO:0060562","epithelial tube morphogenesis",1.71918155684572,1.34014462442308,0.778546032209309,0.53438282,"forebrain development"),
c("GO:0050890","cognition",1.71918155684572,7.65252361795182,0.879325327232864,0.55712925,"forebrain development"),
c("GO:0010001","glial cell differentiation",0.935685699791444,4.59788139919982,0.782347924328214,0.55882229,"forebrain development"),
c("GO:0001764","neuron migration",0.704582605264641,7.44139821665223,0.77330222473527,0.56448898,"forebrain development"),
c("GO:0042063","gliogenesis",1.20060875937095,4.57938258057523,0.779564361296223,0.57566074,"forebrain development"),
c("GO:0031290","retinal ganglion cell axon guidance",0.107096556000225,0.700403367418548,0.746888263116193,0.60068161,"forebrain development"),
c("GO:0001667","ameboidal-type cell migration",1.03714559494955,0.859991903242411,0.962569869819823,0.60397962,"forebrain development"),
c("GO:0021879","forebrain neuron differentiation",0.208556451158334,1.90902283459851,0.775701680771252,0.60638804,"forebrain development"),
c("GO:0008038","neuron recognition",0.276196381263739,3.89462949316971,0.797781526180754,0.60773574,"forebrain development"),
c("GO:0021988","olfactory lobe development",0.1972831294741,1.51696501455632,0.795326084310811,0.60839643,"forebrain development"),
c("GO:0021872","forebrain generation of neurons",0.270559720421622,1.96275283496721,0.773969249431161,0.62429725,"forebrain development"),
c("GO:0031102","neuron projection regeneration",0.169099825263514,1.37840313422818,0.766255839836633,0.62749136,"forebrain development"),
c("GO:0050767","regulation of neurogenesis",1.99537793810946,12.0897615096553,0.678825995906697,0.64015784,"forebrain development"),
c("GO:0021953","central nervous system neuron differentiation",0.935685699791444,4.53142376514512,0.768157981504067,0.64079345,"forebrain development"),
c("GO:0048565","digestive tract development",0.715855926948876,0.677445807252673,0.828736305263492,0.64466793,"forebrain development"),
c("GO:0021772","olfactory bulb development",0.186009807789865,1.06020677902028,0.79632313341782,0.67217984,"forebrain development"),
c("GO:0045137","development of primary sexual characteristics",1.25133870695,0.593098909060364,0.832295808448635,0.67229046,"forebrain development"),
c("GO:0060047","heart contraction",0.512936136632659,1.78779374531323,0.893148465930976,0.67459229,"forebrain development"),
c("GO:0007548","sex differentiation",1.53317174905586,0.778867558089622,0.871790774758417,0.68954044,"forebrain development"),
c("GO:0048511","rhythmic process",1.51062510568739,0.561346916284289,1,0,"rhythmic process"),
c("GO:0050804","modulation of chemical synaptic transmission",2.33357758863649,24.8910225247778,0.938369161710119,0,"modulation of chemical synaptic transmission"),
c("GO:2001236","regulation of extrinsic apoptotic signaling pathway",0.84549912631757,2.30160364307073,0.952646282227857,0.32509473,"modulation of chemical synaptic transmission"),
c("GO:0090092","regulation of transmembrane receptor protein serine/threonine kinase signaling pathway",1.43734851473987,2.5880067594686,0.952978693091897,0.34665742,"modulation of chemical synaptic transmission"),
c("GO:0099177","regulation of trans-synaptic signaling",2.33921424947861,25.2577329649729,0.952983424170117,0.36912841,"modulation of chemical synaptic transmission"),
c("GO:0030111","regulation of Wnt signaling pathway",1.83755143453018,0.8683845431552,0.951742439702257,0.36930583,"modulation of chemical synaptic transmission"),
c("GO:1902043","positive regulation of extrinsic apoptotic signaling pathway via death domain receptors",0.0789132517896398,0.401833487533222,0.950365803134791,0.6437939,"modulation of chemical synaptic transmission"),
c("GO:0050808","synapse organization",1.60644834000338,26.5137239959185,0.944571926526167,0,"synapse organization"),
c("GO:0140014","mitotic nuclear division",0.963869004002029,9.93282867613331,0.923235457231032,0.1952572,"synapse organization"),
c("GO:0051383","kinetochore organization",0.112733216842343,1.5786840531884,0.962167733508965,0.23914669,"synapse organization"),
c("GO:0048285","organelle fission",1.88828138210924,4.88908721400655,0.957641413275123,0.3174156,"synapse organization"),
c("GO:0008608","attachment of spindle microtubules to kinetochore",0.140916521052928,0.763155515840663,0.959105414542044,0.52282203,"synapse organization"),
c("GO:0051306","mitotic sister chromatid separation",0.0225466433684685,0.573724655509574,0.936734003053575,0.60938076,"synapse organization"),
c("GO:0099054","presynapse assembly",0.124006538526577,0.682791479027167,0.807136475064096,0.6510484,"synapse organization"),
c("GO:0099172","presynapse organization",0.140916521052928,0.936572810545152,0.953006701882665,0.65882949,"synapse organization"),
c("GO:0007052","mitotic spindle organization",0.518572797474776,1.71103573910386,0.932494729176346,0.69999203,"synapse organization"),
c("GO:0035249","synaptic transmission, glutamatergic",0.169099825263514,4.69051249193507,0.956502106784844,0.00442864,"synaptic transmission, glutamatergic"),
c("GO:0010644","cell communication by electrical coupling",0.124006538526577,0.511629807031991,0.977813925477655,0.13383208,"synaptic transmission, glutamatergic"),
c("GO:0019932","second-messenger-mediated signaling",1.24570204610789,4.33947887004059,0.918219646561975,0.16655162,"synaptic transmission, glutamatergic"),
c("GO:0035637","multicellular organismal signaling",0.665125979369821,1.60857967225797,0.887854379621507,0.19082063,"synaptic transmission, glutamatergic"),
c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",1.03714559494955,2.84910239901151,0.920969305628984,0.21716381,"synaptic transmission, glutamatergic"),
c("GO:0007188","adenylate cyclase-modulating G protein-coupled receptor signaling pathway",1.28515867200271,1.27992342723186,0.925422287866489,0.22249142,"synaptic transmission, glutamatergic"),
c("GO:0007215","glutamate receptor signaling pathway",0.259286398737388,0.825257144945945,0.932201249893933,0.36252839,"synaptic transmission, glutamatergic"),
c("GO:0007214","gamma-aminobutyric acid signaling pathway",0.152189842737163,1.01361796237566,0.960270688885421,0.41713283,"synaptic transmission, glutamatergic"),
c("GO:0070371","ERK1 and ERK2 cascade",0.208556451158334,0.492678558284518,0.930552077318819,0.4348791,"synaptic transmission, glutamatergic"),
c("GO:0019935","cyclic-nucleotide-mediated signaling",0.281833042105856,0.520863956700953,0.928733365096534,0.44764534,"synaptic transmission, glutamatergic"),
c("GO:0048015","phosphatidylinositol-mediated signaling",0.41147624147455,0.876105911270818,0.926311825303519,0.46479397,"synaptic transmission, glutamatergic"),
c("GO:0048017","inositol lipid-mediated signaling",0.434022884843019,0.663478228513781,0.925957590723683,0.46731755,"synaptic transmission, glutamatergic"),
c("GO:0023061","signal release",0.890592413054507,3.21113477071556,0.940979052542723,0.483506,"synaptic transmission, glutamatergic"),
c("GO:0051932","synaptic transmission, GABAergic",0.107096556000225,0.811259891911177,0.957984550048677,0.55830298,"synaptic transmission, glutamatergic"),
c("GO:0016055","Wnt signaling pathway",1.45989515810834,0.790568099144049,0.902142618026443,0.59998773,"synaptic transmission, glutamatergic"),
c("GO:0198738","cell-cell signaling by wnt",1.45989515810834,0.72433398209028,0.951875175863343,0.6350655,"synaptic transmission, glutamatergic"),
c("GO:0030509","BMP signaling pathway",0.439659545685136,0.962613347417623,0.89425990522277,0.65508043,"synaptic transmission, glutamatergic"),
c("GO:0006029","proteoglycan metabolic process",0.428386224000902,2.99035367316813,0.976800133417806,0.00485817,"proteoglycan metabolic process"),
c("GO:1903510","mucopolysaccharide metabolic process",0.462206189053605,0.946406857555906,0.972129345011715,0.51971584,"proteoglycan metabolic process"),
c("GO:0030166","proteoglycan biosynthetic process",0.293106363790091,0.582892170430238,0.975933924374572,0.69131955,"proteoglycan metabolic process"),
c("GO:0050673","epithelial cell proliferation",0.518572797474776,3.66102186954443,0.988205568572776,0.00495699,"epithelial cell proliferation"),
c("GO:0033002","muscle cell proliferation",0.107096556000225,1.39280455870425,0.989180614247468,0.5736012,"epithelial cell proliferation"),
c("GO:0061351","neural precursor cell proliferation",0.394566258948199,1.75223043361677,0.98837118530035,0.64288058,"epithelial cell proliferation"),
c("GO:0098742","cell-cell adhesion via plasma-membrane adhesion molecules",1.4486218364241,5.95815746089915,0.986578209597046,0.00556569,"cell-cell adhesion via plasma-membrane adhesion molecules"),
c("GO:0031589","cell-substrate adhesion",1.03714559494955,1.38969337726221,0.987915293255946,0.62993776,"cell-cell adhesion via plasma-membrane adhesion molecules"),
c("GO:0007059","chromosome segregation",1.60644834000338,5.86662078807648,0.996774563886741,0.00563535,"chromosome segregation"),
c("GO:0050803","regulation of synapse structure or activity",1.1836987768446,12.940209460426,0.950496472788316,0.02706311,"regulation of synapse structure or activity"),
c("GO:0009914","hormone transport",0.473479510737839,0.499221681531868,0.943765552337645,0.25676263,"regulation of synapse structure or activity"),
c("GO:0060078","regulation of postsynaptic membrane potential",0.52984611915901,4.02784114889728,0.951019101886526,0.25974708,"regulation of synapse structure or activity"),
c("GO:0042391","regulation of membrane potential",2.35048757116284,9.93282867613331,0.946684032039795,0.30701293,"regulation of synapse structure or activity"),
c("GO:0001505","regulation of neurotransmitter levels",1.20060875937095,0.55249002256169,0.950422892156154,0.30754577,"regulation of synapse structure or activity"),
c("GO:0072503","cellular divalent inorganic cation homeostasis",2.57595400484753,3.93100939795285,0.926487088660103,0.33923183,"regulation of synapse structure or activity"),
c("GO:0051098","regulation of binding",2.05174454653064,3.93100939795285,0.966177517681021,0.02901046,"regulation of binding"),
c("GO:0045926","negative regulation of growth",1.34152528042388,2.99724663354876,0.942807455493242,0.0361331,"negative regulation of growth"),
c("GO:0040013","negative regulation of locomotion",1.78682148695113,0.907409171935913,0.95597387400369,0.20956631,"negative regulation of growth"),
c("GO:0051983","regulation of chromosome segregation",0.541119440843244,0.520863956700953,0.975527401078434,0.03857399,"regulation of chromosome segregation"),
c("GO:0051960","regulation of nervous system development",2.4350374837946,18.0725824112872,0.912048603229269,0.03901561,"regulation of nervous system development"),
c("GO:0086064","cell communication by electrical coupling involved in cardiac conduction",0.0958232343159912,1.06020677902028,0.865197551881188,0.30651257,"regulation of nervous system development"),
c("GO:0045664","regulation of neuron differentiation",1.04841891663379,9.93282867613331,0.923583078743264,0.42151012,"regulation of nervous system development"),
c("GO:0045669","positive regulation of osteoblast differentiation",0.372019615579731,0.371901377858339,0.918594100219831,0.44687461,"regulation of nervous system development"),
c("GO:0022604","regulation of cell morphogenesis",1.71354489600361,0.69577452442712,0.92908386150253,0.44812087,"regulation of nervous system development"),
c("GO:0045665","negative regulation of neuron differentiation",0.405839580632433,2.2614461580213,0.908541435640664,0.45073787,"regulation of nervous system development"),
c("GO:0045667","regulation of osteoblast differentiation",0.693309283580407,1.76379538833324,0.926669336274852,0.4760679,"regulation of nervous system development"),
c("GO:0045598","regulation of fat cell differentiation",0.749675892001578,0.316989905394165,0.92610442439144,0.48000516,"regulation of nervous system development"),
c("GO:0031644","regulation of nervous system process",0.772222535370047,0.859991903242411,0.948203450961432,0.57149731,"regulation of nervous system development"),
c("GO:0010810","regulation of cell-substrate adhesion",1.1836987768446,0.877047265986805,0.973456162928944,0.04225877,"regulation of cell-substrate adhesion"),
c("GO:0090287","regulation of cellular response to growth factor stimulus",1.63463164421397,4.38947717231886,0.965889143555046,0.0439916,"regulation of cellular response to growth factor stimulus"),
c("GO:0050678","regulation of epithelial cell proliferation",1.90519136463559,3.77655078779577,0.962182702214381,0.04486458,"regulation of epithelial cell proliferation"),
c("GO:2000177","regulation of neural precursor cell proliferation",0.473479510737839,0.339915626536469,0.966924738988867,0.46570533,"regulation of epithelial cell proliferation"),
c("GO:0048660","regulation of smooth muscle cell proliferation",0.760949213685813,1.37628564318403,0.965431877431032,0.49070351,"regulation of epithelial cell proliferation"),
c("GO:0050679","positive regulation of epithelial cell proliferation",1.07660222084437,0.852772111467513,0.951040267495026,0.51075468,"regulation of epithelial cell proliferation"),
c("GO:0010975","regulation of neuron projection development",2.38994419705766,18.9213402202731,0.904676539631123,0.04622212,"regulation of neuron projection development"),
c("GO:0051783","regulation of nuclear division",0.783495857054281,3.27612829421319,0.918723703249037,0.41768699,"regulation of neuron projection development"),
c("GO:1901888","regulation of cell junction assembly",1.11042218589707,6.27566103936872,0.918776888320597,0.43537634,"regulation of neuron projection development"),
c("GO:0050807","regulation of synapse organization",1.14987881179189,11.4435069340255,0.886532287233872,0.43723034,"regulation of neuron projection development"),
c("GO:1901890","positive regulation of cell junction assembly",0.580576066738064,1.50346748413527,0.906882375192424,0.57135026,"regulation of neuron projection development"),
c("GO:0007346","regulation of mitotic cell cycle",2.53649737895271,6.8660365847108,0.898568446595363,0.04674236,"regulation of mitotic cell cycle"),
c("GO:1905818","regulation of chromosome separation",0.439659545685136,1.6443952490412,0.890321814078249,0.5877477,"regulation of mitotic cell cycle"),
c("GO:1904029","regulation of cyclin-dependent protein kinase activity",0.552392762527479,1.45189324166974,0.896847708702054,0.60286747,"regulation of mitotic cell cycle"),
c("GO:0007088","regulation of mitotic nuclear division",0.625669353475001,3.21360566926816,0.863576760080601,0.61145122,"regulation of mitotic cell cycle"),
c("GO:0000079","regulation of cyclin-dependent protein serine/threonine kinase activity",0.535482780001127,1.26730876993149,0.897147599668294,0.61512871,"regulation of mitotic cell cycle"),
c("GO:0045787","positive regulation of cell cycle",1.8093681303196,0.529789501068895,0.888603965266758,0.695924,"regulation of mitotic cell cycle"),
c("GO:0034765","regulation of ion transmembrane transport",2.64923059579505,3.04143229462119,0.947292436358974,0.04740139,"regulation of ion transmembrane transport"),
c("GO:0051384","response to glucocorticoid",0.783495857054281,7.08150781535204,0.940774268615489,0.09728083,"response to glucocorticoid"),
c("GO:0050919","negative chemotaxis",0.259286398737388,0.541106936352026,0.935339129466718,0.21610723,"response to glucocorticoid"),
c("GO:0010038","response to metal ion",1.97846795558311,4.97385094753624,0.933390710659767,0.26617116,"response to glucocorticoid"),
c("GO:0009410","response to xenobiotic stimulus",2.24902767600473,2.94256974170973,0.947785473930379,0.3025396,"response to glucocorticoid"),
c("GO:1901654","response to ketone",1.06532889916014,3.45128268128024,0.946337359378681,0.34755845,"response to glucocorticoid"),
c("GO:0071774","response to fibroblast growth factor",0.462206189053605,1.24018867301212,0.945532695768831,0.46985512,"response to glucocorticoid"),
c("GO:0097305","response to alcohol",1.30770531537117,2.13611232956302,0.945250064155302,0.50976358,"response to glucocorticoid"),
c("GO:0048545","response to steroid hormone",1.53880840989798,6.16858650651211,0.938545194971564,0.67382143,"response to glucocorticoid"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PostGruffirevigo_treemap.pdf", width=16, height=16 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Enriched GO-terms post-Gruffi",
  inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

