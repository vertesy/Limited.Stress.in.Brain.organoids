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
revigo.data <- rbind(c("GO:0001666","response to hypoxia",0.0389767702923432,8.23066846929386,0.901939232051765,0,"response to hypoxia"),
c("GO:0048545","response to steroid hormone",0.047888232122576,6.54241328168805,0.875254736601385,0.22121132,"response to hypoxia"),
c("GO:0071496","cellular response to external stimulus",0.451389777894553,2.29176866700372,0.904516216622744,0.26274768,"response to hypoxia"),
c("GO:0048678","response to axon injury",0.00861731315987179,0.664574196162756,0.935391408379281,0.33670817,"response to hypoxia"),
c("GO:0010038","response to metal ion",0.118761489923502,4.14010370078061,0.880883301855707,0.44818848,"response to hypoxia"),
c("GO:0009410","response to xenobiotic stimulus",0.0983823728882093,2.89958464903298,0.892266309338303,0.47181184,"response to hypoxia"),
c("GO:0050919","negative chemotaxis",0.0066701318208623,0.507651757833093,0.887673775963951,0.52202078,"response to hypoxia"),
c("GO:1901654","response to ketone",0.0219285762284622,3.45412038825901,0.888884297124241,0.55576868,"response to hypoxia"),
c("GO:0009612","response to mechanical stimulus",0.0290337166463373,5.72035586613469,0.901000904915214,0.6351101,"response to hypoxia"),
c("GO:0031667","response to nutrient levels",0.118069619107301,1.47330800161712,0.912416791681024,0.64556441,"response to hypoxia"),
c("GO:0070482","response to oxygen levels",0.043749436042426,6.90598681647345,0.913276020505651,0.65181965,"response to hypoxia"),
c("GO:0071560","cellular response to transforming growth factor beta stimulus",0.0337773818232859,1.68734670435434,0.866806852766842,0.65400239,"response to hypoxia"),
c("GO:0097305","response to alcohol",0.0770959522077183,1.52266454616737,0.880512827734105,0.66153339,"response to hypoxia"),
c("GO:0006816","calcium ion transport",0.225015446948133,0.960205296304466,0.995457603644356,0,"calcium ion transport"),
c("GO:0051651","maintenance of location in cell",0.0478592315494418,0.508097875879915,0.991406930568322,0.19757193,"calcium ion transport"),
c("GO:0030900","forebrain development",0.0629312437011791,18.5684689842499,0.557796663938344,0,"forebrain development"),
c("GO:0007611","learning or memory",0.0377794609158033,7.75113207676207,0.672740549445922,0.46736878,"forebrain development"),
c("GO:0001503","ossification",0.0455060421865537,3.21268407820508,0.684558852232777,0.47314125,"forebrain development"),
c("GO:0031099","regeneration",0.0301274525473979,4.11841792964102,0.671677131815256,0.48388327,"forebrain development"),
c("GO:0060560","developmental growth involved in morphogenesis",0.0402486525712281,10.647354481331,0.635110908211931,0.49311063,"forebrain development"),
c("GO:0045165","cell fate commitment",0.0740633208456866,6.06027113856749,0.645837352516532,0.49643632,"forebrain development"),
c("GO:0014706","striated muscle tissue development",0.0598986123391473,6.40864929829986,0.636379706141099,0.50636513,"forebrain development"),
c("GO:0048736","appendage development",0.0463884881976368,1.25334790040431,0.614565113025925,0.54432547,"forebrain development"),
c("GO:0055123","digestive system development",0.0342289621763754,1.40702425416595,0.615033745101213,0.54705315,"forebrain development"),
c("GO:0045444","fat cell differentiation",0.0164433249670823,2.05019102813787,0.674180211187617,0.550163,"forebrain development"),
c("GO:0001649","osteoblast differentiation",0.0164474679061014,0.842812293684679,0.64979180176725,0.55017171,"forebrain development"),
c("GO:0060541","respiratory system development",0.0393744924381834,1.55222199070538,0.612080434810397,0.55208991,"forebrain development"),
c("GO:0001655","urogenital system development",0.0671694703177891,10.931791867313,0.60040166409242,0.5721812,"forebrain development"),
c("GO:0003015","heart process",0.0212491342293185,0.483897152817864,0.690827925779954,0.58721922,"forebrain development"),
c("GO:0001763","morphogenesis of a branching structure",0.036793441429241,1.06667167884898,0.626642946009592,0.59075552,"forebrain development"),
c("GO:0035107","appendage morphogenesis",0.0376593156842474,1.7926753764288,0.603580194196211,0.59163065,"forebrain development"),
c("GO:0048608","reproductive structure development",0.148317216886255,6.22272328528146,0.572314517952297,0.60767777,"forebrain development"),
c("GO:0003002","regionalization",0.100105835520184,4.83853915852377,0.597434137254626,0.61029571,"forebrain development"),
c("GO:0003012","muscle system process",0.0463180582343109,4.016416107137,0.677192610232858,0.6166288,"forebrain development"),
c("GO:0006936","muscle contraction",0.0389063403290173,2.07846386323552,0.678137257507149,0.61779588,"forebrain development"),
c("GO:0007389","pattern specification process",0.139086748751546,5.1486609167129,0.589652323996862,0.62562698,"forebrain development"),
c("GO:0048880","sensory system development",0.108097564888161,5.31592161166456,0.589409148363091,0.62934699,"forebrain development"),
c("GO:0050890","cognition",0.043774293676541,6.22272328528146,0.677785458845214,0.63702157,"forebrain development"),
c("GO:0021675","nerve development",0.0181253582088649,1.92088028336241,0.610247963917522,0.6384626,"forebrain development"),
c("GO:0001822","kidney development",0.0578685722197544,10.4499356087666,0.550960081829903,0.64153934,"forebrain development"),
c("GO:0061458","reproductive system development",0.148884799531881,6.43875069507304,0.581679272435921,0.64481756,"forebrain development"),
c("GO:0050954","sensory perception of mechanical stimulus",0.0529757612381157,1.06500794564352,0.674285892626738,0.65111842,"forebrain development"),
c("GO:0090596","sensory organ morphogenesis",0.0741171790529357,6.65662883000878,0.554200686082447,0.65242985,"forebrain development"),
c("GO:0007517","muscle organ development",0.0642321265531982,3.0817745892221,0.57567019675429,0.65334684,"forebrain development"),
c("GO:0048568","embryonic organ development",0.120580240152918,4.1038724361218,0.566223540044198,0.68289084,"forebrain development"),
c("GO:0048762","mesenchymal cell differentiation",0.0359524248083497,2.91147056550276,0.576618725546415,0.68529782,"forebrain development"),
c("GO:0048732","gland development",0.0804475898742262,3.34463774983014,0.576841866617654,0.6869331,"forebrain development"),
c("GO:0051216","cartilage development",0.0417069671059756,4.73419669218065,0.579458543251758,0.69199061,"forebrain development"),
c("GO:0060485","mesenchyme development",0.0485801029387772,4.69590172071705,0.577580356490406,0.69900459,"forebrain development"),
c("GO:0050804","modulation of chemical synaptic transmission",0.0782808327672007,23.8863239881029,0.841437732956799,0,"modulation of chemical synaptic transmission"),
c("GO:0050803","regulation of synapse structure or activity",0.0410979550701578,10.5845208266048,0.876145217836766,0.1489172,"modulation of chemical synaptic transmission"),
c("GO:0051098","regulation of binding",0.0596334642419205,3.53600703032315,0.884131021055436,0.15272007,"modulation of chemical synaptic transmission"),
c("GO:0040013","negative regulation of locomotion",0.0511984403988921,1.87894269176976,0.855841108553822,0.16198096,"modulation of chemical synaptic transmission"),
c("GO:0051960","regulation of nervous system development",0.0716397015194726,16.8410287161565,0.839842548679803,0.16576144,"modulation of chemical synaptic transmission"),
c("GO:0010959","regulation of metal ion transport",0.0762632214648653,1.25575944838503,0.865184173406813,0.16648482,"modulation of chemical synaptic transmission"),
c("GO:0010810","regulation of cell-substrate adhesion",0.0350824076143242,0.873449767171385,0.893695575727177,0.17138616,"modulation of chemical synaptic transmission"),
c("GO:0050678","regulation of epithelial cell proliferation",0.052665040811678,4.12188741530969,0.886758455140535,0.17612211,"modulation of chemical synaptic transmission"),
c("GO:0010975","regulation of neuron projection development",0.0629685301523516,19.3441482557237,0.827034662102146,0.17828903,"modulation of chemical synaptic transmission"),
c("GO:0090287","regulation of cellular response to growth factor stimulus",0.0685407831331341,6.56719420128714,0.865412606354019,0.17933611,"modulation of chemical synaptic transmission"),
c("GO:0007346","regulation of mitotic cell cycle",0.187049553776467,4.12930178728912,0.852159770705674,0.19273801,"modulation of chemical synaptic transmission"),
c("GO:1900120","regulation of receptor binding",0.00377007450744391,1.86540856528912,0.898217636661931,0.43186369,"modulation of chemical synaptic transmission"),
c("GO:0042391","regulation of membrane potential",0.105470941550008,8.07016588291466,0.868633710025167,0.47098707,"modulation of chemical synaptic transmission"),
c("GO:0060078","regulation of postsynaptic membrane potential",0.0486795334752373,4.01044651579887,0.874857927960821,0.47649303,"modulation of chemical synaptic transmission"),
c("GO:0055074","calcium ion homeostasis",0.116507731097074,2.37964529062717,0.858041041793604,0.50704789,"modulation of chemical synaptic transmission"),
c("GO:0032409","regulation of transporter activity",0.052739613714023,0.38102152204197,0.885008544604124,0.50796265,"modulation of chemical synaptic transmission"),
c("GO:0033674","positive regulation of kinase activity",0.113669817868943,0.734993394339413,0.840168512616746,0.53542471,"modulation of chemical synaptic transmission"),
c("GO:0045664","regulation of neuron differentiation",0.037477026367404,6.90598681647345,0.848028323885379,0.5663978,"modulation of chemical synaptic transmission"),
c("GO:2001236","regulation of extrinsic apoptotic signaling pathway",0.0246546301030755,3.79390305529343,0.849818977881905,0.58654343,"modulation of chemical synaptic transmission"),
c("GO:0030510","regulation of BMP signaling pathway",0.0271859658437878,0.585812433262909,0.853829007700851,0.59027493,"modulation of chemical synaptic transmission"),
c("GO:1901888","regulation of cell junction assembly",0.0377546032816883,6.23819765925774,0.839611084010033,0.60982193,"modulation of chemical synaptic transmission"),
c("GO:0050807","regulation of synapse organization",0.0395029235477777,10.1401348829945,0.814167616063508,0.61163723,"modulation of chemical synaptic transmission"),
c("GO:0030178","negative regulation of Wnt signaling pathway",0.0563273989046235,0.702721089766887,0.813050642819932,0.61965807,"modulation of chemical synaptic transmission"),
c("GO:0090092","regulation of transmembrane receptor protein serine/threonine kinase signaling pathway",0.0603377638751792,3.83969056925386,0.846662048135204,0.62258402,"modulation of chemical synaptic transmission"),
c("GO:0099177","regulation of trans-synaptic signaling",0.0783844062426799,23.2408195938654,0.850722999022684,0.63397318,"modulation of chemical synaptic transmission"),
c("GO:0001558","regulation of cell growth",0.0786329825838301,6.126917831412,0.821458360484402,0.64063909,"modulation of chemical synaptic transmission"),
c("GO:0051783","regulation of nuclear division",0.0677453388414536,2.74718093751468,0.835039390243669,0.64388193,"modulation of chemical synaptic transmission"),
c("GO:0030111","regulation of Wnt signaling pathway",0.0992026748140048,1.19976765037011,0.841839121887445,0.64674368,"modulation of chemical synaptic transmission"),
c("GO:0051271","negative regulation of cellular component movement",0.0483232407195887,0.889959396384282,0.836631492933262,0.66291436,"modulation of chemical synaptic transmission"),
c("GO:1901890","positive regulation of cell junction assembly",0.0168410471129225,2.05901361571145,0.825228650011543,0.66398224,"modulation of chemical synaptic transmission"),
c("GO:0090257","regulation of muscle system process",0.0407085188023559,1.96989029932178,0.868180089196678,0.68708604,"modulation of chemical synaptic transmission"),
c("GO:0050808","synapse organization",0.05749156476901,25.0945884109352,0.924584324478273,0,"synapse organization"),
c("GO:0140014","mitotic nuclear division",0.139604616128942,7.77192060390262,0.89322504140156,0.34515148,"synapse organization"),
c("GO:0030198","extracellular matrix organization",0.0637598315050129,1.75054725726032,0.936924069518868,0.34771524,"synapse organization"),
c("GO:0043062","extracellular structure organization",0.0654750082589489,2.05882559931908,0.936817061655706,0.34837902,"synapse organization"),
c("GO:0051383","kinetochore organization",0.0243190520425228,1.59213002113856,0.934382627314875,0.41094704,"synapse organization"),
c("GO:0045229","external encapsulating structure organization",1.03756179089974,1.59708208239495,0.923454098973737,0.43476586,"synapse organization"),
c("GO:0048285","organelle fission",0.274088559630191,4.81852800705369,0.925274070962433,0.49074705,"synapse organization"),
c("GO:0035249","synaptic transmission, glutamatergic",0.00565511176116586,5.34803812313582,0.96852316192097,0.00631746,"synaptic transmission, glutamatergic"),
c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.0792544234367055,4.3044430289778,0.809954488642578,0.28568846,"synaptic transmission, glutamatergic"),
c("GO:0019932","second-messenger-mediated signaling",0.0821751954452197,4.27366815441415,0.813775616060893,0.34247203,"synaptic transmission, glutamatergic"),
c("GO:0019935","cyclic-nucleotide-mediated signaling",0.00551839477353328,0.483897152817864,0.840643516742705,0.37042008,"synaptic transmission, glutamatergic"),
c("GO:0048015","phosphatidylinositol-mediated signaling",0.049578551242397,0.674033036298575,0.819459698160631,0.42569245,"synaptic transmission, glutamatergic"),
c("GO:0097191","extrinsic apoptotic signaling pathway",0.0147571487862804,2.67830424214533,0.814256285369624,0.57862275,"synaptic transmission, glutamatergic"),
c("GO:0070997","neuron death",0.0102744887675394,1.13631022684037,0.980709759593408,0.65739968,"synaptic transmission, glutamatergic"),
c("GO:0051402","neuron apoptotic process",0.00854274025752675,1.11634480964854,0.980563000382439,0.66948443,"synaptic transmission, glutamatergic"),
c("GO:0016055","Wnt signaling pathway",0.122854713674441,3.99396204500268,0.796687097105806,0.68328751,"synaptic transmission, glutamatergic"),
c("GO:0050673","epithelial cell proliferation",0.0169446205884017,3.63769097608044,0.984973129571601,0.00674689,"epithelial cell proliferation"),
c("GO:0048659","smooth muscle cell proliferation",0.000546867950530325,0.38102152204197,0.986565047538728,0.66532961,"epithelial cell proliferation"),
c("GO:0006029","proteoglycan metabolic process",0.0323770684348068,1.0796642763117,0.995643686237139,0.00702879,"proteoglycan metabolic process"),
c("GO:0098742","cell-cell adhesion via plasma-membrane adhesion molecules",0.167416165764625,3.24440104704001,0.989436676679469,0.0078624,"cell-cell adhesion via plasma-membrane adhesion molecules"),
c("GO:0007059","chromosome segregation",0.444595357903116,2.88098508559726,0.994799361104687,0.0092244,"chromosome segregation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PREGruffirevigo_treemap.pdf", width=16, height=16 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap Pre-Gruffi",
  inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

