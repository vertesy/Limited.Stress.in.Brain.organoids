# RUn WGCNA on Retina integration
# Oliver Eichmueller
# Sun May 22 14:58:32 2022 ------------------------------

library(gridExtra)
library(igraph);
library(RColorBrewer);
library(tidyverse)
library(WGCNA)
library(gplots)
library(scWGCNA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(UCell)

OutDir <- 'Gruffi_revision/Sridhar2020/WGCNA/'
dir.create(OutDir)

sc.obj.integration <- ReadRDS(file = "gruffi_revision/Sridhar2020/d205_d125P_d125C_integration_progeny.Rds")
meta.data <- read.csv('Gruffi_revision/Sridhar2020/paper metadata/F6/cca_fetalvsorg_125CP_205_metadata.csv', row.names = 1)
colors_ret <- colorRampPalette(c(pals::tol(12), "grey"))(14)
names(colors_ret) <- meta.data$type %>% unique

DimPlot(sc.obj.integration, group.by = 'type', label = T)

# Construct meta cells for "type" annotation
sc.obj.integration$metacell_group <- paste0(
  "Meta_", 
  as.character(sc.obj.integration$type)
)

sc.obj.integration@active.assay <- 'integrated'
genes.keep <- VariableFeatures(sc.obj.integration)

# use k=20 cells per granule
seurat_list <- list()
for(group in unique(sc.obj.integration$metacell_group)){
  print(group)
  cur_seurat <- subset(sc.obj.integration, metacell_group == group)
  cur_seurat <- cur_seurat[genes.keep,]
  cur_metacell_seurat <- scWGCNA::construct_metacells(
    cur_seurat, name=group,
    k=20, reduction='umap',
    assay='RNA', slot='data'
  )
  cur_metacell_seurat$type <- as.character(unique(cur_seurat$type))
  
  seurat_list[[group]] <- cur_metacell_seurat
}

# merge all of the metacells objects
metacell_seurat <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])

# re-run dim reduc
metacell_seurat <- ScaleData(metacell_seurat, features = rownames(metacell_seurat))
metacell_seurat <- RunPCA(metacell_seurat, features=rownames(metacell_seurat))
metacell_seurat <- RunUMAP(metacell_seurat, dims=1:15)

pl1 <- DimPlot(sc.obj.integration, group.by='type', reduction='umap', label=TRUE)
pl2 <- DimPlot(metacell_seurat, group.by='type', reduction='umap', label=TRUE)
pl1+pl2

saveRDS(metacell_seurat, file= paste0(OutDir, 'metacell_Sridhar_type_filtered_k20.rds'))

# Run WGCNA
# enableWGCNAThreads(nThreads = 1)
disableWGCNAThreads()

# how many groups are there
nclusters <- length(unique(metacell_seurat$type))

# which genes are we using ?
genes.use <- rownames(metacell_seurat)

# cell meta-data table
targets <- metacell_seurat@meta.data[,]

# format the expression matrix for WGCNA
datExpr <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[genes.use,])
datExpr <- as.data.frame(t(datExpr))

# only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));

# Call the network topology analysis function for each set in turn
powerTable = list(
  data = pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose = 100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
)

# Plot the results:
pdf(paste0(OutDir,"1_Power_k20.pdf"), height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (col in 1:length(plotCols)){
  ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
  ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;

for (col in 1:length(plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
       xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
       main = colNames[col]);
  addGrid();
  
  if (col==1){
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
         labels=powers,cex=cex1,col=colors[1]);
  } else
    text(powerTable$data[,1], powerTable$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[1]);
  if (col==1){
    legend("bottomright", legend = 'Metacells', col = colors, pch = 20) ;
  } else
    legend("topright", legend = 'Metacells', col = colors, pch = 20) ;
}
dev.off()


softPower=10

nSets = 1
setLabels = 'Meta'
shortLabels = setLabels

multiExpr <- list()
multiExpr[['Meta']] <- list(data=datExpr)

checkSets(multiExpr) # check data size

# construct network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 10000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")


consMEs = net$multiMEs;
moduleLabels = net$colors;

# Convert the numeric labels to color labels
moduleColors = as.character(moduleLabels)
consTree = net$dendrograms[[1]];

# module eigengenes
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"

# intramodular connectivity
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

# compile into a module metadata table
geneInfo=as.data.frame(cbind(colnames(datExpr),moduleColors, KMEs))

# how many modules did we get?
nmodules <- length(unique(moduleColors))

# merged gene symbol column
colnames(geneInfo)[1]= "GeneSymbol"
colnames(geneInfo)[2]= "Initially.Assigned.Module.Color"

# save info
write.csv(geneInfo,file=paste0(OutDir,'geneInfoSigned_type_k20.csv'))

PCvalues=MEs

# Visualisation
pdf(paste0(OutDir,"SignedDendro_k20.pdf"),height=5, width=8)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("Gene dendrogram and module colors"))
dev.off()

plot_df <- cbind(dplyr::select(targets, c(type)), PCvalues)
plot_df <- reshape2::melt(plot_df, id.vars = c('type'))
plot_df$type <- factor(plot_df$type)

colors <- sub('ME', '', as.character(levels(plot_df$variable)))
p <- ggplot(plot_df, aes(x=type, y=value, fill=type)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('Module Eigengene') + xlab('') +
  facet_grid(~variable, scales='free_y') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )

# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf(paste0(OutDir,'ME_trajectory_Plot_modulevsCondition_k20.pdf'),width=w,height=8,useDingbats=F)
p 
dev.off()

p <- ggplot(plot_df, aes(x=type, y=value, fill=type)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('Module Eigengene') + xlab('') +
  facet_grid(variable~type, scales='free') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )
# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf(paste0(OutDir, 'ME_trajectory_Plot_modulevstype_k20.pdf'),width=w,height=8,useDingbats=F)
p 
dev.off()


write.csv(plot_df, file = paste0(OutDir,'plot_df_k20.csv'))

# Network Plots

load("ConsensusTOM-block.1.rda")

TOM.matrix = as.matrix(datExpr);
uniquemodcolors = unique(geneInfo$Initially.Assigned.Module.Color);
uniquemodcolors=uniquemodcolors

pdf(paste0(OutDir, 'ModuleNetworks_type_k20.pdf'),height=9,width=10);

for (mod in uniquemodcolors)  {
  numgenesingraph = 25;
  numconnections2keep = 500;
  cat('module:',mod,'\n');
  geneInfo=geneInfo[geneInfo$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo)==paste('kME',mod,sep=''));
  rowind = which(geneInfo$Initially.Assigned.Module.Color==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  matchind = match(submatrix$GeneSymbol,colnames(datExpr));
  reducedTOM = TOM.matrix[matchind,matchind];
  
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  gB <- graph.adjacency(as.matrix(reducedTOM[11:25,11:25]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))
  
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  
  plot(g1,
       edge.color=adjustcolor(mod, alpha.f=0.25),
       edge.alpha=0.25,
       vertex.color=adjustcolor(mod, alpha.f=0.75),
       vertex.label=as.character(submatrix$GeneSymbol),
       vertex.label.cex=2.2,
       vertex.label.dist=1.1,
       vertex.label.degree=-pi/4,
       vertex.label.color="black",
       #vertex.frame.color='black',
       layout= jitter(layoutCircle),
       vertex.size=6,
       main=paste(mod,"module")
  )
}
dev.off();

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
geneTree = hclust(as.dist(dissTOM), method = "average")
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^10;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
pdf(paste0(OutDir, "TOMplot_k20.pdf"))
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col = myheatcol)
dev.off()
png(paste0(OutDir, "TOMplot_k20.png"))
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col = myheatcol)
dev.off()

saveRDS(datExpr, paste0(OutDir,'datExpr_WGCNA_type_k20_220521.Rds'))
saveRDS(net, paste0(OutDir,'net_WGCNA_type_k20_220521.Rds'))
saveRDS(metacell_seurat, paste0(OutDir,'metacell_seurat_WGCNA_type_k20_220521.Rds'))
write.csv(geneInfo,file= paste0(OutDir,'geneInfo_signed_WGCNA_type_k20_220521.Rds'))

# Use modules for plotting in seurat -------------------------------------------

gene_mod_split = list()

for (i in 1:length(uniquemodcolors)) {
  color_oi = uniquemodcolors[i]
  gene_mod_split[[color_oi]] <- geneInfo %>% filter(Initially.Assigned.Module.Color == color_oi) %>%
    dplyr::arrange(dplyr::desc(paste0("kME", color_oi))) %>%
    dplyr::slice(1:30) %>% magrittr::use_series(GeneSymbol)
}


names(gene_mod_split) <- paste0(names(gene_mod_split), "_type_k20")

sc.obj.integration <- UCell::AddModuleScore_UCell(sc.obj.integration, features = gene_mod_split, assay = "RNA")

png(paste0(OutDir, "Score_UMAPs_type_k20.png"), width = 5000, height = 3000, res = 200)
FeaturePlot(sc.obj.integration, features = colnames(sc.obj.integration@meta.data)[62:72]
            , pt.size = 1, cols = c(alpha("grey", .5), "red")) & NoAxes()
dev.off()

png(paste0(OutDir, "Score_UMAPs_type_ident_k20.png"), width = 5000, height = 8000, res = 200)
FeaturePlot(sc.obj.integration, features = colnames(sc.obj.integration@meta.data)[62:72]
            , pt.size = 1, cols = c(alpha("grey", .5), "red"), split.by = "orig.ident") & NoAxes()
dev.off()
png(paste0(OutDir, "BLUE_UMAPs_type_ident_k20.png"), width = 3000, height = 1000, res = 200)
FeaturePlot(sc.obj.integration, features = "blue_type_k20_UCell", min.cutoff = 'q01', max.cutoff = 'q99',
            pt.size = 1, cols = c(alpha("grey", .5), "red"), split.by = "orig.ident") & NoAxes() &coord_fixed(0.6)
dev.off()
VlnPlot(sc.obj.integration, features = "blue_type_k20_UCell", group.by = "type")

pdf(paste0(OutDir, "PctperCluster_ident.pdf"), width = 15, height = 5)
sc.obj.integration@meta.data %>%
  group_by(orig.ident) %>%
  count(type) %>%
  mutate(rel = round(n/sum(n)*100,2)) %>%
  ungroup() %>% group_by(type) %>%
  mutate(rel2 = round(rel/sum(rel)*100,2)) %>%
  ggplot(aes(x = type, y = rel2, fill = orig.ident)) + geom_col() +
  scale_fill_manual(values = brewer.pal(9,"Set1")[c(1,2,9)]) +
  theme_pubr(x.text.angle = 45) + ylab("Percentage of Cluster") + coord_fixed(0.05)
dev.off()

pdf(paste0(OutDir, "PctperCluster_stress.pdf"), width = 15, height = 5)
sc.obj.integration@meta.data %>%
  group_by(type) %>%
  count(is.Stressed.RGC.Hypoxia) %>%
  mutate(rel = round(n/sum(n)*100,2)) %>%
  ggplot(aes(x = type, y = rel, fill = is.Stressed.RGC.Hypoxia)) + geom_col() +
  scale_fill_manual(values = c(alpha("grey", 1), "dark red")) +
  theme_pubr(x.text.angle = 45) + ylab("Percentage of Cluster")+ coord_fixed(0.05)
dev.off()

# Plot enrichments -------------------------------------------------------------

library(clusterProfiler)
library(enrichplot)

blue_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='blue') %>%
                             dplyr::arrange(desc(kMEblue)) %>% 
                             # dplyr::slice(1:50) %>%
                             magrittr::use_series(GeneSymbol),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")


pdf(paste0(OutDir,"TreePlot_blue.pdf"), width = 10, height = 10)
treeplot(pairwise_termsim(blue_module), showCategory = 20)
dev.off()
pdf(paste0(OutDir,"emapPlot_blue.pdf"), width = 10, height = 10)
emapplot(pairwise_termsim(blue_module), showCategory = 20)
dev.off()

# Use Hypoxia to re-calculate Stress -------------------------------------------

# hypoxia		GO:0071456
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = "GO:0071456", save.UMAP = F, 
                                          new_GO_term_computation = T, 
                                          clustering = "integrated_snn_res.18.reassigned", plot.each.gene = F)

FeaturePlot(sc.obj.integration, features = "Score.GO.0071456", min.cutoff = 'q05', max.cutoff = 'q95',
            pt.size = 1, cols = c(alpha("grey", .5), "red"), split.by = "orig.ident") & NoAxes() & coord_fixed(0.6)



sc.obj.integration$is.Stressed.Merge.RGC <- sc.obj.integration$is.Stressed
sc.obj.integration <- Shiny.GO.thresh(obj = sc.obj.integration,
                                      stress.ident1 = "integrated_snn_res.18.reassigned_cl.av_GO:0071456",
                                      stress.ident2 = "integrated_snn_res.18.reassigned_cl.av_GO:0034976",
                                      notstress.ident3 = "integrated_snn_res.18.reassigned_cl.av_GO:0042063",
                                      notstress.ident4 = "integrated_snn_res.18.reassigned_cl.av_RGC1",
                                      plot.cluster.shiny = "orig.ident")


pl1 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'is.Stressed', legend = T, label = F,
                            save.plot = F, axes = F, cols = c(alpha("grey", 0.5), "dark red"), splitby = "orig.ident")
pl2 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'type', label =F, legend = T,
                            save.plot = F, axes = F, cols = colors_ret, splitby = "orig.ident")
png(paste0(OutDir, "Summary_UMAP_int_stress_hypoxia.png"),width = 5000, height = 2000, res = 300)
cowplot::plot_grid(pl2,pl1, ncol =1 )
dev.off()
sc.obj.integration$is.Stressed.RGC.Hypoxia <- sc.obj.integration$is.Stressed

png(paste0(OutDir, "Stress_RGC_Hypoxia_perident.png"), width = 3000, height = 1000, res = 200 )
DimPlot(sc.obj.integration, group.by = c("is.Stressed.RGC.Hypoxia"),
        split.by = "orig.ident", cols = c(alpha("grey", 0.5), "dark red")) & NoAxes() &coord_fixed(.6)
dev.off()

pdf(paste0(OutDir,"Pct_Stressed.pdf"), width = 5, height = 10)
sc.obj.integration@meta.data %>%
  group_by(orig.ident) %>%
  count(is.Stressed.RGC.Hypoxia) %>%
  mutate(rel = round(n/sum(n)*100,2), source = "Integration") %>%
  ggplot(aes(x = orig.ident, y = rel, fill = is.Stressed.RGC.Hypoxia)) +
  geom_col() + coord_fixed(0.1) + scale_fill_manual(values = c(alpha("grey", 1), "dark red")) +
  theme_pubr(x.text.angle = 45, ) + xlab("") + ylab("Percentage")
dev.off()


plot_df_int_all <- sc.obj.integration@meta.data %>%
  group_by(orig.ident, integrated_snn_res.18.reassigned) %>%
  summarize(across(.cols = c("Score.GO.0071456", 
                             "Score.GO.0034976", 
                             "Score.GO.0042063"), 
                   .fns = function(x) x=mean(x)),
            is.Stressed = unique(is.Stressed),
            number = length(orig.ident))

pdf(paste0(OutDir,"Scores_perGranule_Sridhar.pdf"), width = 10, height = 10)
ggplot(plot_df_int_all,aes(x = Score.GO.0071456,
                           y = Score.GO.0034976, 
                           color = is.Stressed, size = number)) + 
  geom_point(alpha = 0.5) + theme_minimal() + coord_fixed(1) + scale_color_manual(values = c("grey", "dark red")) +
  facet_grid(cols = vars(orig.ident)) +
  xlab("GO:0071456 Score") + ylab("GO:0034976 Score") + 
 # scale_size(name = "GO:0042063 score", limits = c(-0.1,0.15), range = c(1,5)) +
  xlim(-0.1,0.2) + ylim(-0.1,0.2)
dev.off()

# Run Progeny on Integration ---------------------------------------------------

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 

sc.obj.integration <- progeny(sc.obj.integration, scale=FALSE, organism="Human", top=500, perm=1, 
                              return_assay = TRUE)
sc.obj.integration <- Seurat::ScaleData(sc.obj.integration, assay = "progeny")

FeaturePlot(sc.obj.integration, features = "Hypoxia", min.cutoff = 'q05', max.cutoff = 'q95', 
            split.by = "orig.ident",
            cols = c(alpha("grey", 0.5), "red"), ncol = 1, combine = T) & NoAxes() & coord_fixed(0.6)

sc.obj.integration@meta.data %>%
  mutate(bc = row.names(.)) %>%
  left_join(sc.obj.integration@reductions$umap@cell.embeddings %>%
              as.data.frame() %>%
              mutate(bc = row.names(.))) %>%
  left_join(data.frame(Hypoxia = GetAssayData(sc.obj.integration, assay = "progeny")["Hypoxia",])%>%
              mutate(bc = row.names(.))) %>%
  mutate(Hypoxia = clip.outliers(Hypoxia)) %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2, color = type)) +
  geom_point(size = 0.5, alpha = 0.5) + facet_grid(rows = vars(is.Stressed), cols = vars(orig.ident)) +
  scale_color_gradient(low = alpha("grey", 0.5), high = "red") + ggpubr::theme_pubr()

CellsClusters <- data.frame(Cell = row.names(sc.obj.integration@meta.data), 
                            CellType = ifelse(sc.obj.integration$is.Stressed.RGC.Hypoxia == TRUE, "Stressed",
                                              sc.obj.integration$type),
                            orig.ident = sc.obj.integration$orig.ident,
                            stringsAsFactors = FALSE)

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sc.obj.integration, slot = "scale.data", 
                               assay = "progeny"))) %>%
  tibble::rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  mutate(CellType.orig.ident = paste(CellType, orig.ident, sep = "_")) %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


pdf(paste0(OutDir,"progeny_heatmap_filtering.pdf"), width = 10, height = 10)
ggplot(CellsClusters %>% 
         count(CellType), aes(x= CellType, y = n, fill = n>quantile(n, probs =c(.02)))) + 
  geom_col() + coord_fixed(0.001) + scale_fill_manual("> 2nd percentile", values = c(alpha("dark blue", 1), "#EEC61F")) +
  geom_hline(yintercept = quantile(table(CellsClusters$CellType), probs =c(.02))) +
  theme_pubr(x.text.angle = 45, ) + xlab("") + ylab("Percentage")
dev.off()
CellType_pass <- CellsClusters %>% 
  count(CellType) %>%
  filter(n>quantile(n, probs =c(.02))) %>%
  magrittr::use_series(CellType)

paletteLength = 100
myColor = colorRampPalette(c("white", "black"))(paletteLength)


col_labs <- colnames(summarized_progeny_scores_df)[c(4,2:3,5:14)]
pdf(paste0(OutDir,"progeny_heatmap.pdf"), width = 10, height = 10)
pheatmap(t(summarized_progeny_scores_df[CellType_pass,col_labs]),fontsize=14, 
         fontsize_row = 10, cluster_rows = F,
         color=myColor, cellwidth = 20, cellheight = 20,
         scale = "none",
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)
dev.off()

png(paste0(OutDir, "Type_UMAP_overview.png"), width = 2000, height = 1000, res = 200)
DimPlot(sc.obj.integration, group.by = "type", cols = colors_ret) + 
  NoAxes() + coord_fixed(0.6)
dev.off()
png(paste0(OutDir, "Ident_UMAP_overview.png"), width = 2000, height = 1000, res = 200)
DimPlot(sc.obj.integration, group.by = "orig.ident", cols = alpha(brewer.pal(9,"Set1")[c(1,2,9)], 0.5)) + 
  NoAxes() + coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "Scores_UMAP_overview.png"), width = 2000, height = 2000, res = 200)
FeaturePlot(sc.obj.integration, features = c(scores, "RGC1"), cols = c(alpha("grey", .5), "red"),
            max.cutoff = 'q90', min.cutoff = 'q10', ncol = 2) &
  NoAxes() & coord_fixed(0.6)
dev.off()
png(paste0(OutDir, "Hypoxia_Score_UMAP_overview.png"), width = 1000, height = 1000, res = 200)
FeaturePlot(sc.obj.integration, features = c("Score.GO.0071456"), cols = c(alpha("grey", .5), "red"),
            max.cutoff = 'q95', min.cutoff = 'q05', ncol = 1) &
  NoAxes() & coord_fixed(0.6)
dev.off()


saveRDS(sc.obj.integration, paste0(OutDir,'processed_Sridhar_integration_scores_220521.Rds'))

# Filter out and re-integrate + re-cluster!!! ==================================
sc.obj.integration_PRE <- sc.obj.integration
cellIDs.keep <- which_names(!sc.obj.integration_PRE$is.Stressed.RGC.Hypoxia)
subset.obj <- subset(x = sc.obj.integration_PRE, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj, save.plot = F)

subset.obj@active.assay <- "RNA"
subset.obj.split <- SplitObject(subset.obj, split.by = "orig.ident")

object.list_int_anchors <- FindIntegrationAnchors(object.list = subset.obj.split, 
                                                  dims = 1:50)


sc.obj.integration <- IntegrateData(anchorset = object.list_int_anchors, 
                                    dims = 1:50)

sc.obj.integration <- ScaleData(sc.obj.integration, verbose = FALSE)
sc.obj.integration <- RunPCA(sc.obj.integration, npcs = 50, features = VariableFeatures(sc.obj.integration))

sc.obj.integration <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = sc.obj.integration, nPCs = 50, dimensions=3:2, reduction="umap")

sc.obj.integration <- FindNeighbors(sc.obj.integration,dims = 1:50)
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = c(.1,.2,.3,.5))
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = c(.25))

png(paste0(OutDir, "PostGruffi_UMAP_type.png"), width = 1000, height = 1000, res = 200)
DimPlot(sc.obj.integration, group.by = "type", cols = colors_ret) + NoAxes() + coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "PostGruffi_UMAP_type_lab.png"), width = 1000, height = 1000, res = 200)
DimPlot(sc.obj.integration, group.by = "type", cols = colors_ret,
        label = T) + NoAxes() + coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "PostGruffi_UMAP_ident.png"), width = 1000, height = 1000, res = 200)
DimPlot(sc.obj.integration, group.by = "orig.ident", cols = alpha(brewer.pal(9,"Set1")[c(1,2,9)], 0.5)) +
  NoAxes()+ coord_fixed(0.6)
dev.off()

subset.obj@active.assay <- "integrated"
setdiff(VariableFeatures(subset.obj),VariableFeatures(sc.obj.integration))

setdiff(row.names(sc.obj.integration@assays$integrated), row.names(subset.obj@assays$integrated))

saveRDS(sc.obj.integration, file = paste0(OutDir, "Sridhar_2020_postGRUFFIintegration.Rds"))

dev.off()


