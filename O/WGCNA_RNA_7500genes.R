## scWGCNA on Stress Project #2
## Fri Nov 26 12:42:19 2021 ------------------------------


## Oliver Eichmueller
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
library(conflicted)

OutDir <- 'WGCNA_test_7.5k_variable/'
dir.create(OutDir)
gc()

metacell_seurat <- readRDS('WGCNA_test_RNA/metacell_seurat_ds.Rds')

library(CodeAndRoll2)
library(Seurat.utils)

metacell_seurat <- FindVariableFeatures(metacell_seurat, nfeatures = 7500)

metacell_seurat <- metacell_seurat[VariableFeatures(metacell_seurat),]

# Run WGCNA
gc()
disableWGCNAThreads()

# how many groups are there
nclusters <- length(unique(metacell_seurat$integrated_snn_res.0.3.Manual.short))

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
);

# Plot the results:
pdf("WGCNA_test_7.5k_variable/1_Power.pdf", height=10, width=18)

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

enableWGCNAThreads(4)
softPower=7

nSets = 1
setLabels = 'Meta'
shortLabels = setLabels

multiExpr <- list()

# seems like datExpr needs to be a data.frame, maybe ask why
multiExpr[['Meta']] <- list(data=as.data.frame(datExpr))

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
write.csv(geneInfo,file=paste0('WGCNA_test_7.5k_variable/geneInfoSigned_ds_7.5kgenes_res03.csv'))

PCvalues=MEs

# Visualisation
pdf("WGCNA_test_7.5k_variable/SignedDendro.pdf",height=5, width=8)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("Gene dendrogram and module colors"))
dev.off()

plot_df <- cbind(dplyr::select(targets, c(integrated_snn_res.0.3.Manual.short)), PCvalues)
plot_df <- reshape2::melt(plot_df, id.vars = c('integrated_snn_res.0.3.Manual.short'))
plot_df$integrated_snn_res.0.3.Manual.short <- factor(plot_df$integrated_snn_res.0.3.Manual.short)

colors <- sub('ME', '', as.character(levels(plot_df$variable)))
p <- ggplot(plot_df, aes(x=integrated_snn_res.0.3.Manual.short, y=value, fill=integrated_snn_res.0.3.Manual.short)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('Module Eigengene') + xlab('') +
  facet_grid(~variable, scales='free_y') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )

# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf('WGCNA_test_7.5k_variable/ME_trajectory_Plot_modulevsCondition.pdf',width=w,height=8,useDingbats=F)
p 
dev.off()

p <- ggplot(plot_df, aes(x=integrated_snn_res.0.3.Manual.short, y=value, 
                         fill=integrated_snn_res.0.3.Manual.short)) +
  geom_boxplot(notch=FALSE) +
  RotatedAxis() + ylab('Module Eigengene') + xlab('') +
  facet_grid(rows = vars(variable), scales='free') +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
  )
# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf('WGCNA_test_7.5k_variable/MCells_ME_trajectory_Plot_modulevsres03.pdf',width=w,height=8,useDingbats=F)
p 
dev.off()


write.csv(plot_df, file = 'WGCNA_test_7.5k_variable/plot_df.csv')

# Network Plots

load("ConsensusTOM-block.1.rda")

TOM.matrix = as.matrix(datExpr);
uniquemodcolors = unique(geneInfo$Initially.Assigned.Module.Color);
uniquemodcolors=uniquemodcolors

pdf(paste0('WGCNA_test_7.5k_variable/ModuleNetworks_res03.pdf'),height=9,width=10);

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
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
png('WGCNA_test_7.5k_variable/NetworkHeatmapPlot_ds_7.5k.png')
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, 7500 variable genes", col = myheatcol)
dev.off()

saveRDS(datExpr, 'WGCNA_test_7.5k_variable/datExpr_WGCNA_res03_211124.Rds')
saveRDS(net, 'WGCNA_test_7.5k_variable/net_WGCNA_res03_211124.Rds')
#saveRDS(metacell_seurat, 'WGCNA_test_7.5k_variable/metacell_seurat_WGCNA_res03_211124.Rds')
write.csv(geneInfo,file='WGCNA_test_7.5k_variable/geneInfo_signed_WGCNA_res03_211124.Rds')


# Use modules for plotting in seurat -------------------------------------------

sc.obj <- 
  readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/ReAnalyze.alternatives.R_2021_11_23-13h/combined.obj_before.Stress.filtering_2021.11.23_13.29.Rds.gz')


gene_mod_split = list()

for (i in 1:length(uniquemodcolors)) {
  color_oi = uniquemodcolors[i]
  gene_mod_split[[color_oi]] <- geneInfo %>% filter(Initially.Assigned.Module.Color == color_oi) %>%
    dplyr::arrange(dplyr::desc(paste0("kME", color_oi))) %>%
    magrittr::use_series(GeneSymbol)
}


names(gene_mod_split) <- paste0(names(gene_mod_split), "_res03_7.5k")

sc.obj <- AddModuleScore(sc.obj, features = gene_mod_split, assay = "RNA")

colnames(sc.obj@meta.data)[169:181] <- names(gene_mod_split)

pl1 <- FeaturePlot(sc.obj, features = "magenta_res03_7.5k", max.cutoff = 'q99', min.cutoff = 'q01',
                   order = F
                   , pt.size = 0.1, cols = c(alpha("grey", .5), "red"))



VlnPlot(sc.obj, features = "magenta_res03_7.5k", group.by = "integrated_snn_res.0.3.Manual.short", pt.size = 0)

png('WGCNA_test_7.5k_variable/Modules_res03_WGCNA.png', width = 2000, height = 1500)
FeaturePlot(sc.obj, features = names(gene_mod_split),
            max.cutoff = 'q99', min.cutoff = 'q01',
            order = F
            , pt.size = 0.1
            , cols = c(alpha("grey", .5), "red"))
dev.off()

metacell_seurat <- AddModuleScore(metacell_seurat, features = gene_mod_split, assay = "RNA")

colnames(metacell_seurat@meta.data)[5:17] <- names(gene_mod_split)

png('WGCNA_test_7.5k_variable/Modules_metacells_res03_WGCNA.png', width = 1000, height = 1000)
FeaturePlot(metacell_seurat, features = names(gene_mod_split),
            ncol = 3, pt.size = .1, cols = c(alpha("grey", .5), "red"))
dev.off()


saveRDS(metacell_seurat, 'WGCNA_test_7.5k_variable/metacell_seurat_res03_5k_scores.rds')

# Perform String analysis

library(STRINGdb)
string_db <- STRINGdb$new(species = 9606)

example_string <- string_db$map(geneInfo %>% filter(Initially.Assigned.Module.Color=='pink'),
                                "GeneSymbol")
string_hits <- example_string %>% arrange(desc(kMEpink)) %>% magrittr::use_series(STRING_id)


string_db$plot_network(string_hits[!is.na(string_hits)])
enrichment <- string_db$get_enrichment(string_hits[!is.na(string_hits)])

enrichment$category %>% table
enrichment %>%
  dplyr::select(category, term, description, p_value, fdr, number_of_genes) %>%
  filter(category %in% "Process") %>%
  arrange((p_value)) %>% dplyr::slice(1:20) %>%
  grid.table()
dev.off()

enrichment %>% 
  dplyr::select(category, term, description, p_value, fdr, number_of_genes) %>%
  filter(term %in% "GO.0034976")

library(clusterProfiler)
library(enrichplot)

magenta_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='magenta') %>%
                          dplyr::arrange(desc(kMEmagenta)) %>% 
                          # dplyr::slice(1:50) %>%
                          magrittr::use_series(GeneSymbol),
                        OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(magenta_module, showCategory = 30)+
  scale_color_viridis_c(direction = -1)
  ggplot2::labs(title = "Stress Module") +
  theme(plot.title = element_text(hjust= 0.5))

ggplot(magenta_module@result %>% 
         mutate(Count_term = stringr::str_split(GeneRatio, "/", simplify = T)[,2]) %>% 
         mutate(pct = round(as.numeric(Count)/as.numeric(Count_term)*100, 2)) %>%
         filter(Count>6) %>%
         arrange(p.adjust) %>% dplyr::slice(1:24) %>%
         mutate(Description = factor(Description, levels = rev(.$Description))), 
       aes(x = pct, y = Description, size = Count, color = -log(p.adjust)))+
  geom_point() + 
  scale_size(breaks = c(3,5,7,10,15,20), limits = c(3,18)) +
  xlim(1,30) + scale_color_gradient() + theme_linedraw()
  
pdf('WGCNA_test_7.5k_variable/cnet_magenta.pdf', width = 20, height = 20)
cnetplot(magenta_module, showCategory = 40)
dev.off()
magenta_module_pw <- pairwise_termsim(magenta_module)

pdf('WGCNA_test_7.5k_variable/tree_magenta.pdf', width = 15, height = 15)
treeplot(magenta_module_pw, showCategory = 40)
dev.off()

pink_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='pink') %>%
                          dplyr::arrange(desc(kMEpink)) %>% 
                          # dplyr::slice(1:50) %>%
                          magrittr::use_series(GeneSymbol),
                        OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(pink_module, showCategory = 20)+
  ggplot2::labs(title = "Pink Module") +
  theme(plot.title = element_text(hjust= 0.5))

pink_module_pw <- pairwise_termsim(pink_module)

treeplot(pink_module_pw, showCategory = 20)

brown_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='brown') %>%
                           magrittr::use_series(GeneSymbol),
                         OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(brown_module, showCategory = 30)+
  ggplot2::labs(title = "Brown Module") +
  theme(plot.title = element_text(hjust= 0.5))

brown_module_pw <- pairwise_termsim(brown_module)

treeplot(brown_module_pw, showCategory = 50)

turquoise_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='turquoise') %>%
                               dplyr::arrange(desc(kMEturquoise)) %>% dplyr::slice(1:50) %>%
                               magrittr::use_series(GeneSymbol),
                             OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(turquoise_module, showCategory = 20)+
  ggplot2::labs(title = "Turquoise Module") +
  theme(plot.title = element_text(hjust= 0.5))

turquoise_module_pw <- pairwise_termsim(turquoise_module)

treeplot(turquoise_module_pw, showCategory = 20)


greenyellow_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='greenyellow') %>%
                            magrittr::use_series(GeneSymbol),
                          OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(greenyellow_module, showCategory = 30)+
  ggplot2::labs(title = "Greenyellow Module") +
  theme(plot.title = element_text(hjust= 0.5))

greenyellow_module_pw <- pairwise_termsim(greenyellow_module)

treeplot(greenyellow_module_pw, showCategory = 20)


median_scores <- sc.obj@meta.data %>%
  dplyr::group_by(integrated_snn_res.0.3.Manual.short) %>%
  summarise(across(.cols = names(gene_mod_split)[-5], .fns = function(x = .x) x = median(x)))

mean_scores <- sc.obj@meta.data %>%
  dplyr::group_by(integrated_snn_res.0.3.Manual.short) %>%
  summarise(across(.cols = names(gene_mod_split)[-5], .fns = function(x = .x) x = mean(x)))

median_scores_df <- data.frame(median_scores[,-1], row.names = median_scores$integrated_snn_res.0.3.Manual.short)
mean_scores_df <- data.frame(mean_scores[,-1], row.names = mean_scores$integrated_snn_res.0.3.Manual.short)


pheatmap::pheatmap(clip.outliers(scale(median_scores_df)))
pheatmap::pheatmap(clip.outliers(scale(mean_scores_df)))


# Test GO.0001666, GO.0036293, GO.0070482
source("https://raw.githubusercontent.com/vertesy/TheCorvinas/0819f413f7fff4b3e446eed44b314022931bb364/R/GO-scoring/Seurat.gene.sets.and.GO.terms.R")
sc.obj@active.assay <- "RNA"
top_GOs <- magenta_module@result %>% 
  filter(Count>6) %>%
  arrange(p.adjust) %>% dplyr::slice(1:24) %>% magrittr::use_series(ID)

GO_add <- top_GOs[!(stringr::str_replace(top_GOs, ":", ".") %in% names(sc.obj@misc$GO))]

for (i in 1:length(GO_add)) {
  sc.obj <- GetGOTerms(sc.obj, GO = GO_add[i], web.open = F)
  sc.obj <- AddGOScore(sc.obj, GO = GO_add[i])
  
}

sc.obj <- GetGOTerms(sc.obj, GO = "GO:0001666", web.open = F)
sc.obj <- AddGOScore(sc.obj, GO = "GO:0001666")

sc.obj <- GetGOTerms(sc.obj, GO = "GO:0036293", web.open = F)
sc.obj <- AddGOScore(sc.obj, GO = "GO:0036293")

sc.obj <- GetGOTerms(sc.obj, GO = "GO:0070482", web.open = F)
sc.obj <- AddGOScore(sc.obj, GO = "GO:0070482")

FeaturePlot(sc.obj, features = c("Score.GO.0071456", "Score.GO.0001666"), min.cutoff = 'q01', max.cutoff = 'q99')
dev.off()
library(gridExtra)

plot_list <- list()  

for (i in 1:(length(top_GOs))) {
  title = magenta_module[top_GOs[i],"Description"]
  plot_list[[i]] <- FeaturePlot(sc.obj, features = paste0("Score.", stringr::str_replace(top_GOs[i], ":", ".")), 
              min.cutoff = 'q01', max.cutoff = 'q99') + ggtitle(title)
  
}

library(patchwork)
library(cowplot)
library(tidyverse)

pdf('WGCNA_test_7.5k_variable/scores_GO_magenta.pdf', width = 15, height = 15,onefile = T)

for (i in seq(1,24, by = 4)) {
  
  print(plot_grid(plotlist = plot_list[i:(i+3)]))
  
}
dev.off()


goi <- names(sc.obj@misc$expr.q90)[sc.obj@misc$expr.q90>0]

sc.obj@misc$GO$GO.0071456_q90 <- intersect(sc.obj@misc$GO$GO.0071456, goi)
sc.obj@misc$GO$GO.0001666_q90 <- intersect(sc.obj@misc$GO$GO.0001666, goi)
sc.obj@misc$GO$GO.0071456_q90_kickout <- setdiff(sc.obj@misc$GO$GO.0071456, goi)
sc.obj@misc$GO$GO.0001666_q90_kickout <- setdiff(sc.obj@misc$GO$GO.0001666, goi)

sc.obj <- AddGOScore(sc.obj, GO = "GO.0071456_q90")
sc.obj <- AddGOScore(sc.obj, GO = "GO.0001666_q90")
sc.obj <- AddGOScore(sc.obj, GO = "GO.0071456_q90_kickout")
sc.obj <- AddGOScore(sc.obj, GO = "GO.0001666_q90_kickout")

FeaturePlot(sc.obj, features = c("Score.GO.0071456", "Score.GO.0001666",
                                 "Score.GO.0071456_q90", "Score.GO.0001666_q90"), 
            min.cutoff = 'q01', max.cutoff = 'q99')

intersect(geneInfo %>% filter(Initially.Assigned.Module.Color == "magenta") %>% magrittr::use_series(GeneSymbol),
          sc.obj@misc$GO$GO.0071456_q90)

saveRDS(sc.obj, 'WGCNA_test_7.5k_variable/scobj_scores_WGCNA.Rds')
saveRDS(metacell_seurat, 'WGCNA_test_7.5k_variable/metacell_seurat_res03_5k_scores.rds')

