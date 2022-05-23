# Sridhar et al. 2020 Retina organoids check for Stress
# Oliver Eichmueller
# Thu May 19 14:59:49 2022 ------------------------------

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer);
require(qgraph); require(enrichplot)

OutDir <- 'Gruffi_revision/Sridhar2020/'
dir.create(OutDir)

# load matrices ----------------------------------------------------------------
d205        <- Seurat::Read10X('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Sridhar_2020/GSE142526_RAW/ORG_D205_filtered_feature_bc_matrix/')
d125_fetalP  <- Seurat::Read10X('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Sridhar_2020/GSE142526_RAW/D125Pfetal_filtered_gene_bc_matrices/GRCh38/')
d125_fetalC  <- Seurat::Read10X('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Sridhar_2020/GSE142526_RAW/D125Cfetal_filtered_gene_bc_matrices/GRCh38/')
meta.data <- read.csv('Gruffi_revision/Sridhar2020/paper metadata/F6/cca_fetalvsorg_125CP_205_metadata.csv', row.names = 1)

meta.data.d205 <- meta.data %>% filter(orig.ident == "D205")
row.names(meta.data.d205) <- gsub("_", "-", row.names(meta.data.d205))
meta.data.d125P <- meta.data %>% filter(orig.ident == "D125P")
row.names(meta.data.d125P) <- gsub("_2", "-1", row.names(meta.data.d125P))
meta.data.d125C <- meta.data %>% filter(orig.ident == "D125C")
row.names(meta.data.d125C) <- gsub("_3", "-1", row.names(meta.data.d125C))

sc.obj.d205         <- CreateSeuratObject(d205, min.cells = 10, min.features = 500, project = "Sridhar_d205", 
                                          meta.data = meta.data.d205)
sc.obj.d125_fetalP  <- CreateSeuratObject(d125_fetalP, min.cells = 10, min.features = 500, project = "Sridhar_d125_fetalP",
                                          meta.data = meta.data.d125P)
sc.obj.d125_fetalC  <- CreateSeuratObject(d125_fetalC, min.cells = 10, min.features = 500, project = "Sridhar_d125_fetalC",
                                          meta.data = meta.data.d125C)
obj.list <- list(sc.obj.d205, sc.obj.d125_fetalP, sc.obj.d125_fetalC)
# percentageo of mito reads

for (i in 1:3) {
  obj.list[[i]][["percent.mito"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
  obj.list[[i]][["percent.ribo"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^RP[SL]")
}

# Preprocess -------------------------------------------------------------------
for (i in 1:3) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], 
                                 normalization.method = "LogNormalize", scale.factor = 10000)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
  obj.list[[i]] = ScaleData(obj.list[[i]], verbose = T)
  
  obj.list[[i]] = RunPCA(obj.list[[i]], npcs = 50, verbose = T)

  obj.list[[i]] = SetupReductionsNtoKdimensions(obj = obj.list[[i]], nPCs = 50, dimensions=3:2, reduction="umap")
  
  obj.list[[i]] = FindNeighbors(obj.list[[i]], reduction = "pca", dims = 1:50)
  
  obj.list[[i]] = FindClusters(obj.list[[i]], resolution = c(0.1,0.2,0.3,0.5))
}

# Plot -------------------------------------------------------------------------
colors_ret <- colorRampPalette(c(pals::tol(12), "grey"))(14)
names(colors_ret) <- meta.data$type %>% unique
png(paste0(OutDir, "UMAP_type.png"), width = 5000, height = 2000, res = 300)
cowplot::plot_grid(
  clUMAP(ident = "type", obj = obj.list[[1]], save.plot = F, title = "D205", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)),
  clUMAP(ident = "type", obj = obj.list[[2]], save.plot = F, title = "D125P", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)),
  clUMAP(ident = "type", obj = obj.list[[3]], save.plot = F, title = "D125C", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)), ncol = 3)
dev.off()
pdf(paste0(OutDir, "UMAP_type.pdf"), width = 15, height = 7.5)
cowplot::plot_grid(
  clUMAP(ident = "type", obj = obj.list[[1]], save.plot = F, title = "D205", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)),
  clUMAP(ident = "type", obj = obj.list[[2]], save.plot = F, title = "D125P", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)),
  clUMAP(ident = "type", obj = obj.list[[3]], save.plot = F, title = "D125C", cols = colors_ret, legend = F, axes = F) +
    theme(title = element_text(size = 30)), ncol = 3)
dev.off()

# Run Gruffi -------------------------------------------------------------------

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

for (i in 1:3) {
  obj.list[[i]] <- aut.res.clustering(obj = obj.list[[i]])
  
  granule.res.4.gruffi <- obj.list[[i]]@misc$gruffi$'optimal.granule.res'	
  
  obj.list[[i]] <- reassign.small.clusters(obj.list[[i]], ident = granule.res.4.gruffi)
  
  granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')
  
  
  
  # Glycolytic process	GO:0006096
  obj.list[[i]] <- GO_score_evaluation(obj = obj.list[[i]], GO_term = go1, save.UMAP = F, 
                                       new_GO_term_computation = T, 
                                       clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  # ER stress 	GO:0034976
  obj.list[[i]] <- GO_score_evaluation(obj = obj.list[[i]], GO_term = go2, save.UMAP = F, 
                                       new_GO_term_computation = T, 
                                       clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  # Gliogenesis		GO:0042063
  obj.list[[i]] <- GO_score_evaluation(obj = obj.list[[i]], GO_term = go3, save.UMAP = F, 
                                       new_GO_term_computation = T, 
                                       clustering = granule.res.4.gruffi, plot.each.gene = F)
}

# Shiny for D205
granule.res.4.gruffi <- paste0(obj.list[[1]]@misc$gruffi$optimal.granule.res, '.reassigned')
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
obj.list[[1]] <- Shiny.GO.thresh(obj = obj.list[[1]],
                                 stress.ident1 = i1,
                                 stress.ident2 = i2,
                                 notstress.ident3 = i3,
                                 plot.cluster.shiny = "orig.ident")

# Shiny for D125P
granule.res.4.gruffi <- paste0(obj.list[[2]]@misc$gruffi$optimal.granule.res, '.reassigned')
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
obj.list[[2]] <- Shiny.GO.thresh(obj = obj.list[[2]],
                                 stress.ident1 = i1,
                                 stress.ident2 = i2,
                                 notstress.ident3 = i3,
                                 plot.cluster.shiny = "orig.ident")

# Shiny for D125C
granule.res.4.gruffi <- paste0(obj.list[[3]]@misc$gruffi$optimal.granule.res, '.reassigned')
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
obj.list[[3]] <- Shiny.GO.thresh(obj = obj.list[[3]],
                                 stress.ident1 = i1,
                                 stress.ident2 = i2,
                                 notstress.ident3 = i3,
                                 plot.cluster.shiny = "orig.ident")


pl1 <- Seurat.utils::clUMAP(obj = obj.list[[1]], 'is.Stressed', label =F, 
                            save.plot = F, title = "D205", axes = F, cols = c(alpha("grey", 0.5), "dark red"))
pl2 <- Seurat.utils::clUMAP(obj = obj.list[[2]], 'is.Stressed', label =F, 
                            save.plot = F, title = "D125P", axes = F, cols = c(alpha("grey", 0.5), "dark red"))
pl3 <- Seurat.utils::clUMAP(obj = obj.list[[3]], 'is.Stressed', label =F, 
                            save.plot = F, title = "D125C", axes = F, cols = c(alpha("grey", 0.5), "dark red"))
png(paste0(OutDir, "UMAP_Stress.png"), width = 5000, height = 2000, res = 300)
cowplot::plot_grid(pl1,pl2, pl3, ncol = 3)
dev.off()
pdf(paste0(OutDir, "UMAP_stress.pdf"), width = 15, height = 7.5)
cowplot::plot_grid(pl1,pl2, pl3, ncol = 3)
dev.off()

scores <- c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063")

for (i in 1:3) {
   pl <- FeaturePlot(obj.list[[i]], features = scores, min.cutoff = 'q10', max.cutoff = 'q90', 
              cols = c(alpha("grey", 0.5), "red"), ncol = 3, combine = T) & NoAxes() & coord_fixed(0.6)
  ggsave(paste0(OutDir, unique(obj.list[[i]]$orig.ident), "_Scores.png"), pl, width = 15, height = 7.5)
}

FeaturePlot(obj.list[[1]], features = scores[1], min.cutoff = 'q10', max.cutoff = 'q90', 
            cols = c(alpha("grey", 0.5), "red"), ncol = 1, combine = T) + NoAxes() + coord_fixed(0.6) +
  scale_color_gradient(low = alpha("grey", 0.5),
                       high = "red", limits = round(quantile(obj.list[[1]][[scores[1]]][,1], probs = c(.1,.9)),2))
hist(obj.list[[1]][[scores[3]]][,1])

FeaturePlot(obj.list[[2]], features ="Hypoxia", min.cutoff = 'q05', max.cutoff = 'q95', 
            cols = c(alpha("grey", 0.5), "red"), ncol = 1, combine = T) + NoAxes() + coord_fixed(0.6)

FeaturePlot(obj.list[[1]], features =c("percent.mito", "percent.ribo", "DDIT4"), min.cutoff = 'q05', max.cutoff = 'q95', 
            cols = c(alpha("grey", 0.5), "red"), ncol = 1, combine = T) & NoAxes() & coord_fixed(0.6)

# Run Progeny ------------------------------------------------------------------

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
for (i in 1:3) {
  obj.list[[i]] <- progeny(obj.list[[i]], scale=FALSE, organism="Human", top=500, perm=1, 
                    return_assay = TRUE)
  obj.list[[i]] <- Seurat::ScaleData(obj.list[[i]], assay = "progeny")
}

for (i in 1:3) {
  pl <- FeaturePlot(obj.list[[i]], features = "Hypoxia", min.cutoff = 'q10', max.cutoff = 'q90', 
                    cols = c(alpha("grey", 0.5), "red"), ncol = 1, combine = T) + NoAxes() + coord_fixed(0.6) +
    scale_color_gradient(low = alpha("grey", 0.5),
                         high = "red", limits = quantile(obj.list[[1]]@assays$progeny@data["Hypoxia",], probs = c(0.0,0.95)))
  ggsave(paste0(OutDir, unique(obj.list[[i]]$orig.ident), "_Hypoxia_Score.png"), pl, width = 5, height = 3)
}



CellsClusters <- data.frame(Cell = row.names(sc.obj@meta.data), 
                            CellType = ifelse(sc.obj$is.Stressed == TRUE, "Stressed",
                                              sc.obj$RNA_snn_res.0.1),
                            stringsAsFactors = FALSE)

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sc.obj, slot = "scale.data", 
                               assay = "progeny"))) %>%
  tibble::rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

# Integrate the data -----------------------------------------------------------

obj.list_int_anchors <- 
  FindIntegrationAnchors(object.list = obj.list, dims = 1:50, assay = rep("RNA", 3))

sc.obj.integration <- IntegrateData(anchorset = obj.list_int_anchors, 
                                    dims = 1:50)

sc.obj.integration = ScaleData(sc.obj.integration, verbose = T)

sc.obj.integration = RunPCA(sc.obj.integration, npcs = 50, verbose = T)

sc.obj.integration = SetupReductionsNtoKdimensions(obj = sc.obj.integration, nPCs = 50, dimensions=3:2, reduction="umap")

sc.obj.integration = FindNeighbors(sc.obj.integration, reduction = "pca", dims = 1:50)

sc.obj.integration = FindClusters(sc.obj.integration, resolution = c(0.1,0.2,0.3,0.5))

clUMAP(ident = "type", obj = sc.obj.integration, save.plot = F, 
       #title = "D205", 
       cols = colors_ret, 
       legend = F, axes = F) +
  theme(title = element_text(size = 30))

sc.obj.integration$is.Stressed.Single <- sc.obj.integration$is.Stressed


# Run Gruffi again
ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

sc.obj.integration <- aut.res.clustering(obj = sc.obj.integration)

granule.res.4.gruffi <- sc.obj.integration@misc$gruffi$'optimal.granule.res'	

sc.obj.integration <- reassign.small.clusters(sc.obj.integration, ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')



# Glycolytic process	GO:0006096
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = go1, save.UMAP = F, 
                                     new_GO_term_computation = T, 
                                     clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
sc.obj.integration <- GO_score_evaluation(obj =sc.obj.integration, GO_term = go2, save.UMAP = F, 
                                     new_GO_term_computation = T, 
                                     clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = go3, save.UMAP = F, 
                                     new_GO_term_computation = T, 
                                     clustering = granule.res.4.gruffi, plot.each.gene = F)

# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
sc.obj.integration <- Shiny.GO.thresh(obj = sc.obj.integration,
                                 stress.ident1 = i1,
                                 stress.ident2 = i2,
                                 notstress.ident3 = i3,
                                 plot.cluster.shiny = "orig.ident")

# Add Score for RGC
sc.obj.integration <- AddModuleScore(sc.obj.integration, assay = "RNA", 
                                     features = list(RGC = c("ISL1", "SNCG", "RBPMS")),name = "RGC")
sc.obj.integration@misc$GO$RGC <- c("ISL1", "SNCG", "RBPMS")
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = "RGC1", save.UMAP = F, 
                                          new_GO_term_computation = F, 
                                          clustering = granule.res.4.gruffi, plot.each.gene = F)

sc.obj.integration$integrated_snn_res.18.reassigned_cl.av_RGC1.vec <- 
  as.numeric(as.vector(sc.obj.integration$integrated_snn_res.18.reassigned_cl.av_RGC1))

sc.obj.integration$is.Stressed.Merge <- sc.obj.integration$is.Stressed
sc.obj.integration <- Shiny.GO.thresh(obj = sc.obj.integration,
                                      stress.ident1 = i1,
                                      stress.ident2 = i2,
                                      notstress.ident3 = i3,
                                      notstress.ident4 = "integrated_snn_res.18.reassigned_cl.av_RGC1",
                                      plot.cluster.shiny = "orig.ident")


pl1 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'is.Stressed', legend = T, label = F,
                     save.plot = F, axes = F, cols = c(alpha("grey", 0.5), "dark red"), splitby = "orig.ident")
pl2 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'type', label =F, legend = T,
                     save.plot = F, axes = F, cols = colors_ret, splitby = "orig.ident")
png(paste0(OutDir, "Summary_UMAP_int_stress.png"),width = 5000, height = 2000, res = 300)
cowplot::plot_grid(pl2,pl1, ncol =1 )
dev.off()
png(paste0(OutDir, "Summary_UMAP_int_stress_nolab.png"),width = 5000, height = 2000, res = 300)
cowplot::plot_grid(pl2,pl1, ncol =1 )
dev.off()

png(paste0(OutDir, "Scores_UMAP_int_stress.png"),width = 3500, height = 3500, res = 300)
FeaturePlot(sc.obj.integration, features = c(scores, "RGC1"), min.cutoff = 'q05', max.cutoff = 'q95', split.by = "orig.ident",
            cols = c(alpha("grey", 0.5), "red"), combine = T) & NoAxes() & coord_fixed(0.6)
dev.off()

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
  ggplot(aes(x=UMAP_1,y=UMAP_2, color = Hypoxia)) +
  geom_point(size = 0.5, alpha = 0.5) + facet_grid(rows = vars(is.Stressed), cols = vars(orig.ident)) +
  scale_color_gradient(low = alpha("grey", 0.5), high = "red") + ggpubr::theme_pubr()

CellsClusters <- data.frame(Cell = row.names(sc.obj.integration@meta.data), 
                            CellType = ifelse(sc.obj.integration$is.Stressed == TRUE, "Stressed",
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
  group_by(Pathway, CellType.orig.ident) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
         fontsize_row = 10, 
         color=myColor, 
         scale = "row",
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)
dev.off()

progeny_scores_df %>%
  filter(Pathway == "Hypoxia") %>%
  ggplot(aes(x=CellType, y = Activity, fill = CellType)) + 
  geom_violin() + facet_grid(rows = vars(orig.ident))


saveRDS(sc.obj.integration, file = paste0(OutDir, "d205_d125P_d125C_integration_progeny.Rds"))
