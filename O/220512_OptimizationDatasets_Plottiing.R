# Figure for organoid optimization datasets
# Oliver Eichmueller
# Thu May 12 17:01:47 2022 ------------------------------
require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer);
require(qgraph); require(enrichplot); require(readr); require(progeny)

OutDir <- 'Gruffi_revision/Figure_OptimizationDatasets/'
dir.create(OutDir)

# load datasets ----------------------------------------------------------------
sc.obj.C <- readRDS("Gruffi_revision/Cakir_2019/analysis/Cakir_object_proc_220510.Rds")
sc.obj.Q <- readRDS("Gruffi_revision/Qian_2020/analysis/Qian_object_proc_220510.Rds")
sc.obj.G <- readRDS("Gruffi_revision/Giandomenico_2019/analysis/Giandomenico_object_proc_220510.Rds")
obj.list <- list(sc.obj.C, sc.obj.G, sc.obj.Q)
# Plots for Cakir et al --------------------------------------------------------

cowplot::plot_grid(
clUMAP(obj = obj.list[[1]], ident = "predicted.cell.type.N", save.plot = F),
clUMAP(obj = obj.list[[1]], ident = "RNA_snn_res.0.1", save.plot = F))
qUMAP(obj = obj.list[[1]],feature = "DCX", save.plot = F)
label.Cakir <- c("RG 1", "EN 1", "Endothelial 1", "RG 2", "Endothelial 2", "immature EN", "EN 2",
                 "Dividing", "Outlier 1", "Outlier 2", "Outlier 3")
obj.list[[1]]$RNA_snn_res.0.1.named <-
  factor(obj.list[[1]]$RNA_snn_res.0.1, labels = label.Cakir)

png(paste0(OutDir, "Cakir_UMAP_named.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[1]], ident = "RNA_snn_res.0.1.named", save.plot = F, cols = pals::tol(11),
       legend = T, label = F) + NoAxes()
dev.off()

png(paste0(OutDir, "Cakir_UMAP_Scores.png"), width = 750, height = 500)
FeaturePlot(obj.list[[1]], features = c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063", "Score.GO.0001935"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red")) & NoAxes() & coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "Cakir_UMAP_Stress.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[1]], ident = "is.Stressed", save.plot = F, cols = c("grey", "dark red"),
       legend = T, label = F) + NoAxes()
dev.off()

sc.obj.C.post <- readRDS('Gruffi_revision/Cakir_2019/analysis/Cakir_object_POSTGruffi_proc_220510.Rds')

sc.obj.C.post$RNA_snn_res.0.1.named <- obj.list[[1]]$RNA_snn_res.0.1.named[names(sc.obj.C.post$RNA_snn_res.0.1)]

png(paste0(OutDir, "Cakir_UMAP_named_post.png"), width = 750, height = 500)
clUMAP(obj = sc.obj.C.post, ident = "RNA_snn_res.0.1.named", save.plot = F, cols = pals::tol(11),
       legend = T, label = F) + NoAxes()
dev.off()

# Plots for Giandomenico et al -------------------------------------------------

cowplot::plot_grid(
clUMAP(obj = obj.list[[2]], ident = "predicted.cell.type.N", save.plot = F),
clUMAP(obj = obj.list[[2]], ident = "RNA_snn_res.0.3", save.plot = F))
qUMAP(obj = obj.list[[2]],feature = "DCX", save.plot = F)
label.Giandomenico <- c("EN 1", "immature EN", "RG 1", "EN 2", "IN", "EN 3",
                 "Dividing", "RG 2")
obj.list[[2]]$RNA_snn_res.0.1.named <-
  factor(obj.list[[2]]$RNA_snn_res.0.1, labels = label.Giandomenico)

png(paste0(OutDir, "Giandomenico_UMAP_named.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[2]], ident = "RNA_snn_res.0.1.named", save.plot = F, cols = pals::tol(11),
       legend = T, label = F) + NoAxes()
dev.off()

png(paste0(OutDir, "Giandomenico_UMAP_Scores.png"), width = 750, height = 500)
FeaturePlot(obj.list[[2]], features = c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red")) & NoAxes() & coord_fixed(0.6)
dev.off()

FeaturePlot(obj.list[[2]], features = c("Hypoxia"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red"), ncol = 3) & 
  NoAxes() & coord_fixed(0.6)


png(paste0(OutDir, "Giandomenico_UMAP_Stress.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[2]], ident = "is.Stressed", save.plot = F, cols = c("grey", "dark red"),
       legend = T, label = F) + NoAxes()
dev.off()

sc.obj.G.post <- readRDS('Gruffi_revision/Giandomenico_2019/analysis/Giandomenico_object_POSTGruffi_proc_220510.Rds')

sc.obj.G.post$RNA_snn_res.0.1.named <- obj.list[[2]]$RNA_snn_res.0.1.named[names(sc.obj.G.post$RNA_snn_res.0.1)]

png(paste0(OutDir, "Giandomenico_UMAP_named_post.png"), width = 750, height = 500)
clUMAP(obj = sc.obj.G.post, ident = "RNA_snn_res.0.1.named", save.plot = F, cols = pals::tol(11),
       legend = T, label = F) + NoAxes()
dev.off()


# Plots for Qian et al -------------------------------------------------

cowplot::plot_grid(
  clUMAP(obj = obj.list[[3]], ident = "predicted.cell.type.N", save.plot = F),
  clUMAP(obj = obj.list[[3]], ident = "RNA_snn_res.0.1", save.plot = F))
qUMAP(obj = obj.list[[2]],feature = "DCX", save.plot = F)
label.Qian <- c("EN 1", "EN 2", "RG", "Outlier", "EN 3")
obj.list[[3]]$RNA_snn_res.0.1.named <-
  factor(obj.list[[3]]$RNA_snn_res.0.1, labels = label.Qian)

png(paste0(OutDir, "Qian_UMAP_named.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[3]], ident = "RNA_snn_res.0.1.named", save.plot = F, cols = pals::tol(5),
       legend = T, label = F) + NoAxes()
dev.off()

png(paste0(OutDir, "Qian_UMAP_Scores.png"), width = 750, height = 500)
FeaturePlot(obj.list[[3]], features = c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red")) & NoAxes() & coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "Qian_UMAP_Hypoxia.png"), width = 375, height = 500)
FeaturePlot(obj.list[[3]], features = c("Hypoxia"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red"), ncol = 1) & 
  NoAxes() & coord_fixed(0.6)
dev.off()

png(paste0(OutDir, "Qian_UMAP_Stress.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[3]], ident = "is.Stressed", save.plot = F, cols = c("grey", "dark red"),
       legend = T, label = F) + NoAxes()
dev.off()

png(paste0(OutDir, "Qian_UMAP_Stress_Manual.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[3]], ident = "is.Stressed.Manual", save.plot = F, cols = c("grey", "dark red"),
       legend = T, label = F) + NoAxes()
dev.off()

# Try out different granule size -----------------------------------------------

obj.list[[3]] <- aut.res.clustering(obj.list[[3]], lower.median = 75, upper.median = 150, upper.res = 20)


granule.res.4.gruffi <- obj.list[[3]]@misc$gruffi$'optimal.granule.res'	

obj.list[[3]] <- reassign.small.clusters(obj.list[[3]], ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

# Glycolytic process	GO:0006096
obj.list[[3]] <- GO_score_evaluation(obj = obj.list[[3]], GO_term = go1, save.UMAP = F, 
                                     new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
obj.list[[3]] <- GO_score_evaluation(obj = obj.list[[3]], GO_term = go2, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
obj.list[[3]] <- GO_score_evaluation(obj = obj.list[[3]], GO_term = go3, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)


# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))
obj.list[[3]]$is.Stressed.old <- obj.list[[3]]$is.Stressed
# Call Shiny app
obj.list[[3]] <- Shiny.GO.thresh(obj = obj.list[[3]],
                          stress.ident1 = i1,
                          stress.ident2 = i2,
                          notstress.ident3 = i3,
                          plot.cluster.shiny = "orig.ident")


obj.list[[3]]@meta.data %>% count(is.Stressed) %>% mutate(pct = round(n/sum(n)*100,2))


png(paste0(OutDir, "Qian_UMAP_Stress_GranuleOpt.png"), width = 750, height = 500)
clUMAP(obj = obj.list[[3]], ident = "is.Stressed", save.plot = F, cols = c("grey", "dark red"),
       legend = T, label = F) + NoAxes()
dev.off()


# Progeny for all of them ------------------------------------------------------

# Check progeny for stressed cells ---------------------------------------------

summarized_progeny_scores_df.list <- list()
# Run Progeny on all three datasets
for (i in 1:3) {
  
  CellsClusters <- data.frame(Cell = row.names(obj.list[[i]]@meta.data), 
                              CellType = ifelse(obj.list[[i]]$is.Stressed == TRUE, "Stressed",
                                                as.vector(obj.list[[i]]$RNA_snn_res.0.1.named)),
                              stringsAsFactors = FALSE)
  
  
  ## We compute the Progeny activity scores and add them to our Seurat object
  ## as a new assay called Progeny. 
  obj.list[[i]] <- progeny(obj.list[[i]], scale=FALSE, organism="Human", top=500, perm=1, 
                           return_assay = TRUE)
  
  ## We can now directly apply Seurat functions in our Progeny scores. 
  ## For instance, we scale the pathway activity scores. 
  obj.list[[i]] <- Seurat::ScaleData(obj.list[[i]], assay = "progeny") 
  
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(obj.list[[i]], slot = "scale.data", 
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
  summarized_progeny_scores_df.list[[i]] <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
  
}

for (i in 1:3) {
  row.names(summarized_progeny_scores_df.list[[i]]) <- 
    paste(row.names(summarized_progeny_scores_df.list[[i]]),
          c("Cakir", "Giandomenico", "Qian")[i], sep = "_")
}


# Plot Progeny scores ----------------------------------------------------------
source("https://raw.githubusercontent.com/vertesy/Rocinante/main/R/Rocinante.R")

# Plot Progeny Cakir
myColor = matlabColors.pheatmap(t(summarized_progeny_scores_df.list[[1]][,-1]), 100)
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[1]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[1]][,-1])))))
annot_df["Stressed",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Cakir_Progeny.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[1]][,-1]),
         fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "none", 
         #breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()

# Plot Progeny Giandomenico
myColor = matlabColors.pheatmap(t(summarized_progeny_scores_df.list[[2]][,-1]), 100)
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[2]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[2]][,-1])))))
annot_df["Stressed",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Giandomenico_Progeny.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[2]][,-1]),
         fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "none", 
         #breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()

# Plot Progeny Qian
myColor = matlabColors.pheatmap(t(summarized_progeny_scores_df.list[[3]][,-1]), 100)
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[3]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[3]][,-1])))))
annot_df["Stressed_Qian",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Qian_Progeny.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[3]][,-1]),
         fontsize=14, 
         fontsize_row = 10, 
         color=myColor,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "none", 
         #breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()

# Use black map and scale
colors2 <- colorRampPalette(c("white", "black"))(100)
# Plot Progeny Cakir
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[1]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[1]][,-1])))))
annot_df["Stressed",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Cakir_Progeny_scale.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[1]][,c(4,8,2:3,5:7,9:14)]),fontsize=14, 
         fontsize_row = 10, 
         color=colors2,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "row", cluster_rows = F,
         breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()

# Plot Progeny Giandomenico
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[2]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[2]][,-1])))))
annot_df["Stressed",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Giandomenico_Progeny_scale.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[2]][,c(4,8,2:3,5:7,9:14)]),fontsize=14, 
         fontsize_row = 10, 
         color=colors2,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "row", cluster_rows = F,
         breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()

# Plot Progeny Qian
annot_df <- data.frame(row.names = colnames(t(summarized_progeny_scores_df.list[[3]][,-1])), 
                       `Cell Type` = rep("non stressed", length(colnames(t(summarized_progeny_scores_df.list[[3]][,-1])))))
annot_df["Stressed_Qian",] = "stressed"

pdf('Gruffi_revision/Figure_OptimizationDatasets/Qian_Progeny_scale.pdf')
pheatmap(t(summarized_progeny_scores_df.list[[3]][,c(4,8,2:3,5:7,9:14)]),fontsize=14, 
         fontsize_row = 10, 
         color=colors2,  
         annotation_col = annot_df, annotation_colors = list(`Cell.Type` = c(`non stressed` = "grey", stressed = "dark red")),
         scale = "row", cluster_rows = F,
         breaks = seq(-2.5,2.5, 0.05),
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 10,  treeheight_row = 10,
         border_color = NA, cellwidth = 20, cellheight = 20)
dev.off()


summary_df <- data.frame(
  Pct.Stressed =
    c(obj.list[[1]]@meta.data %>%
        count(is.Stressed) %>%
        summarise(Pct.Stressed = round(n/sum(n)*100,2)) %>%
        magrittr::use_series(Pct.Stressed),
      obj.list[[2]]@meta.data %>%
        count(is.Stressed) %>%
        summarise(Pct.Stressed = round(n/sum(n)*100,2)) %>%
        magrittr::use_series(Pct.Stressed),
      obj.list[[3]]@meta.data %>%
        count(is.Stressed) %>%
        summarise(Pct.Stressed = round(n/sum(n)*100,2)) %>%
        magrittr::use_series(Pct.Stressed)),
  Dataset = c(rep(c("Cakir", "Giandomenico", "Qian"),each=2)),
  is.Stressed = rep(c(FALSE, TRUE), 3)
)

pdf(paste0(OutDir, "PctStressed_bar.pdf"), width = 5, height = 5)
ggbarplot(summary_df, x = "Dataset", y = "Pct.Stressed", fill = "is.Stressed", 
          color = NA, palette = c("grey", "dark red"), label = T, lab.pos = "in", lab.col = "white",
          ylab = "Pct of dataset", width = 0.75, ggtheme = theme_pubr(x.text.angle = 45)) + coord_fixed(0.05)
dev.off()

pdf(paste0(OutDir, "PctStressed_bar_all.pdf"), width = 20, height = 5)

sc.obj.pre@meta.data %>%
  group_by(ShortNames) %>%
  count(is.Stressed) %>%
  mutate(Pct.Stressed = round(n/sum(n)*100,2)) %>%
  ggbarplot(x = "ShortNames", y = "Pct.Stressed", fill = "is.Stressed", 
            color = NA, palette = c("grey", "dark red"), label = T, lab.pos = "in", lab.col = "white",
            ylab = "Pct of dataset", xlab = "",
            width = 0.75, ggtheme = theme_pubr(x.text.angle = 45, legend = "none")) + coord_fixed(0.05) +
  geom_hline(yintercept = c(6,15))
dev.off()


saveRDS(obj.list[[1]], paste0(OutDir, "Cakir_processed_plotted.Rds"))
saveRDS(obj.list[[2]], paste0(OutDir, "Giandomenico_processed_plotted.Rds"))
saveRDS(obj.list[[3]], paste0(OutDir, "Qian_processed_plotted.Rds"))




