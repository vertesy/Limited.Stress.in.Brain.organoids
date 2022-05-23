# Distance between different cell types
# Oliver Eichmueller
# Sun May  1 07:22:29 2022 ------------------------------

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer);
require(qgraph); require(enrichplot)

OutDir <- 'Gruffi_revision/Dist_PostGruffi/'
dir.create(OutDir)

# load datasets ----------------------------------------------------------------
sc.obj.pre  <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration/Before.filtering/combined.obj_w.Gruffi_2022.02.11_11.25.Rds.gz')
sc.obj.post <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration/After.filtering/combined.obj_After.Gruffi.DGEA_2022.02.12_18.42.Rds.gz')

# Plot UMAPs -------------------------------------------------------------------

colors <- getDiscretePalette(
  ident.used = "integrated_snn_res.0.3.Manual.short", 
  obj = sc.obj.pre)
names(colors) <- names(table(sc.obj.pre$integrated_snn_res.0.3.Manual.short))[
  rev(order(table(sc.obj.pre$integrated_snn_res.0.3.Manual.short)))]
pdf(paste0(OutDir, "PrevsPostGruffi_LabelTransfer_UMAP.pdf"), width = 15, height = 5)
cowplot::plot_grid(
  clUMAP(ident = "integrated_snn_res.0.3.Manual.short", 
       obj = sc.obj.pre, save.plot = F, axes = F, cols = colors,
       title = "Pre-Gruffi"),
  clUMAP(ident = "integrated_snn_res.0.3.Manual.short.Before.Gruffi", 
       obj = sc.obj.post, save.plot = F, axes = F, cols = colors,
       title = "Pre-Gruffi"), 
  ncol = 2)
dev.off()


# Calculate Variable features of the RNA assay ---------------------------------
sc.obj.pre@active.assay <- "RNA"
sc.obj.post@active.assay <- "RNA"

sc.obj.pre <- FindVariableFeatures(sc.obj.pre)
sc.obj.post <- FindVariableFeatures(sc.obj.post)

# Average across After.Gruffi annotations --------------------------------------

avg.pre <- as.data.frame(
  AverageExpression(sc.obj.pre, assays = "RNA", 
                    features = VariableFeatures(sc.obj.pre),
                    group.by = "integrated_snn_res.0.3.Manual.short")$RNA)

avg.post <- as.data.frame(
  AverageExpression(sc.obj.post, assays = "RNA", 
                    features = VariableFeatures(sc.obj.post),
                    group.by = "integrated_snn_res.0.3.Manual.short.Before.Gruffi")$RNA)

# calculate euclidean distances ------------------------------------------------

dist.pre  <- as.matrix(dist(t(avg.pre)))
dist.post <- as.matrix(dist(t(avg.post)))

coi <- intersect(colnames(dist.pre), colnames(dist.post))

dist.rel <-
  round((dist.post[coi,coi] - dist.pre[coi,coi])/
          dist.pre[coi,coi]*100,2)

# plot in heatmap --------------------------------------------------------------
# Color palette
coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(100)
colnames(sc.obj.d125_fetal)[[1]]
# pct difference all clusters (adj. color ramp)
pdf(paste0(OutDir, "pct_dist_PrevsPost_all_clust.pdf"), width = 10, height = 10)
pheatmap(dist.rel, 
         cluster_rows = T, cluster_cols = T, 
         display_numbers = T, number_format = "%.1f",
         color = coul[c(1:33,seq(34,100,by=2))], 
         #gaps_row = 5, gaps_col = 5,
         cutree_rows = 3, cutree_cols = 3, 
         width = 2,height = 2, 
         border_color = NA,
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "% Difference of euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# pct difference all clusters (adj. color ramp)
pdf(paste0(OutDir, "pct_dist_PrevsPost_NoUncl_clust.pdf"), width = 10, height = 10)
pheatmap(dist.rel[-c(11, 15, 16),-c(11,15, 16)], 
         cluster_rows = T, cluster_cols = T, 
         display_numbers = T, number_format = "%.1f",
         color = coul[c(1:33,seq(34,100,by=2))], 
         #gaps_row = 5, gaps_col = 5,
         cutree_rows = 3, cutree_cols = 3, 
         width = 2,height = 2, 
         border_color = NA,
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "% Difference of euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()


# pct difference first selection (adj. color ramp)
coi.plot <- colnames(dist.rel)[c(2:10,12:14)][c(3,1,6,11,7,
                                                4,5,8,9,12,2,10)]
pdf(paste0(OutDir, "pct_dist_PrevsPost_select1.pdf"), width = 10, height = 10)
pheatmap(dist.rel[coi.plot,coi.plot], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul[c(1:49,seq(50,100,by=2))], 
         gaps_row = 5, gaps_col = 5,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA,
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "% Difference of euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# pct difference second selection (adj. color ramp)
coi.plot.subset <- colnames(dist.rel)[c(4,7,13,8,6,9,3,14,12)]
pdf(paste0(OutDir, "pct_dist_PrevsPost_select2.pdf"), width = 10, height = 10)
pheatmap(dist.rel[coi.plot.subset,coi.plot.subset], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul[c(1:49,seq(50,100,by=2))], 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "% Difference of euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all pre.Gruffi
pdf(paste0(OutDir, "pct_dist_Pre_all.pdf"), width = 10, height = 10)
pheatmap(dist.pre, 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Pre-Gruffi euclidean distance\n2000 top variable feature", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all pre.Gruffi first selection
pdf(paste0(OutDir, "pct_dist_Pre_select.pdf"), width = 10, height = 10)
pheatmap(dist.pre[coi.plot,coi.plot], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Pre-Gruffi euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all pre.Gruffi second selection
pdf(paste0(OutDir, "pct_dist_Pre_select2.pdf"), width = 10, height = 10)
pheatmap(dist.pre[coi.plot.subset,coi.plot.subset], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Pre-Gruffi euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all post.Gruffi
pdf(paste0(OutDir, "pct_dist_Post_all.pdf"), width = 10, height = 10)
pheatmap(dist.post, 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Post-Gruffi euclidean distance\n2000 top variable feature", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all post.Gruffi first selection
pdf(paste0(OutDir, "pct_dist_Post_select.pdf"), width = 10, height = 10)
pheatmap(dist.post[coi.plot,coi.plot], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Post-Gruffi euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# dist all post.Gruffi second selection
pdf(paste0(OutDir, "pct_dist_Post_select2.pdf"), width = 10, height = 10)
pheatmap(dist.post[coi.plot.subset,coi.plot.subset], 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         #breaks= c(-50:50),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "Post-Gruffi euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()

# Selection only terminal cell types
coi.plot.final <- colnames(dist.rel)[c(3,14,6,9,10)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(30)
plot_matrix <- dist.rel[coi.plot.final,coi.plot.final]
plot_matrix[upper.tri(plot_matrix)] <- NA

pdf(paste0("Gruffi_revision/Dist_PostGruffi/rel_dist_Post_selected_triangle.pdf"), width = 10, height = 10)
pheatmap(plot_matrix, 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = T, number_format = "%.1f",
         color = coul, 
         #gaps_row = 4, gaps_col = 4,
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA, 
         breaks= c(-15:15),
         # annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         main = "% Difference of euclidean distance\n2000 top variable features", 
         cellwidth = 20, cellheight = 20)
dev.off()


# Save Tables as xlsx ----------------------------------------------------------
openxlsx::write.xlsx(list(pre.Gruffi = dist.pre,
                          post.Gruffi = dist.post,
                          pct.diff = dist.rel,
                          genes.used = data.frame(pre = VariableFeatures(sc.obj.pre),
                                                  post = VariableFeatures(sc.obj.post))),
                     file = paste0(OutDir, "Dist_Tables_PrePost.xlsx"), rowNames = T)

