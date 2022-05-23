rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
devtools::load_all("./gruffi_dev/")

library(Seurat)

# Setup ------------------------
OutDir <- OutDirOrig <- paste0("./downsampled_organoid_integration_RDS/")
setwd(OutDirOrig)

combined.obj <- readRDS("./b4.gruffi.ribolow.30K__2021.12.10_19.53.Rds")
sub.combined.obj <- subset(x = combined.obj, cells = sample(x = colnames(combined.obj), size = 24211))
combined.obj <- sub.combined.obj
DefaultAssay(combined.obj) <- "integrated"

combined.obj <- FindVariableFeatures(combined.obj)
combined.obj <- ScaleData(combined.obj)
combined.obj <- RunPCA(combined.obj, npcs = 50, features = VariableFeatures(object = combined.obj))
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:50, n.components = 2L)

combined.obj <- FindNeighbors(combined.obj)
combined.obj <- aut.res.clustering(combined.obj)
combined.obj <- reassign.small.clusters(combined.obj, ident = "seurat_clusters") 

# Computation of (Cellwise and) Granulewise GO Scores
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0006096", stat.av = "normalized.mean", 
                                    new_GO_term_computation = T, clustering = "seurat_clusters.reassigned")#, plot.each.gene = T)
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0034976", stat.av = "normalized.mean", 
                                    new_GO_term_computation = T, clustering = "seurat_clusters.reassigned")#, plot.each.gene = T)
combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0042063", stat.av = "normalized.mean", 
                                    new_GO_term_computation = T, clustering = "seurat_clusters.reassigned")#, plot.each.gene = T)

# Interactive Stress Assignment Based on Granule GO Scores
i1 <- "seurat_clusters.reassigned_cl.av_GO:0006096"
i2 <- "seurat_clusters.reassigned_cl.av_GO:0034976"
i3 <- "seurat_clusters.reassigned_cl.av_GO:0042063"

combined.obj <- Shiny.GO.thresh(stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3, plot.cluster.shiny = "orig.ident")

# Stress Assignment Based on Cellwise GO Scores
i1 <- paste0("Score.GO.0006096")
i2 <- paste0("Score.GO.0034976")
i3 <- paste0("Score.GO.0042063")

thresh.stress.ident1 <- plot_norm_and_skew(combined.obj$Score.GO.0006096, q = .99, tresholding = "fitted", plot.hist = F)
thresh.stress.ident2 <- plot_norm_and_skew(combined.obj$Score.GO.0034976, q = .99, tresholding = "fitted", plot.hist = F)
thresh.notstress.ident3 <- plot_norm_and_skew(combined.obj$Score.GO.0042063, q = .99, tresholding = "fitted", plot.hist = F)

combined.obj$is.Stressed.cellwise <- ((combined.obj$Score.GO.0006096 > thresh.stress.ident1)|(combined.obj$Score.GO.0034976 > thresh.stress.ident2)) & 
  !(combined.obj$Score.GO.0042063 > thresh.notstress.ident3)

# Total Number of Cells Classified as Stressed Cellwise vs Granulewise
combined.obj$is.Stressed.only.cellwise <- combined.obj$is.Stressed.cellwise & !combined.obj$is.Stressed
combined.obj$is.Stressed.only.granule <- combined.obj$is.Stressed & !combined.obj$is.Stressed.cellwise

combined.obj$stress.ident <- FALSE
combined.obj$stress.ident[combined.obj$is.Stressed.only.granule == TRUE] <- "only.granule"
combined.obj$stress.ident[combined.obj$is.Stressed.only.cellwise == TRUE] <- "only.cellwise"
combined.obj$stress.ident[(combined.obj$is.Stressed == TRUE & !combined.obj$is.Stressed.only.granule == TRUE) |
                            (combined.obj$is.Stressed.cellwise == TRUE & !combined.obj$is.Stressed.only.cellwise == TRUE)] <- "both"
table(combined.obj$stress.ident)

# Differential Gene Expression Only Cellwise vs Only Granulewise Stress Annootation
DefaultAssay(combined.obj) <- "RNA"
Idents(combined.obj) <- combined.obj$stress.ident
marker.stress <- FindMarkers(combined.obj, ident.1 = "only.granule", ident.2 = "only.cellwise")
marker.stress$gene <- rownames(marker.stress)
write.table(marker.stress, "marker.granule.cellwise.stress.csv", sep = ";", row.names = F, dec = ",")

saveRDS(combined.obj, "combined.obj.processed.SEO.downsampled.24211.RDS")
