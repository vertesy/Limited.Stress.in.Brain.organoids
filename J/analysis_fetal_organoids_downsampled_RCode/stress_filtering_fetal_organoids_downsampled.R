rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
devtools::load_all("./gruffi_dev/")

library(Seurat)
library(ggplot2)
library(Seurat.utils)

OutDir <- OutDirOrig <- "./analysis_fetal_organoids_downsampled_RDS/" 
setwd(OutDir)

combined.obj <- readRDS("referencemapped.invivo.SEO.downsampled.24211.RDS")

combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:50)
combined.obj <- aut.res.clustering(combined.obj)
combined.obj <- reassign.small.clusters(combined.obj, ident = "seurat_clusters")

# number of granules
length(unique(combined.obj$seurat_clusters.reassigned))
# median number of cells per granule
median(table(combined.obj$seurat_clusters.reassigned))

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

# Total Number & Ratio of Stressed Cells in Integrated Object 
table(combined.obj$is.Stressed)
length(which(combined.obj$is.Stressed))/dim(combined.obj)[2]

combined.obj$orig.ident <- "in vivo"
combined.obj$orig.ident[which(is.na(combined.obj$predicted.cell.type))] <- "SEO"

# Total Number & Ratio of Stressed Cells in Geschwind Fetal Data 
length(which(combined.obj$is.Stressed & combined.obj$orig.ident == "in vivo"))
length(which(combined.obj$is.Stressed & combined.obj$orig.ident == "in vivo"))/length(which(combined.obj$orig.ident == "in vivo"))

# Total Number & Ratio of Stressed Cells in Geschwind Fetal Data 
length(which(combined.obj$is.Stressed & combined.obj$orig.ident == "SEO"))
length(which(combined.obj$is.Stressed & combined.obj$orig.ident == "SEO"))/length(which(combined.obj$orig.ident == "SEO"))

saveRDS(combined.obj, "combined.obj.processed.invivo.SEO.downsampled.24211.RDS")
