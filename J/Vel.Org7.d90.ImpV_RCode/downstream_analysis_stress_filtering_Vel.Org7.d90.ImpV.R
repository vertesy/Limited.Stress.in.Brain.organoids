rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)
devtools::load_all("./gruffi_dev/")

library(Seurat)
library(ggplot2)
library(Seurat.utils)

OutDir <- OutDirOrig <- "./Vel.Org7.d90.ImpV_RDS"
setwd(OutDirOrig)

combined.obj <- readRDS("./Vel.Org7.d90.ImpV.201029.Rds")

#LQ cell is filtered, if less then 1000 genes are detected within a cell
combined.obj <- subset(x = combined.obj, subset = `nFeature_RNA` > 1000)

combined.obj <- NormalizeData(object = combined.obj)
combined.obj <- FindVariableFeatures(object = combined.obj, mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR')
combined.obj <- ScaleData(combined.obj)
combined.obj <- RunPCA(combined.obj, npcs = 50, features = VariableFeatures(object = combined.obj))
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:50, n.components = 2L)

combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:50)
combined.obj <- aut.res.clustering(combined.obj, assay = "RNA")
combined.obj <- reassign.small.clusters(combined.obj, ident = "seurat_clusters")

# number of granules
length(unique(combined.obj$seurat_clusters.reassigned))
# median number of cells per granule
median(table(combined.obj$seurat_clusters.reassigned))

# Stress Assignment for Different Cluster Resolutions
for(res in c(1,6,10,20)) {
  combined.obj <- FindClusters(combined.obj, res = res)
  combined.obj <- reassign.small.clusters(combined.obj, ident = paste0("RNA_snn_res.",res))
  
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0006096", stat.av = "normalized.mean", #plot.each.gene = T,
                                      new_GO_term_computation = T, clustering = paste0("RNA_snn_res.",res,".reassigned"))  
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0034976", stat.av = "normalized.mean", #plot.each.gene = T,
                                      new_GO_term_computation = T, clustering = paste0("RNA_snn_res.",res,".reassigned")) 
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = "GO:0042063", stat.av = "normalized.mean", #plot.each.gene = T,
                                      new_GO_term_computation = T, clustering = paste0("RNA_snn_res.",res,".reassigned")) 
  
  Idents(combined.obj) <- combined.obj$orig.ident
  
  i1 <- paste0("RNA_snn_res.",res,".reassigned_cl.av_GO:0006096")
  i2 <- paste0("RNA_snn_res.",res,".reassigned_cl.av_GO:0034976")
  i3 <- paste0("RNA_snn_res.",res,".reassigned_cl.av_GO:0042063")
  
  # Interactive Stress Assignment Based on Granule GO Scores
  combined.obj <- Shiny.GO.thresh(stress.ident1 = i1, 
                                  stress.ident2 = i2, 
                                  notstress.ident3 = i3)
  
  combined.obj <- AddMetaData(combined.obj, combined.obj$is.Stressed, paste0("is.Stressed.res.",res,".reassigned")) 
}
saveRDS(combined.obj, "Vel.Org7.d90.ImpV.processed.RDS")

Idents(combined.obj) <- combined.obj$is.Stressed.res.6.reassigned
combined.obj.filtered <- subset(combined.obj, ident = FALSE)
rm(combined.obj)

combined.obj.filtered <- FindVariableFeatures(object = combined.obj.filtered, mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR')
combined.obj.filtered <- ScaleData(combined.obj.filtered)
combined.obj.filtered <- RunPCA(combined.obj.filtered, npcs = 50, features = VariableFeatures(object = combined.obj.filtered))
combined.obj.filtered <- RunUMAP(combined.obj.filtered, reduction = "pca", dims = 1:50, n.components = 2L)

saveRDS(combined.obj.filtered, "Vel.Org7.d90.ImpV.stressfiltered.RDS")
