# Integrate Vel5, Vel6, Kan25 and Kan26 to check were non-neuronal cells go

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr)

OutDir <- 'Gruffi_revision/Int_sc.obj.integration_Vel6/'
dir.create(OutDir)
# Load datasets ----------------------------------------------------------------
Vel5 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Vel.Org5.d90.ImpV.201029_non_neuronal.Rds')
Vel6 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Vel.Org6.d90.ImpV.201029_non_neuronal.Rds')
Kan25 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Kan.25.d120.ImpK.201020_non_neuronal.Rds')
Kan26 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Kan.26.d120.ImpK.201015_non_neuronal.Rds')


# preprocess for integration ---------------------------------------------------

object.list <- list(Vel5, Vel6, Kan25, Kan26)

for (i in 1: length(object.list)) {
  object.list[[i]][["percent.mito"]] <- PercentageFeatureSet(object.list[[i]], pattern = "^MT\\.")
  
}


for (i in 1: length(object.list)) {
  object.list[[i]] <- subset(object.list[[i]], subset = nFeature_RNA >500)
  object.list[[i]] <- subset(object.list[[i]], subset = percent.mito <15)
  
}

for (i in 1: length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], 
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = 2000)
}

# perform integration ----------------------------------------------------------
object.list_int_anchors <- FindIntegrationAnchors(object.list = object.list, 
                                                  dims = 1:50)


sc.obj.integration <- IntegrateData(anchorset = object.list_int_anchors, 
                                    dims = 1:50)

saveRDS(sc.obj.integration, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_integration.Rds"))

# scale data and calc. visualis. -----------------------------------------------


sc.obj.integration <- FindVariableFeatures(sc.obj.integration)
sc.obj.integration <- ScaleData(sc.obj.integration, verbose = FALSE)
sc.obj.integration <- RunPCA(sc.obj.integration, npcs = 50, features = VariableFeatures(sc.obj.integration))


sc.obj.integration <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = sc.obj.integration, nPCs = 50, dimensions=3:2, reduction="umap")

sc.obj.integration <- FindNeighbors(sc.obj.integration,dims = 1:50)
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = c(.1,.2,.3,.5))

DimPlot(sc.obj.integration, group.by = "predicted.clusters")
DimPlot(sc.obj.integration, group.by = "predicted.clusters", reduction = "pca")

# Run Gruffi on integration ----------------------------------------------------


ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

sc.obj.integration <- aut.res.clustering(obj = sc.obj.integration)

granule.res.4.gruffi <- sc.obj.integration@misc$gruffi$'optimal.granule.res'	

sc.obj.integration <- reassign.small.clusters(sc.obj.integration, ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

require(MarkdownHelpers)

# Glycolytic process	GO:0006096
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = go1, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = go2, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
sc.obj.integration <- GO_score_evaluation(obj = sc.obj.integration, GO_term = go3, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)


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


pl1 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'is.Stressed', label =F, save.plot = F, splitby = "orig.ident")
pl2 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'predicted.clusters', label =F, save.plot = F)
cowplot::plot_grid(pl1,pl2)
FeaturePlot(sc.obj.integration,features = "DDIT4")

cellIDs.keep <- which_names(!sc.obj.integration$'is.Stressed')
subset.obj <- subset(x = sc.obj.integration, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj, save.plot = F)
subset.obj <- FindVariableFeatures(subset.obj)
subset.obj <- ScaleData(subset.obj)
subset.obj <- RunPCA(subset.obj, features = VariableFeatures(subset.obj))
subset.obj <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = subset.obj, nPCs = 50, dimensions=3:2, reduction="umap")
Seurat.utils::clUMAP('predicted.clusters', label = F, obj = subset.obj, save.plot = F)
Seurat.utils::clUMAP('integrated_snn_res.0.1', label = F, obj = subset.obj, save.plot = F)
DimPlot(subset.obj, reduction = "pca", group.by = "integrated_snn_res.0.1")

saveRDS(sc.obj.integration, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_integration_proc.Rds"))
sc.obj.integration <- readRDS(file = paste0("Gruffi_revision/Int_Vel5_Vel6/Vel5_Vel6_Kan25_Kan26_integration_proc.Rds"))
saveRDS(subset.obj, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_integration.Rds"))

# Compare distance in PCA space ------------------------------------------------

labels = c("immature EN", "DL-EN", "UL-EN", "Glia", "Retina/Adrenocortical", "IN", "Immune Cells")
stdevs_int <- sc.obj.integration@reductions$pca@stdev
names(stdevs_int) <- colnames(sc.obj.integration@reductions$pca@cell.embeddings)

dist_pca_Vel5_int <- sc.obj.integration@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(sc.obj.integration@meta.data %>% mutate(bc = row.names(.))) %>%
  filter(orig.ident == "Org5_Velasco") %>%
  dplyr::group_by(RNA_snn_res.0.1) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_int[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$RNA_snn_res.0.1)) %>% dplyr::select(-1) %>%
  dist() %>% as.matrix()
colnames(dist_pca_Vel5_int) <- rownames(dist_pca_Vel5_int) <- labels
annot_df <- data.frame(row.names= labels, `Cell Type` = c(rep("Telencephalon",4), "Mis-differentiated", "Telencephalon", "Mis-differentiated"))




pl4 <- pheatmap(dist_pca_Vel5_int, cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         #show_rownames = F, show_colnames = F,
         border_color = NA,
         annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         #clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Weighted euclidean distance in top 50 PCs\nVel5 Post-Integration", cellwidth = 20, cellheight = 20)

dev.off()
pdf('Gruffi_revision/Int_Vel5_Vel6/Vel5_PCA_dist_postIntegration.pdf', width = 7.5, height = 5)
pl4
dev.off()

stdevs_Vel5 <- Vel5@reductions$pca@stdev
names(stdevs_Vel5) <- colnames(Vel5@reductions$pca@cell.embeddings)
dist_pca_Vel5 <- Vel5@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(Vel5@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(RNA_snn_res.0.1) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_Vel5[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$RNA_snn_res.0.1)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()

colnames(dist_pca_Vel5) <- rownames(dist_pca_Vel5) <- labels

pl2 <- pheatmap(dist_pca_Vel5, cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         #show_rownames = F, show_colnames = F,
         border_color = NA,
         annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Weighted euclidean distance in top 50 PCs\nVel5 Single", cellwidth = 20, cellheight = 20)

pdf('Gruffi_revision/Int_Vel5_Vel6/Vel5_PCA_dist_preIntegration.pdf', width = 7.5, height = 5)
pl2
dev.off()

Vel5$RNA_snn_res.0.1_named <- factor(Vel5$RNA_snn_res.0.1, labels = labels)

pl1 <-clUMAP(ident = "RNA_snn_res.0.1_named", label = T, cols = pals::tol(7),label.cex = 5,
       obj = Vel5, save.plot = F, title = "res 0.1 clustering Vel5 pre-integration")

png('Gruffi_revision/Int_Vel5_Vel6/Vel5_UMAP_preInt.png', width = 500, height = 500)
pl1
dev.off()


Vel5_fromInt <- sc.obj.integration[,names(sc.obj.integration$orig.ident[sc.obj.integration$orig.ident=="Org5_Velasco"])]
Vel6_fromInt <- sc.obj.integration[,names(sc.obj.integration$orig.ident[sc.obj.integration$orig.ident=="Org6_Velasco"])]

Vel5_fromInt$RNA_snn_res.0.1_named <- factor(Vel5_fromInt$RNA_snn_res.0.1, labels = labels)

pl3 <- clUMAP(ident = "RNA_snn_res.0.1_named", label = F, cols = pals::tol(7), label.cex = 5,
       obj = Vel5_fromInt, save.plot = F, title = "res 0.1 clustering Vel5 post-integration")
png('Gruffi_revision/Int_Vel5_Vel6/Vel5_UMAP_postInt_nolab.png', width = 500, height = 500)
pl3
dev.off()


png('Gruffi_revision/Int_Vel5_Vel6/Int_origident_nolab.png', width = 500, height = 500)
clUMAP(ident = "orig.ident", label = F, cols = pals::tol(7), label.cex = 5,
       obj = sc.obj.integration, save.plot = F)
dev.off()

png('Gruffi_revision/Int_Vel5_Vel6/Int_Kan25_Kan26_Vel5_Vel6_nonNDistortion.png', width = 1000, height = 1000)
cowplot::plot_grid(pl1,pl2[[4]],pl3,pl4[[4]],nrow = 2,ncol=2)
dev.off()
pdf('Gruffi_revision/Int_Vel5_Vel6/Int_Kan25_Kan26_Vel5_Vel6_nonNDistortion.pdf', width = 15, height = 15)
cowplot::plot_grid(pl1,pl2[[4]],pl3,pl4[[4]],nrow = 2,ncol=2)
dev.off()


# Plot Stress pathway ----------------------------------------------------------
pl1 <- qUMAP(feature = c("Score.GO.0006096"),obj = Vel5_fromInt, save.plot = F,
      nr.cols = 1, cols = c(alpha("grey", 0.5), "red"))+NoAxes() +NoLegend()
pl2 <- qUMAP(feature = c("Score.GO.0034976"),obj = Vel5_fromInt, save.plot = F,
      nr.cols = 1, cols = c(alpha("grey", 0.5), "red"))+NoAxes() +NoLegend()
pl3 <- qUMAP(feature = c("Score.GO.0042063"),obj = Vel5_fromInt, save.plot = F,
      nr.cols = 1, cols = c(alpha("grey", 0.5), "red")) +NoAxes() +NoLegend()

pl4 <- clUMAP(ident = "is.Stressed", obj = Vel5_fromInt, save.plot = F,label=F, cols = c("dark grey","dark red"))+NoAxes() +NoLegend()
png(filename = "Gruffi_revision/Int_Vel5_Vel6/Vel5_nonN_Score_onInt.png", width = 3000,height = 500, res = 200)
cowplot::plot_grid(pl1,pl2,pl3,pl4,ncol=4)
dev.off()

pl5 <- Vel5_fromInt@meta.data %>%
  group_by(RNA_snn_res.0.1_named) %>%
  count(is.Stressed) %>%
  ggbarplot(x = "RNA_snn_res.0.1_named", y = "n", fill = "is.Stressed",label = T, position = position_fill(),color = NA,
            lab.pos = "in",
            ylab = "Fraction of cells per cluster", xlab = "", ggtheme = theme_pubr(x.text.angle = 30))+scale_fill_manual(values = c("dark grey","dark red"))
  
pdf("Gruffi_revision/Int_Vel5_Vel6/Vel5_nonN_Score_onInt_stat.pdf", width = 5,height = 5)
pl5
dev.off()

