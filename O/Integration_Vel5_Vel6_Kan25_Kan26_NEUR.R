# Integrate only neuronal of Vel5, Vel6, Kan25 and Kan26 

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr)

OutDir <- 'Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/'
dir.create(OutDir)
# Load datasets ----------------------------------------------------------------
InputDir <- '/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/NeuronalFiltered/'
files <- list.files(InputDir)
object.list <- list()
for (i in 1:length(files)) {
  object.list[[i]] <- readRDS(paste0(InputDir, files[[i]]))
}

# preprocess for integration ---------------------------------------------------

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

saveRDS(sc.obj.integration, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_NEUR_integration.Rds"))

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


pl1 <- Seurat.utils::clUMAP(obj = sc.obj.integration, 'is.Stressed', label =F, save.plot = F, cols = c("dark grey","dark red"))+NoAxes() +NoLegend()
pl2 <- Seurat.utils::clUMAP(ident = "orig.ident", label = F, cols = pals::tol(7), label.cex = 5,
              obj = sc.obj.integration, save.plot = F) +NoAxes() +NoLegend()
cowplot::plot_grid(pl1,pl2)

png(filename = "Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/Int_neur_origident.png", width = 1000,height = 500, res = 200)
pl2
dev.off()
png(filename = "Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/Int_neur_Stress.png", width = 1000,height = 500, res = 200)
pl1
dev.off()

cellIDs.keep <- which_names(!sc.obj.integration$'is.Stressed')
subset.obj <- subset(x = sc.obj.integration, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj, save.plot = F)

subset.obj <- FindVariableFeatures(subset.obj)
subset.obj <- ScaleData(subset.obj)
subset.obj <- RunPCA(subset.obj, features = VariableFeatures(subset.obj))
subset.obj <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = subset.obj, nPCs = 50, dimensions=3:2, reduction="umap")

subset.obj <- FindNeighbors(subset.obj,dims = 1:50)
subset.obj <- FindClusters(subset.obj, resolution = c(.1,.2,.3,.5))

Seurat.utils::clUMAP('predicted.clusters.cleaned', label = F, obj = subset.obj, save.plot = F)
Seurat.utils::clUMAP('integrated_snn_res.0.3', label = T, obj = subset.obj, save.plot = F)

cols1 <- c(colors[["immature EN"]], colors[["UL-EN"]], colors[["DL-EN"]],  "grey", "grey",colors[["IN"]], "grey", "grey", "grey", "grey", "grey", "grey")

pl1 <- Seurat.utils::clUMAP(ident = "integrated_snn_res.0.3", label = F, cols = cols1, label.cex = 5,
                            obj = subset.obj, save.plot = F, title = "Post-Gruffi old integration") +NoAxes() +NoLegend()
pl0 <- Seurat.utils::clUMAP(ident = "integrated_snn_res.0.3_POSTGRUFFI", label = F, cols = c(cols1,"dark red"), label.cex = 5,
                            obj = sc.obj.integration.PRE, save.plot = F, title = "Pre-Gruffi old integration") +NoAxes() +NoLegend()

DimPlot(subset.obj, reduction = "pca", group.by = "predicted.cell.type.neuronal")

saveRDS(sc.obj.integration, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_NEUR_integration_proc.Rds"))
sc.obj.integration.PRE <- readRDS(file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_NEUR_integration_proc.Rds"))
saveRDS(subset.obj, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_NEUR_integration_subset.Rds"))


sc.obj.integration.PRE$integrated_snn_res.0.3_POSTGRUFFI <- sc.obj.integration.PRE$is.Stressed
sc.obj.integration.PRE$integrated_snn_res.0.3_POSTGRUFFI[names(subset.obj$integrated_snn_res.0.3)] <-
  as.vector(subset.obj$integrated_snn_res.0.3)
sc.obj.integration.PRE$integrated_snn_res.0.3_POSTGRUFFI <-
  factor(sc.obj.integration.PRE$integrated_snn_res.0.3_POSTGRUFFI, levels = c(0:11,"TRUE"))


# Re integrate after Stress annotation =========================================
subset.obj@active.assay <- "RNA"
subset.obj.split <- SplitObject(subset.obj, split.by = "orig.ident")

object.list_int_anchors <- FindIntegrationAnchors(object.list = subset.obj.split, 
                                                  dims = 1:50)


sc.obj.integration <- IntegrateData(anchorset = object.list_int_anchors, 
                                    dims = 1:50)

sc.obj.integration <- FindVariableFeatures(sc.obj.integration)
sc.obj.integration <- ScaleData(sc.obj.integration, verbose = FALSE)
sc.obj.integration <- RunPCA(sc.obj.integration, npcs = 50, features = VariableFeatures(sc.obj.integration))


sc.obj.integration <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = sc.obj.integration, nPCs = 50, dimensions=3:2, reduction="umap")

sc.obj.integration <- FindNeighbors(sc.obj.integration,dims = 1:50)
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = c(.1,.2,.3,.5))
sc.obj.integration <- FindClusters(sc.obj.integration, resolution = c(.25))

DimPlot(sc.obj.integration, group.by = "orig.ident")
DimPlot(sc.obj.integration, group.by = "integrated_snn_res.0.25", label = T)
DimPlot(sc.obj.integration, group.by = "predicted.cell.type.neuronal", reduction = "umap")

sc.obj.integration$integrated_snn_res.0.3_POSTGRUFFI <-
  subset.obj$integrated_snn_res.0.3[names(sc.obj.integration$integrated_snn_res.0.3)]


cols2 <- c(colors[["UL-EN"]], colors[["immature EN"]], colors[["DL-EN"]],  "grey",colors[["IN"]], "grey", "grey", "grey", "grey", "grey")

pl2 <- Seurat.utils::clUMAP(ident = "integrated_snn_res.0.3_POSTGRUFFI", label = F, cols = cols1, label.cex = 5,
                            obj = sc.obj.integration, save.plot = F, title = "Post-Gruffi re-integrated") +NoAxes() +NoLegend()
pl1+pl2
subset.obj@active.assay <- "integrated"
setdiff(VariableFeatures(subset.obj),VariableFeatures(sc.obj.integration))

setdiff(row.names(sc.obj.integration@assays$integrated), row.names(subset.obj@assays$integrated))

FeaturePlot()

saveRDS(sc.obj.integration, file = paste0(OutDir, "Vel5_Vel6_Kan25_Kan26_NEUR_postGRUFFIintegration.Rds"))




stdevs_int_RE <- sc.obj.integration@reductions$pca@stdev
names(stdevs_int_RE) <- colnames(sc.obj.integration@reductions$pca@cell.embeddings)
dist_pca_int_RE <- sc.obj.integration@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(sc.obj.integration@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(integrated_snn_res.0.3_POSTGRUFFI) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_int_RE[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$integrated_snn_res.0.3_POSTGRUFFI)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()
coi <- c("Cl_0","Cl_1", "Cl_2", "Cl_5")
dist_pca_int_RE_select <- dist_pca_int_RE[coi,coi]
colnames(dist_pca_int_RE_select) <- rownames(dist_pca_int_RE_select) <- c("immature EN","UL-EN", "DL-EN", "IN")

stdevs_int_post <- subset.obj@reductions$pca@stdev
names(stdevs_int_post) <- colnames(subset.obj@reductions$pca@cell.embeddings)
dist_pca_int_post <- subset.obj@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(subset.obj@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(integrated_snn_res.0.3) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_int_post[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$integrated_snn_res.0.3)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()

dist_pca_int_post_select <- dist_pca_int_post[coi,coi]
colnames(dist_pca_int_post_select) <- rownames(dist_pca_int_post_select) <- c("immature EN","UL-EN", "DL-EN", "IN")

stdevs_int_pre <- sc.obj.integration.PRE@reductions$pca@stdev
names(stdevs_int_pre) <- colnames(sc.obj.integration.PRE@reductions$pca@cell.embeddings)
dist_pca_int_pre <- sc.obj.integration.PRE@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(sc.obj.integration.PRE@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(integrated_snn_res.0.3_POSTGRUFFI) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_int_pre[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$integrated_snn_res.0.3_POSTGRUFFI)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()

dist_pca_int_pre_select <- dist_pca_int_pre[coi,coi]
colnames(dist_pca_int_pre_select) <- rownames(dist_pca_int_pre_select) <- c("immature EN","UL-EN", "DL-EN", "IN")



diff_dist_post_REintegration <- round((dist_pca_int_RE_select[colnames(dist_pca_int_RE_select),colnames(dist_pca_int_RE_select)]-
  dist_pca_int_post_select[colnames(dist_pca_int_RE_select),colnames(dist_pca_int_RE_select)])/
  dist_pca_int_post_select[colnames(dist_pca_int_RE_select),colnames(dist_pca_int_RE_select)]*100,2)
coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(25)
pl3 <- pheatmap(diff_dist_post_REintegration, cluster_rows = T, cluster_cols = T, display_numbers = T, color = coul, 
                cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, number_format = "%.1f",number_color = "black",
                #show_rownames = F, show_colnames = F,
                border_color = NA, breaks = (-12:12),treeheight_row = 20, treeheight_col = 20,
                #annotation_row = annot_df, annotation_col = annot_df,
                #annotation_colors = list(predicted.clusters = colors),
                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                main = "% difference re-integration\nweighted euclidean distance\ntop 50 PCs", cellwidth = 20, cellheight = 20)

diff_dist_pre_post <- round((dist_pca_int_post_select-
                                        dist_pca_int_pre_select)/
                                       dist_pca_int_pre_select*100,2)
diff_dist_pre_REintegration <- round((dist_pca_int_RE_select-
                                        dist_pca_int_pre_select)/
                                       dist_pca_int_pre_select*100,2)
coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(25)
pl3 <- pheatmap(diff_dist_post_REintegration, 
                cluster_rows = F, cluster_cols = F, 
                display_numbers = T, color = coul, 
                cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, number_format = "%.1f",number_color = "black",
                #show_rownames = F, show_colnames = F,
                border_color = NA, breaks = (-16:16),treeheight_row = 20, treeheight_col = 20,
                #annotation_row = annot_df, annotation_col = annot_df,
                #annotation_colors = list(predicted.clusters = colors),
                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                main = "% difference\nre-integration vs post-Gruffi", cellwidth = 20, cellheight = 20)

coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(33)
pl4 <- pheatmap(diff_dist_pre_REintegration, 
                cluster_rows = F, cluster_cols = F, 
                display_numbers = T, color = coul, 
                cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, number_format = "%.1f",number_color = "black",
                #show_rownames = F, show_colnames = F,
                border_color = NA, breaks = (-16:16),treeheight_row = 20, treeheight_col = 20,
                #annotation_row = annot_df, annotation_col = annot_df,
                #annotation_colors = list(predicted.clusters = colors),
                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                main = "% difference\nre-integration vs pre-Gruffi", cellwidth = 20, cellheight = 20)

pl5 <- pheatmap(diff_dist_pre_post, 
                cluster_rows = F, cluster_cols = F, 
                display_numbers = T, color = coul, 
                cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, number_format = "%.1f",number_color = "black",
                #show_rownames = F, show_colnames = F,
                border_color = NA, breaks = (-16:16),treeheight_row = 20, treeheight_col = 20,
                #annotation_row = annot_df, annotation_col = annot_df,
                #annotation_colors = list(predicted.clusters = colors),
                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
                main = "% difference\npost-Gruffi vs pre-Gruffi", cellwidth = 20, cellheight = 20)

pdf("Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/Diff_Distance_REvsPOST.pdf", width = 5,height = 5)
pl3
dev.off()
pdf("Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/Diff_Distance_REvsPRE.pdf", width = 5,height = 5)
pl4
dev.off()
pdf("Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/Diff_Distance_POSTvsPRE.pdf", width = 5,height = 5)
pl5
dev.off()

png(filename = "Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PREGruffi_Res03_PostDefined.png", width = 1000,height = 500, res = 200)
pl0
dev.off()
png(filename = "Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/POSTGruffi_Res03_PostDefined.png", width = 1000,height = 500, res = 200)
pl1
dev.off()
png(filename = "Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/REintegrated_Res03_PostDefined.png", width = 1000,height = 500, res = 200)
pl2
dev.off()



# Compare integration features
library(clusterProfiler)
library(enrichplot)
install.packages("GOCompare")
require(gprofiler2)

PRE_int <- enrichGO(row.names(sc.obj.integration.PRE@assays$integrated),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

POST_int <- enrichGO(row.names(sc.obj.integration@assays$integrated),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

POST_int_new <- enrichGO(setdiff(row.names(sc.obj.integration@assays$integrated),
                                 row.names(sc.obj.integration.PRE@assays$integrated)),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

PRE_int_new <- enrichGO(setdiff(row.names(sc.obj.integration.PRE@assays$integrated),
                                row.names(sc.obj.integration@assays$integrated)),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")


openxlsx::write.xlsx(list(Pre = PRE_int@result %>% filter(p.adjust <0.05),
                      Post = POST_int@result %>% filter(p.adjust <0.05)),
                file = "/Users/Oliver.Eichmueller/Desktop/GO_export.xlsx")

goi_post <- HVFInfo(sc.obj.integration) %>% arrange(desc(mean)) %>% mutate(gene=row.names(.)) %>%
  select(gene,mean, variance.standardized)

xlsx::write.xlsx(goi_post,
                file = "/Users/Oliver.Eichmueller/Desktop/GOIpost_export.xlsx")

goi_pre <- HVFInfo(sc.obj.integration.PRE) %>% arrange(desc(mean)) %>% mutate(gene=row.names(.)) %>%
  select(gene,mean, variance.standardized)

xlsx::write.xlsx(goi_pre,
                file = "/Users/Oliver.Eichmueller/Desktop/GOIpre_export.xlsx")

cluster <- compareCluster(geneClusters = list(pre=PRE_int@gene, post = POST_int@gene), fun = "enrichGO",
                          OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP")
enrichplot::dotplot(cluster, showCategory=40)
library(DOSE)
cluster2 <- pairwise_termsim(cluster, showCategory = 2000)
emapplot(cluster2, vertex.label.cex=1.2, showCategory = 30)
emapplot(pairwise_termsim(PRE_int_new), vertex.label.cex=1.2)
emapplot(pairwise_termsim(PRE_int), vertex.label.cex=1.2)+emapplot(pairwise_termsim(POST_int), vertex.label.cex=1.2)

ID_unique <- cluster2@compareClusterResult %>% dplyr::select(Cluster, ID, GeneRatio, p.adjust) %>% group_by(ID) %>%
  count(ID) %>% ungroup() %>% filter(n==1) %>% magrittr::use_series(ID)



Desc_found <- cluster2@compareClusterResult %>% dplyr::select(Cluster, ID, GeneRatio, p.adjust,Description)%>%
  filter(ID %in% ID_unique) %>% magrittr::use_series(Description)
labs <- cluster2@compareClusterResult %>% dplyr::select(Cluster, ID, GeneRatio, p.adjust,Description)%>%
  filter(ID %in% ID_unique) %>% magrittr::use_series(Cluster)

cluster3 <- cluster2
cluster3@termsim <- cluster2@termsim[Desc_found,Desc_found]
cluster3@compareClusterResult <- cluster2@compareClusterResult %>% 
  filter(ID %in% ID_unique)
pdf("Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PrevsPost_Gruffi_GOenrich.pdf")
emapplot(cluster3, showCategory = 30) + scale_fill_manual(values = c(pre="dark red", post="dark blue"))
dev.off()

saveRDS(cluster2, 'Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PRE_vs_POST_GO.Rds')
saveRDS(cluster3, 'Gruffi_revision/Int_neuronal_Kan25_Kan26_Vel5_Vel6/PRE_vs_POST_GO_selection.Rds')

PRE_int_new@result %>% filter(Description == "response to unfolded protein")
cluster2@compareClusterResult %>% filter(ID == "GO:0006986")


openxlsx::write.xlsx(list(PRE = PRE_int@result %>%dplyr::select(ID, p.adjust),
                          POST = POST_int@result %>%dplyr::select(ID, p.adjust),
                          PRE_ONLY = cluster3@compareClusterResult %>% filter(Cluster == "pre")%>%dplyr::select(ID, p.adjust),
                          POST_ONLY = cluster3@compareClusterResult %>% filter(Cluster == "post")%>%dplyr::select(ID, p.adjust)), 
                     '/Users/Oliver.Eichmueller/Desktop/GO_export.xlsx', overwrite = T)
