# Test Gruffi on non-neuronal cell datasets

# # Install CRAN & Bioconductor dependencies (to ensure everything is up to date)
# install.packages('pdist')
# install.packages('raster')
# install.packages('rgl')
# 
# install.packages('BiocManager')
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# 
# # Install custom dependencies --------------------------------------------------
# install.packages('devtools')
# devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
# devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
# devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
# devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
# devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
# devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
# devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)


# Install gruffi
# devtools::install_github(repo = "jn-goe/gruffi", force = TRUE)


require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer)

OutDir <- 'Gruffi_revision/'
dir.create(OutDir)
# Load datasets ----------------------------------------------------------------
Vel5 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Vel.Org5.d90.ImpV.201029_non_neuronal.Rds')
Vel6 <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Vel.Org6.d90.ImpV.201029_non_neuronal.Rds')
coul <- colorRampPalette(RColorBrewer::brewer.pal(8, "PiYG"))(25)
# Run Gruffi for Vel5 ----------------------------------------------------------

Vel5 <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = Vel5, nPCs = 50, dimensions=3:2, reduction="umap")
labels = c("immature EN", "DL-EN", "UL-EN", "Glia", "Retina/Adrenocortical", "IN", "Immune Cells")
Vel5$RNA_snn_res.0.1_named <- factor(Vel5$RNA_snn_res.0.1, labels = labels)
colors <- pals::tol(7)
names(colors) <- labels

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

Vel5 <- aut.res.clustering(obj = Vel5)

granule.res.4.gruffi <- Vel5@misc$gruffi$'optimal.granule.res'	

Vel5 <- reassign.small.clusters(Vel5, ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

require(MarkdownHelpers)
# Glycolytic process	GO:0006096
Vel5 <- GO_score_evaluation(obj = Vel5, GO_term = go1, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
Vel5 <- GO_score_evaluation(obj = Vel5, GO_term = go2, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
Vel5 <- GO_score_evaluation(obj = Vel5, GO_term = go3, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)


# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
Vel5 <- Shiny.GO.thresh(obj = Vel5,
                        stress.ident1 = i1,
                        stress.ident2 = i2,
                        notstress.ident3 = i3,
                        plot.cluster.shiny = "orig.ident")

saveRDS(Vel5, 'Gruffi_revision/Int_Vel5_Vel6/Vel5_stress_proc.Rds')

cellIDs.keep <- which_names(!Vel5$'is.Stressed')
subset.obj <- subset(x = Vel5, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj, save.plot = F)
subset.obj <- FindVariableFeatures(subset.obj)
subset.obj <- ScaleData(subset.obj)
subset.obj <- RunPCA(subset.obj, features = VariableFeatures(subset.obj))
subset.obj <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = subset.obj, nPCs = 50, dimensions=3:2, reduction="umap")
Seurat.utils::clUMAP('predicted.clusters', label = F, obj = subset.obj, save.plot = F)
Seurat.utils::clUMAP('RNA_snn_res.0.1', label = F, obj = subset.obj, save.plot = F)
saveRDS(subset.obj, 'Gruffi_revision/Int_Vel5_Vel6/Vel5_stress_subset.Rds')
subset.obj.Vel5 <- readRDS('Gruffi_revision/Int_Vel5_Vel6/Vel5_stress_subset.Rds')


Cao_ref <- readRDS("/Users/Oliver.Eichmueller/Dropbox (VBC)/scRNA_22/Gruffi_Review/Cao_subset_easy_label.Rds")

feat_found <- intersect(VariableFeatures(Cao_ref), VariableFeatures(Vel5))
VariableFeatures(Vel5)
avg_found_Vel5 <- AverageExpression(Vel5, features = feat_found, group.by = "RNA_snn_res.0.1_named")
avg_found_Cao <- AverageExpression(Cao_ref, features = feat_found, group.by = "label_easy")
colnames(avg_found_Vel5$RNA)[7] <- "Immune Cells "
Cor_cao <- as.matrix(cor(cbind(avg_found_Cao$RNA, avg_found_Vel5$RNA[,])))
pdf("Gruffi_revision/Int_Vel5_Vel6/nonN_Vel5_CorCao_select.pdf", width = 10,height = 10)
pheatmap(Cor_cao[1:length(colnames(avg_found_Cao$RNA)),c("Retina/Adrenocortical", "Immune Cells ")], 
         cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 4, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         cellwidth = 10, cellheight = 10)
dev.off()
pdf("Gruffi_revision/Int_Vel5_Vel6/nonN_Vel5_CorCao_full.pdf", width = 10,height = 10)
pheatmap(Cor_cao[1:length(colnames(avg_found_Cao$RNA)),(length(colnames(avg_found_Cao$RNA))+1):length(colnames(Cor_cao))], 
         cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 4, cutree_cols = 2, width = 2,height = 2, 
         border_color = NA,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         cellwidth = 10, cellheight = 10)
dev.off()

# Run Gruffi for Vel6 ----------------------------------------------------------

Vel6 <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = Vel6, nPCs = 50, dimensions=3:2, reduction="umap")

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

Vel6 <- aut.res.clustering(obj = Vel6)

granule.res.4.gruffi <- Vel6@misc$gruffi$'optimal.granule.res'	

Vel6 <- reassign.small.clusters(Vel6, ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

# Glycolytic process	GO:0006096
Vel6 <- GO_score_evaluation(obj = Vel6, GO_term = go1, save.UMAP = F,
                            new_GO_term_computation = T, 
                            clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
Vel6 <- GO_score_evaluation(obj = Vel6, GO_term = go2, save.UMAP = F, 
                            new_GO_term_computation = T, 
                            clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
Vel6 <- GO_score_evaluation(obj = Vel6, GO_term = go3, save.UMAP = F, 
                            new_GO_term_computation = T, 
                            clustering = granule.res.4.gruffi, plot.each.gene = F)


# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
Vel6 <- Shiny.GO.thresh(obj = Vel6,
                        stress.ident1 = i1,
                        stress.ident2 = i2,
                        notstress.ident3 = i3,
                        plot.cluster.shiny = "orig.ident")

labels = c("immature EN", "UL-EN", "Stressed EN", "Glia", "DL-EN", "Retina/Smooth Muscle", "Glia2")
Vel6$RNA_snn_res.0.1_named <- factor(Vel6$RNA_snn_res.0.1, labels = labels)

Seurat.utils::clUMAP(obj = Vel6, 'RNA_snn_res.0.1_named', label =F, save.plot = F)


feat_found <- intersect(VariableFeatures(Cao_ref), VariableFeatures(Vel6))

avg_found_Vel6 <- AverageExpression(Vel6, features = feat_found, group.by = "RNA_snn_res.0.1_named")
avg_found_CaovsVel6 <- AverageExpression(Cao_ref, features = feat_found, group.by = "label_easy")

Cor_cao <- as.matrix(cor(cbind(avg_found_CaovsVel6$RNA, avg_found_Vel6$RNA[,])))
pdf("Gruffi_revision/Int_Vel5_Vel6/nonN_Vel6_CorCao_select.pdf", width = 10,height = 10)
pheatmap(Cor_cao[1:length(colnames(avg_found_CaovsVel6$RNA)),c("Retina/Smooth Muscle", "Glia")], 
         cluster_rows = T, cluster_cols = F, display_numbers = F, color = coul, 
         cutree_rows = 4, cutree_cols = 2, width = 2,height = 2, show_rownames = T,show_colnames = T,
         border_color = NA,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         cellwidth = 10, cellheight = 10)
dev.off()
pdf("Gruffi_revision/Int_Vel5_Vel6/nonN_Vel6_CorCao_full.pdf", width = 10,height = 10)
pheatmap(Cor_cao[1:length(colnames(avg_found_CaovsVel6$RNA)),length(colnames(avg_found_CaovsVel6$RNA)):length(colnames(Cor_cao))], 
         cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 4, cutree_cols = 2, width = 2,height = 2, show_rownames = T,show_colnames = T,
         border_color = NA,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         cellwidth = 10, cellheight = 10)
dev.off()


cellIDs.keep <- which_names(!Vel6$'is.Stressed')
subset.obj_vel6 <- subset(x = Vel6, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = subset.obj_vel6, save.plot = F)
subset.obj_vel6 <- FindVariableFeatures(subset.obj_vel6)
subset.obj_vel6 <- ScaleData(subset.obj_vel6)
subset.obj_vel6 <- RunPCA(subset.obj_vel6)
subset.obj_vel6 <- Seurat.utils::SetupReductionsNtoKdimensions(
  obj = subset.obj_vel6, nPCs = 50, dimensions=3:2, reduction="umap")

Seurat.utils::clUMAP('predicted.clusters', label = F, obj = subset.obj_vel6, save.plot = F)
Seurat.utils::clUMAP('RNA_snn_res.0.1', label = F, obj = subset.obj_vel6, save.plot = F)



# Test other vis ================
install.packages("qgraph")
library(qgraph)

stdevs_Vel5 <- Vel5@reductions$pca@stdev
names(stdevs_Vel5) <- colnames(Vel5@reductions$pca@cell.embeddings)

dist_pca_Vel5 <- Vel5@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(Vel5@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(RNA_snn_res.0.1_named) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_Vel5[cur_column()]))) %>%
  data.frame(row.names = paste0(.$RNA_snn_res.0.1_named)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()

stdevs_Vel5_int <- Vel5_fromInt@reductions$pca@stdev
names(stdevs_Vel5_int) <- colnames(Vel5_fromInt@reductions$pca@cell.embeddings)

dist_pca_Vel5_int <- Vel5_fromInt@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(Vel5_fromInt@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(RNA_snn_res.0.1_named) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_Vel5_int[cur_column()]))) %>%
  data.frame(row.names = paste0(.$RNA_snn_res.0.1_named)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()

Vel5_pca_raw <- Vel5@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(Vel5@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(RNA_snn_res.0.1) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_Vel5[cur_column()]))) %>%
  data.frame(row.names = paste0("Cl_", .$RNA_snn_res.0.1)) %>% dplyr::select(-1)


qgraph(1/dist_pca_Vel5_int, layout='spring', vsize=5, groups = colnames(dist_pca_Vel5_int),
       colors = colors, details = T
       ,nNodes = 10
       #, graph = "glasso"
       )
pl <- qgraph(1/dist_pca_Vel5, layout='spring', vsize=10, groups = colnames(dist_pca_Vel5),
       colors = colors, details = T, minimum = 0.005, theme = "gray", 
       title = "Pre-integration Distances", title.cex =2,label.cex=1.5, 
       curveAll = F, edge.labels = F, repulsion=1)
pl2 <- qgraph(1/dist_pca_Vel5_int, layout='spring', vsize=10, groups = colnames(dist_pca_Vel5),
       colors = colors, details = T, minimum = 0.005, theme = "gray", 
       title = "Post-integration Distances", title.cex =2,label.cex=1.5, 
       curveAll = F, edge.labels = F, repulsion=1)

cowplot::plot_grid(plot(pl),plot(pl2))
pdf('Gruffi_revision/Int_Vel5_Vel6/preInt_distances_qgraph.pdf', width = 10, height = 5)
plot(pl)
dev.off()

pdf('Gruffi_revision/Int_Vel5_Vel6/postInt_distances_qgraph.pdf', width = 10, height = 5)
plot(pl2)
dev.off()

brewer.pal.info
cols2 <- colorRampPalette(brewer.pal(9, "Greys"))(25)
annot_df <- data.frame(row.names= labels, `Cell Type` = c(rep("Telencephalon",4), "Mis-differentiated", "Telencephalon", "Mis-differentiated"))
pheatmap(round((dist_pca_Vel5/dist_pca_Vel5_int),2)[,c(5,7)], 
         cluster_rows = F, 
         cluster_cols = F, 
         display_numbers = T, color = cols2, 
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         #show_rownames = F, show_colnames = F,
         border_color = NA,
         annotation_row = annot_df, 
         #annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         clustering_distance_rows = "euclidean", na_col = "white",
         #clustering_distance_cols = "euclidean",
         main = "Weighted euclidean distance in top 50 PCs\nVel5 Single", cellwidth = 20, cellheight = 20)



Vel5$RNA_snn_res.0.1_named_stress <- factor(ifelse(as.vector(Vel5$RNA_snn_res.0.1_named) %in% c("Retina/Adrenocortical", "Immune Cells")|
                                              Vel5$is.Stressed==FALSE, 
                                            as.vector(Vel5$RNA_snn_res.0.1_named),
                                            "Stressed"), levels = c(levels(Vel5$RNA_snn_res.0.1_named), "Stressed"))

DimPlot(Vel5, reduction = "umap", group.by = "RNA_snn_res.0.1_named_stress")


dist_pca_Vel5_stress <- Vel5@reductions$pca@cell.embeddings %>%
  as.data.frame() %>%
  mutate(bc = row.names(.)) %>%
  left_join(Vel5@meta.data %>% mutate(bc = row.names(.))) %>%
  dplyr::group_by(RNA_snn_res.0.1_named_stress) %>%
  summarise((across(.cols = 1:50, .fns = function(x) x = mean(x)))) %>%
  mutate(across(.cols = 2:51, .fns = function(x) x = x * as.vector(stdevs_Vel5[cur_column()]))) %>%
  data.frame(row.names = paste0(.$RNA_snn_res.0.1_named_stress)) %>% dplyr::select(-1) %>%
  dist()%>% as.matrix()
pdf('Gruffi_revision/Int_Vel5_Vel6/preInt_distances_stressed_qgraph.pdf', width = 10, height = 5)
qgraph(1/dist_pca_Vel5_stress, layout='spring', vsize=10, groups = colnames(dist_pca_Vel5),
       colors = c(colors, Stressed="black"), minimum = 0.005, theme = "gray", 
       label.cex=1.5, details = F,
       curveAll = F, edge.labels = F, repulsion=1)
dev.off()

pheatmap(dist_pca_Vel5_stress, cluster_rows = T, cluster_cols = T, display_numbers = F, color = coul, 
         cutree_rows = 2, cutree_cols = 2, width = 2,height = 2, 
         #show_rownames = F, show_colnames = F,
         border_color = NA,
         annotation_row = annot_df, annotation_col = annot_df,
         #annotation_colors = list(predicted.clusters = colors),
         #clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         main = "Weighted euclidean distance in top 50 PCs\nVel5 Post-Integration", cellwidth = 20, cellheight = 20)


