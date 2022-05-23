# Eze et al stressed cells
# Oliver Eichmueller
# Wed May 11 16:41:12 2022 ------------------------------

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer);
require(qgraph); require(enrichplot); require(readr)

OutDir <- 'Gruffi_revision/Eze_2021/analysis/'
dir.create(OutDir)

# load datasets ----------------------------------------------------------------

sc.data <- read_table('Gruffi_revision/Eze_2021/downloaded/counts_exprMatrix.tsv.gz')
sc.meta <- read_table('Gruffi_revision/Eze_2021/downloaded/meta.tsv')
row.names(sc.meta) = sc.meta$Cell
sc.umap <- read_table('Gruffi_revision/Eze_2021/downloaded/UMAP.coords.tsv.gz', col_names = F)

sc.obj <- CreateSeuratObject(counts = sc.data[,-1], project = "Eze_2021", row.names = sc.data$gene)

sc.obj@meta.data <- 
  sc.obj@meta.data %>%
  mutate(Cell = row.names(.)) %>%
  left_join(sc.meta %>% select(-nCount_RNA, - nFeature_RNA), by = "Cell") %>%
  data.frame(row.names = .$Cell)
sc.obj@meta.data$percent.mt
row.names(sc.obj)[stringr::str_detect(row.names(sc.obj), pattern = "^MT\\-")]
sc.obj[["percent.mito"]] <- PercentageFeatureSet(sc.obj, pattern = "^MT\\-")
sc.obj = subset(sc.obj, nFeature_RNA>500&percent.mito<15)

sc.obj = NormalizeData(object = sc.obj, normalization.method = "LogNormalize", scale.factor = 10000)

RowsNA<-names(which(rowSums(is.infinite(sc.obj@assays$RNA@counts))>0))

sc.obj <- sc.obj[setdiff(row.names(sc.obj), RowsNA),]

sc.obj = FindVariableFeatures(object = sc.obj)

sc.obj = ScaleData(sc.obj, verbose = T)

sc.obj = RunPCA(sc.obj, npcs = 50, verbose = T)

sc.obj = SetupReductionsNtoKdimensions(obj = sc.obj, nPCs = 25, dimensions=3:2, reduction="umap")

sc.obj@misc$reductions.backup$umap_paper <- CreateDimReducObject(as.matrix(data.frame(UMAP_1 = sc.umap$X2, UMAP_2 = sc.umap$X3, row.names = sc.umap$X1)))
sc.obj@reductions$umap <- sc.obj@misc$reductions.backup$umap2d

sc.obj@reductions$umap <- sc.obj@misc$reductions.backup$umap_paper

DimPlot(sc.obj, group.by = "Cluster", raster = T, reduction = "pca")
DimPlot(sc.obj, group.by = "Cluster", raster = F, reduction = "umap") + NoLegend()
clUMAP(ident = "Cluster", obj = sc.obj)

sc.obj = FindNeighbors(sc.obj, reduction = "pca", dims = 1:25)
sc.obj = FindClusters(sc.obj, resolution = c(0.1,0.2,0.3,0.5, 0.6))

clUMAP("RNA_snn_res.0.1", obj = sc.obj, save.plot = F, label = T)

# Run Predictions of cell types ------------------------------------------------

sc.obj <- predict_non_neuronal_Seurat(sc.obj, model.files = '~/R/Git_Projects/OrganoidClassifier/Classifiers/xgboost_Cao_model.Rds')
sc.obj$predicted.cell.type.nonN <- sc.obj$predicted.cell.type

sc.obj <- predict_neuronal_Seurat(sc.obj, model.files = '~/R/Git_Projects/OrganoidClassifier/Classifiers/xgb_brain_model_final.Rds')
sc.obj$predicted.cell.type.N <- sc.obj$predicted.cell.type

sc.obj <- predict_neuronal_Seurat(sc.obj, model.files = '~/R/Git_Projects/OrganoidClassifier/Classifiers/Polioudakis_fine_labels.Rds')
sc.obj$predicted.cell.type.N.fine <- sc.obj$predicted.cell.type

clUMAP("predicted.cell.type.nonN", obj = sc.obj, save.plot = F)
clUMAP("predicted.cell.type.N", obj = sc.obj, save.plot = F, legend = T, cols = pals::tol(8))

qUMAP("DDIT4", obj = sc.obj, save.plot = F)

# Run Gruffi on Cakir ----------------------------------------------------

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

sc.obj <- aut.res.clustering(obj = sc.obj)

granule.res.4.gruffi <- sc.obj@misc$gruffi$'optimal.granule.res'	

sc.obj <- reassign.small.clusters(sc.obj, ident = granule.res.4.gruffi)

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

# Glycolytic process	GO:0006096
sc.obj <- GO_score_evaluation(obj = sc.obj, GO_term = go1, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# ER stress 	GO:0034976
sc.obj <- GO_score_evaluation(obj = sc.obj, GO_term = go2, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Gliogenesis		GO:0042063
sc.obj <- GO_score_evaluation(obj = sc.obj, GO_term = go3, save.UMAP = F, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)


# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

# Call Shiny app
sc.obj <- Shiny.GO.thresh(obj = sc.obj,
                          stress.ident1 = i1,
                          stress.ident2 = i2,
                          notstress.ident3 = i3,
                          plot.cluster.shiny = "orig.ident")



pl1 <- Seurat.utils::clUMAP(obj = sc.obj, 'is.Stressed',label = F, legend = T, save.plot = F, cols = c("dark grey","dark red"))+
  NoAxes()
pl2 <- Seurat.utils::clUMAP(ident = "Age_.Carnegie_Stage.", label = F, legend = T, cols = pals::tol(7), label.cex = 5,
                            obj = sc.obj, save.plot = F, title = "Age (Carnegie Stage)") +NoAxes() 

png(filename = "Gruffi_revision/Eze_2021/analysis/summary_v1.png", width = 1500,height = 1500, res = 200)
cowplot::plot_grid(pl2,pl1, ncol = 2)
dev.off()
rm(pl1,pl2,pl3)
png(filename = "Gruffi_revision/Eze_2021/analysis/GO_Scores.png", width = 2500,height = 1000, res = 200)
FeaturePlot(sc.obj, features = c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red"), ncol = 3) & 
  NoAxes() & coord_fixed(0.6)
dev.off()


saveRDS(sc.obj, file = "Gruffi_revision/Eze_2021/analysis/Eze_object_proc_220511.Rds")
# sc.obj <- readRDS(file = "Gruffi_revision/Eze_2021/analysis/Eze_object_proc_220511.Rds")

sc.obj.pre  <- readRDS('/Users/Oliver.Eichmueller/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration/Before.filtering/combined.obj_w.Gruffi_2022.02.11_11.25.Rds.gz')

plot_df_int_all <- sc.obj.pre@meta.data %>%
  group_by(integrated_snn_res.120.reassigned) %>%
  summarize(across(.cols = c("integrated_snn_res.120.reassigned_cl.av_GO:0006096", 
                             "integrated_snn_res.120.reassigned_cl.av_GO:0034976", 
                             "integrated_snn_res.120.reassigned_cl.av_GO:0042063",
                             "is.Stressed"), 
                   .fns = function(x) x=unique(x)))

pdf("Gruffi_revision/Eze_2021/analysis/Scores_perGranule_org.pdf", width = 10, height = 10)
ggplot(plot_df_int_all,aes(x = as.numeric(as.vector(`integrated_snn_res.120.reassigned_cl.av_GO:0006096`)),
             y = as.numeric(as.vector(`integrated_snn_res.120.reassigned_cl.av_GO:0034976`)), 
             color = is.Stressed, size = as.numeric(as.vector(`integrated_snn_res.120.reassigned_cl.av_GO:0042063`)) )) + 
  geom_point(alpha = 0.5) + theme_minimal() + coord_fixed(1) + scale_color_manual(values = c("grey", "dark red")) +
  xlab("GO:0006096 Score") + ylab("GO:0034976 Score") + 
  scale_size(name = "GO:0042063 score", breaks = c(seq(0,10,2)/10), range = c(0,10)) +
  xlim(-2,4) + ylim(-1.5,2.5)
dev.off()
pdf("Gruffi_revision/Eze_2021/analysis/Scores_perGranule_fetal.pdf", width = 10, height = 10)
sc.obj@meta.data %>%
  group_by(RNA_snn_res.25.reassigned) %>%
  summarize(across(.cols = c("RNA_snn_res.25.reassigned_cl.av_GO:0006096", 
                             "RNA_snn_res.25.reassigned_cl.av_GO:0034976", 
                             "RNA_snn_res.25.reassigned_cl.av_GO:0042063",
                             "is.Stressed"), 
                   .fns = function(x) x=unique(x))) %>%
  ggplot(aes(x = as.numeric(as.vector(`RNA_snn_res.25.reassigned_cl.av_GO:0006096`)),
             y = as.numeric(as.vector(`RNA_snn_res.25.reassigned_cl.av_GO:0034976`)), 
             color = is.Stressed, size = as.numeric(as.vector(`RNA_snn_res.25.reassigned_cl.av_GO:0042063`)) )) + 
  geom_point(alpha = 0.5) + theme_minimal() + coord_fixed(1) + scale_color_manual(values = c("grey", "dark red"))+
  xlab("GO:0006096 Score") + ylab("GO:0034976 Score") + 
  scale_size(name = "GO:0042063 score", breaks = c(seq(0,10,2)/10), range = c(0,10)) +
  xlim(-2,4) + ylim(-1.5,2.5)
dev.off()

pdf("Gruffi_revision/Eze_2021/analysis/Venn_comparison_GO_found.pdf", width = 20, height = 5)
cowplot::plot_grid(
ggvenn::ggvenn(list(fetal.GO.0006096 = sc.obj@misc$GO$GO.0006096,
                    org.GO.0006096 = sc.obj.pre@misc$GO$GO.0006096), auto_scale = T, 
               set_name_size = 4, fill_color = c(pals::tol(2))),
ggvenn::ggvenn(list(fetal.GO.0034976 = sc.obj@misc$GO$GO.0034976,
                    org.GO.0034976 =sc.obj.pre@misc$GO$GO.0034976), auto_scale = T, 
               set_name_size = 4, fill_color = c(pals::tol(2))),
ggvenn::ggvenn(list(fetal.GO.0042063 = sc.obj@misc$GO$GO.0042063,
                    org.GO.0042063 =sc.obj.pre@misc$GO$GO.0042063), auto_scale = T, 
               set_name_size = 4, fill_color = c(pals::tol(2))),
ncol = 3)
dev.off()
