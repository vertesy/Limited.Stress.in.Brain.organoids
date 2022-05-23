# Qian et al stressed cells
# Oliver Eichmueller
# Tue May 10 19:03:49 2022 ------------------------------

require(Stringendo);require(CodeAndRoll2); require(ReadWriter); require(MarkdownHelpers);
require(MarkdownReports); require(ggExpress); require(Seurat.utils); require(clusterProfiler);
require(pdist); require(raster); require(rgl); require(tidyr); require(corrplot); require(pheatmap)

# devtools::load_all("/Users/Oliver.Eichmueller/Library/R/x86_64/4.1/library/gruffi")
source("/Users/Oliver.Eichmueller/Documents/GitHub/gruffi/R/gruffi.R")
library(gruffi)
require(Seurat); require(dplyr); require(ggplot2); require(ggpubr);  require(RColorBrewer);
require(qgraph); require(enrichplot); require(readr)

OutDir <- 'Gruffi_revision/Qian_2020/analysis/'
dir.create(OutDir)

# load datasets ----------------------------------------------------------------

sc.data <- read_table('Gruffi_revision/Qian_2020/GSM4094681_C1_150d_fseq12BC77_S2.deg.txt')
colnames(sc.data)[3]
row.names(sc.data) <- sc.data$GENE

sc.obj <- CreateSeuratObject(counts = sc.data[,-1], project = "Qian_2020", row.names = sc.data$GENE)

sc.obj[["percent.mito"]] <- PercentageFeatureSet(sc.obj, pattern = "^MT\\-")

sc.obj = subset(sc.obj, nFeature_RNA>500&percent.mito<15)

sc.obj = NormalizeData(object = sc.obj, normalization.method = "LogNormalize", scale.factor = 10000)

sc.obj = FindVariableFeatures(object = sc.obj)

sc.obj = ScaleData(sc.obj, verbose = T)

sc.obj = RunPCA(sc.obj, npcs = 50, verbose = T)

sc.obj = SetupReductionsNtoKdimensions(obj = sc.obj, nPCs = 25, dimensions=3:2, reduction="umap")
sc.obj = RunTSNE(obj = sc.obj, dims = 1:25)

sc.obj = FindNeighbors(sc.obj, reduction = "pca", dims = 1:25)
sc.obj = FindClusters(sc.obj, resolution = c(0.1,0.2,0.3,0.5, 0.6))

clUMAP("RNA_snn_res.0.1", obj = sc.obj, save.plot = F)

# Run Predictions of cell types ------------------------------------------------

sc.obj <- predict_non_neuronal_Seurat(sc.obj, model.files = '~/R/Git_Projects/OrganoidClassifier/Classifiers/xgboost_Cao_model.Rds')
sc.obj$predicted.cell.type.nonN <- sc.obj$predicted.cell.type

sc.obj <- predict_neuronal_Seurat(sc.obj, model.files = '~/R/Git_Projects/OrganoidClassifier/Classifiers/xgb_brain_model_final.Rds')
sc.obj$predicted.cell.type.N <- sc.obj$predicted.cell.type

clUMAP("predicted.cell.type.nonN", obj = sc.obj, save.plot = F, reduction = "tsne")
clUMAP("predicted.cell.type.N", obj = sc.obj, save.plot = F, reduction = "tsne")

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


pl1 <- Seurat.utils::clUMAP(obj = sc.obj, 'is.Stressed', label =F, save.plot = F, cols = c("dark grey","dark red"))+NoAxes() +NoLegend()
pl2 <- Seurat.utils::clUMAP(ident = "RNA_snn_res.0.6", label = T, cols = pals::tol(12), label.cex = 5,
                            obj = sc.obj, save.plot = F) +NoAxes() +NoLegend()

pl3 <- Seurat.utils::clUMAP(ident = "predicted.cell.type.N", label = T, cols = pals::tol(12), label.cex = 5,
                            obj = sc.obj, save.plot = F) +NoAxes() +NoLegend()
pl4 <- clUMAP(ident = "RNA_snn_res.0.6", obj = sc.obj, save.plot = F, aspect.ratio = 1,
              cols = pals::tol(12), reduction = "tsne", title = "tSNE\nresolution as paper") + NoAxes()
qUMAP(feature = "RELN", obj = sc.obj, save.plot = F, 
      cols = c(alpha("grey", alpha = 0.5), "red")) + NoAxes()

png(filename = "Gruffi_revision/Qian_2020/analysis/summary_v1.png", width = 1500,height = 1500, res = 200)
cowplot::plot_grid(pl2,pl4, pl3,pl1, ncol = 2)
dev.off()

saveRDS(sc.obj, file = "Gruffi_revision/Qian_2020/analysis/Qian_object_proc_220510.Rds")
# sc.obj <- readRDS(file = "Gruffi_revision/Qian_2020/analysis/Qian_object_proc_220510.Rds")


png("Gruffi_revision/Qian_2020/analysis/Scores_Qian.png", width = 1000, height = 1000)
FeaturePlot(sc.obj, features = c("Score.GO.0006096", "Score.GO.0034976", "Score.GO.0042063"), 
            min.cutoff = 'q10', max.cutoff = 'q90', cols = c(alpha("grey", 0.5), "red")) & NoAxes() & coord_fixed(0.6)
dev.off()

# Check progeny for stressed cells ---------------------------------------------


CellsClusters <- data.frame(Cell = row.names(sc.obj@meta.data), 
                            CellType = ifelse(sc.obj$is.Stressed == TRUE, "Stressed",
                                              sc.obj$RNA_snn_res.0.1),
                            stringsAsFactors = FALSE)


## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
sc.obj <- progeny(sc.obj, scale=FALSE, organism="Human", top=500, perm=1, 
                  return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
sc.obj <- Seurat::ScaleData(sc.obj, assay = "progeny") 

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


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

pdf('Gruffi_revision/Qian_2020/analysis/PROGENy_Qian.pdf')
pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
         fontsize_row = 10, 
         color=myColor, 
         scale = "row",
         main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 0,  border_color = NA)
dev.off()



# Subset data and re-cluster ---------------------------------------------------

sc.obj.subset <- subset(sc.obj, is.Stressed==FALSE&nFeature_RNA>1000)


sc.obj.subset = NormalizeData(object = sc.obj.subset, 
                              normalization.method = "LogNormalize", scale.factor = 10000)

sc.obj.subset = FindVariableFeatures(object = sc.obj.subset)

sc.obj.subset = ScaleData(sc.obj.subset, verbose = T)

sc.obj.subset = RunPCA(sc.obj.subset, npcs = 50, verbose = T)

sc.obj.subset = SetupReductionsNtoKdimensions(obj = sc.obj.subset, 
                                              nPCs = 50, dimensions=3:2, reduction="umap")

sc.obj.subset = RunTSNE(sc.obj.subset, dims = 1:50)

sc.obj.subset = FindNeighbors(sc.obj.subset, reduction = "pca", dims = 1:50)
sc.obj.subset = FindClusters(sc.obj.subset, resolution = c(0.1,0.2,0.3,0.5))

clUMAP("RNA_snn_res.0.1", obj = sc.obj.subset, save.plot = F)

png('Gruffi_revision/Qian_2020/analysis/predicted.cell.type_POSTGruffi.png')
clUMAP("predicted.cell.type.N", obj = sc.obj.subset, save.plot = F) +NoAxes() +NoLegend()
dev.off()

saveRDS(sc.obj.subset, file = "Gruffi_revision/Qian_2020/analysis/Qian_object_POSTGruffi_proc_220510.Rds")

qUMAP(feature = "DLX1", obj = sc.obj.subset, save.plot = F)
clUMAP(ident = "predicted.cell.type.nonN", obj = sc.obj.subset, save.plot = F)
qUMAP(feature = "DDIT4", obj = sc.obj, save.plot = F)
