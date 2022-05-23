######################################################################
# Fig.4CD.downampled.v4.R
######################################################################
# source('~/GitHub/Projects/SEO/Analysis/Fig.4CD.downampled.v4.R')
# stop();
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)


# Setup and Parameters ------------------------
source('~/GitHub/Packages/Seurat.pipeline/Load.packages.local.R')
source('~/GitHub/Projects/SEO/Analysis/Gene.Lists.SEO.2022.R')
source('~/GitHub/Projects/SEO/Analysis/Parameters.SEO.2022.R')
set.seed(p$'seed')


p$'min.granule.size' = 50
p$'snn_res' <- c(.1, .2, .3, .4, .5)

# Setup ------------------------
OutDir <- OutDirOrig <- "~/Dropbox (VBC)/Abel.IMBA/AnalysisD/SEO/Fig.4CD.downsample.v4/"
setup_MarkdownReports(OutDir = OutDir, scriptname = 'Fig.4CD.downampled.v4.R')
md.LogSettingsFromList(p)

combined.obj <- read_rds('~/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration.small.Fig.4.CD/obj.Fig.4.CD__2021.01.31_16.03.Rds.gz')


{
  DefaultAssay(combined.obj) <- 'integrated'
  combined.obj <- DietSeurat(combined.obj)
  tic(); combined.obj <- FindVariableFeatures(combined.obj); toc()
  tic(); combined.obj <- ScaleData(combined.obj, verbose = T); toc()
  tic(); combined.obj <- RunPCA(combined.obj, npcs = p$'n.PC', verbose = T); toc()
  tic(); combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions = 3:2, reduction = "umap"); toc()
  combined.obj <- FlipReductionCoordinates(flip = 'x', dim = 2, reduction = 'umap', obj = combined.obj)

  tic(); combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
  tic(); combined.obj <- FindClusters(combined.obj, resolution = p$'snn_res'); toc()

  clUMAP(obj = combined.obj, prefix = "Fig.4C", suffix = 'before', 'integrated_snn_res.0.3', aspect.ratio = 1, h = 6)
  clUMAP("Phase", legend = T, label =F)
}

cells.per.project <- table(combined.obj$'project')
qbarplot(cells.per.project)

{
  ManualNames.DS.before.0.4 <- c(
    '0' = "UL-EN",
    '1' = "IN",
    '2' = "DL-EN",
    '3' = "Progenitors",
    '4' = "Stressed Neurons",
    '5' = "LQ-cells",
    '6' = "oRG",
    '7' = "Immature neurons",
    '8' = "Glia",
    '9' = "Dividing, S-phase",
    '10' = "Stressed Progenitors",
    '11' = "Dividing, G2M-phase",
    '12' = "IPC",
    '13' = "MSN")
  combined.obj$'ManualNames.0.4.before' <-  translate(vec = as.character(combined.obj$'integrated_snn_res.0.4')
                                                      , oldvalues = names(ManualNames.DS.before.0.4)
                                                      , newvalues = ManualNames.DS.before.0.4)

  clUMAP(obj = combined.obj, prefix = "Fig.4C", suffix = 'before.names', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, axes = F)
  clUMAP(obj = combined.obj, prefix = "Fig.4C", suffix = 'before.names.pure', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, label = F, legend = F, axes = F)
}

# gruffi ------------------------


{
  # devtools::load_all("~/gruffi/")
  devtools::load_all(path = "~/GitHub/Packages/gruffi/")
  # ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

  go1 <- "GO:0006096" # Glyco
  go2 <- "GO:0034976" # ER
  go3 <- "GO:0042063" # Glia

}

DefaultAssay(combined.obj) <- 'integrated'

clUMAP('integrated_snn_res.0.3', suffix = 'Before')
clUMAP('integrated_snn_res.0.5', suffix = 'Before')


combined.obj <- aut.res.clustering(combined.obj, assay = 'integrated', upper.res = 20) # in the following "integrated_snn_res.25.reassigned"
granule.res.4.gruffi.raw <- "integrated_snn_res.10"

plot.clust.size.distr(category = granule.res.4.gruffi.raw, vline = p$'min.granule.size', filtercol = T)

combined.obj <- reassign.small.clusters(combined.obj, ident = granule.res.4.gruffi.raw, cell.num = p$'min.granule.size')

granule.res.4.gruffi <- sppp(granule.res.4.gruffi.raw, "reassigned")
plot.clust.size.distr(category = granule.res.4.gruffi, suffix = "reclassified", vline = p$'min.granule.size')



create_set_OutDir(OutDirOrig)

# ------------------------
{
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go1,
                                      new_GO_term_computation = T, clustering = granule.res.4.gruffi)
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go2,
                                      new_GO_term_computation = T, clustering = granule.res.4.gruffi)
  combined.obj <- GO_score_evaluation(obj = combined.obj, GO_term = go3,
                                      new_GO_term_computation = T, clustering = granule.res.4.gruffi)

  PlotGoTermScores(GO = go1, title_ = 'glycolytic process', only.draw.plot = T) # glycolytic process
  PlotGoTermScores(GO = go2, title_ = 'response to endoplasmic reticulum stress', only.draw.plot = T) #
  PlotGoTermScores(GO = go3, title_ = 'gliogenesis', only.draw.plot = T)

  i1 <- kppu(granule.res.4.gruffi, 'cl.av', go1)
  i2 <- kppu(granule.res.4.gruffi, 'cl.av', go2)
  i3 <- kppu(granule.res.4.gruffi, 'cl.av', go3)
}
# ------------------------
combined.obj <- Shiny.GO.thresh(obj = combined.obj,
                                stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3,
                                plot.cluster.shiny = "orig.ident")
clUMAP('is.Stressed', suffix = 'w.Gruffi')
isave.RDS(combined.obj, inOutDir = T, suffix = 'SEO.small.w.Gruffi' )


# combined.obj <- read_rds('/Users/abel.vertesy/Dropbox (VBC)/Abel.IMBA/AnalysisD/SEO/Downsample/Analysis.2022.02.14.Fig.4CD.downsample.v2/combined.obj_SEO.small.w.Gruffi_2022.02.17_18.27.Rds.gz')
combined.obj.before <- combined.obj



# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
{
  combined.obj$'is.HQ' <- combined.obj$'nFeature_RNA' > 1000

  cellIDs.keep <- which_names(!combined.obj$'is.Stressed' & combined.obj$'is.HQ')
  combined.obj <- subset(x = combined.obj, cells = cellIDs.keep)
  clUMAP('is.Stressed', suffix = 'removed')

  p$'n.PC' = 50
  tic(); combined.obj <- FindVariableFeatures(combined.obj); toc()
  tic(); combined.obj <- ScaleData(combined.obj, verbose = T); toc()
  tic(); combined.obj <- RunPCA(combined.obj, npcs = p$'n.PC', verbose = T); toc()
  tic(); combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions = 3:2, reduction = "umap"); toc()
  tic(); combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:p$'n.PC'); toc()
  tic(); combined.obj <- FindClusters(combined.obj, resolution = p$'snn_res'); toc()
  plot3D.umap(category = 'integrated_snn_res.0.3')
  qMarkerCheck.BrainOrg()
}


{

  clUMAP(obj = combined.obj, prefix = "Fig.4D", suffix = 'after', 'integrated_snn_res.0.5', aspect.ratio = 1, h = 6, axes = F)
  clUMAP(obj = combined.obj, prefix = "Fig.4D", suffix = 'after.names', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, axes = F)
  clUMAP(obj = combined.obj, prefix = "Fig.4D", suffix = 'after.names.pure', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, label = F, legend = F, axes = F)

  clUMAP(obj = combined.obj.before, prefix = "Fig.4C", suffix = 'before.names', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, axes = F)
  clUMAP(obj = combined.obj.before, prefix = "Fig.4C", suffix = 'before.names.pure', ident = 'ManualNames.0.4.before', aspect.ratio = 1, h = 6, PNG = F, label = F, legend = F, axes = F)

}
  scBarplot.CellsPerCluster(ident = 'ManualNames.0.4.before')
  scBarplot.CellsPerCluster(ident = 'ManualNames.0.4.before', obj = combined.obj.before)



{
  combined.obj
  Ident.for.DEG = 'new.cluster.ids'
  clUMAP(ident = Ident.for.DEG)
  Idents(combined.obj)  <- Ident.for.DEG
  df.markers <- FindAllMarkers(object = combined.obj)
  write.simple.tsv(df.markers)
  Seurat.utils::barplot.cells.per.cluster(ident = Ident.for.DEG)
  qUMAP('nFeature_RNA')
  qUMAP('POLR2A')
}

combined.obj.after <- combined.obj
isave.RDS(combined.obj.after, inOutDir = T, suffix = 'SEO.small.AFTER.Gruffi' )



clUMAP(obj = combined.obj.before, prefix = "Fig.4C", ident = 'integrated_snn_res.0.4', aspect.ratio = 1, h = 6, PNG = F, axes = F)
