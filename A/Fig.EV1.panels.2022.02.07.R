######################################################################
# Fig.EV1.panels.2022.02.07.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.EV1.panels.2022.02.07.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try(source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
require('MarkdownReportsDev')
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')
# try(source("~/GitHub/Packages/Seurat.multicore/00.Load.Seurat3.Multicore.LOCAL.R"));

# Setup ------------------------
Dir.Fig.EV1 <- p0(OutDirOrig, "/Fig.EV1/")
setup_MarkdownReports(OutDir = Dir.Fig.EV1, scriptname = "Fig.EV1.R")

# Metadata ------------------------
# Parameters ------------------------


# "Fig.EV1.A" ------------------------
qUMAP('nFeature_RNA')

# "Fig.EV1.B" ------------------------
"Oliver"

# "Fig.EV1.C" ------------------------
{
  "Fig.EV1.C"
  create_set_OutDir(Dir.Fig.EV1, "/C/")
  HiCited.ER.stress.genes <- c('HERPUD1', 'AARS', 'CEBPB', 'BBC3', 'ATF3', 'ATF4'
                               , 'ATF6', 'EIF2AK3', 'STC2', 'DNAJB9', 'CREB3', 'HSPA5')

  multiFeaturePlot.A4(list.of.genes = HiCited.ER.stress.genes, raster = T, suffix = 'raster'
                      , nr.Col = 6, nr.Row = 2, layout = F, w = 12, h=4
                      , format = 'png')

  multiFeaturePlot.A4(list.of.genes = HiCited.ER.stress.genes, raster = F
                      , nr.Col = 6, nr.Row = 2, layout = F, w = 24, h=8
                      , format = 'png') # noraster

  multiFeaturePlot.A4(list.of.genes = HiCited.ER.stress.genes, raster = F
                      , nr.Col = 6, nr.Row = 2, layout = F, w = 24, h=8
                      , format = 'jpg')


}; create_set_OutDir(OutDirOrig)



# "Fig.EV1.D" ------------------------
{
  "Fig.EV1.D"
  create_set_OutDir(Dir.Fig.EV1, "/B/")
  # link_google(vector_of_gene_symbols = Other.Glycolytic.genes, suffix = 'Glycolysis', Open = T)
  HiCited.Glycolytic.genes <- c('SLC2A1','HK2', 'PFKFB3', 'GAPDH', 'TOP2A'
                                , 'PGK1', 'TPI1', 'ENO2', 'PKM', 'TOP2B')
  "TOP2A and TOP2B are just placeholders!!"

  multiFeaturePlot.A4(list.of.genes = HiCited.Glycolytic.genes, raster = T, suffix = 'raster'
                      , nr.Col = 5, nr.Row = 2, layout = F, w = 24, h=8
                      , format = 'png')
  # noraster
  multiFeaturePlot.A4(list.of.genes = HiCited.Glycolytic.genes, raster = F, cex = .05
                      , nr.Col = 5, nr.Row = 2, layout = F, w = 24, h=8
                      , format = 'png')

  multiFeaturePlot.A4(list.of.genes = HiCited.Glycolytic.genes, raster = F, cex = .05
                      , nr.Col = 5, nr.Row = 2, layout = F, w = 24, h=8
                      , format = 'jpg')

}; create_set_OutDir(OutDirOrig)


# Fig.EV1.E" ------------------------
{
  qUMAP('VIM', prefix = "Fig.EV1.E")
}


# "Fig.EV1.F-G" ------------------------
{
  "Fig.EV1.F-G"
  create_set_OutDir(Dir.Fig.EV1, "/Fig.EV1/G/")

  scores.Fig.EV1.G <- c('Score.GO.0042594', 'Score.GO.0071456')

  # qUMAP('Score.GO.0042594', title = 'GO-score response to starvation', sub = 'Score.GO.0042594')
  # qUMAP('Score.GO.0071456', title = 'cellular response to hypoxia', sub = 'Score.GO.0071456')
  try(combined.obj <- PlotGoTermScores(GO = "GO:0042594", desc = "GO-score response to starvation", obj = combined.obj, plot.each.gene = F), silent = T)
  try(combined.obj <- PlotGoTermScores(GO = "GO:0071456", desc = "cellular response to hypoxia", obj = combined.obj, plot.each.gene = F), silent = T)

}


# "Fig.EV1.H-J" ------------------------
{
  "Fig.EV1.H-J ShinyGO KEGG pathways"
  create_set_OutDir(Dir.Fig.EV1, "/Fig.EV1/C-E.ShinyGO")
  all.genes.detected <- rownames(combined.obj@assays$RNA)
  write.simple.tsv(all.genes.detected)
  "top.150.coding.2022.02.tsv from string list"

}; create_set_OutDir(OutDirOrig)






# ------------------------
# ------------------------
# ------------------------


