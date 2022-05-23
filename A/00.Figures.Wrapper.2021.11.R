######################################################################
# 000.Figures.Wrapper.2021.11.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/00.Figures.Wrapper.2021.11.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
source('~/GitHub/Packages/Seurat.pipeline/Load.packages.local.R')

# Setup ------------------------

OutDirOrig  = '~/Dropbox (VBC)/Group Folder Knoblich/papers_in_progress/2021 Stress paper/2022.01.28.EMBO.Abel/Figures.02.11/'
setup_MarkdownReports(OutDir = OutDirOrig, scriptname = "000.Figures.Wrapper.2021.11.R")

# Metadata ------------------------
# Parameters ------------------------

# Read In ------------------------
InputFolder  = '~/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration/Before.filtering/'
combined.obj <- read_rds(p0(InputFolder, '/combined.obj_w.Gruffi_2022.02.11_11.25.Rds.gz'))

{
  recall.all.genes()
  recall.parameters(overwrite = T)
  recall.meta.tags.n.datasets()
  set.mm()
  clUMAP('integrated_snn_res.0.3')
  say()
}


range(combined.obj$percent.ribo)
range(combined.obj$percent.mito)
oo()

# QC ------------------------


if (F) source('~/GitHub/Projects/SEO/Analysis/Analysis.2022.02.01.R')

{
  # Fig. 1 ------------------------
  source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.1.panels.2022.02.01.R')
  oo()

  # Fig. 2 ------------------------
  source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.2.panels.2021.12.06.R')

  source('~/GitHub/Projects/SEO/Figures.Scripts/supplementary/Fig.EV2.panels.2022.02.07.R')

  # Fig. 3 ------------------------
  source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.3.panels.2021.12.07.R')


  # Fig. 4 ------------------------
  source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.4.panels.2022.02.13.R')



}




# Supplementary and EV figures ------------------------
{
  "Supplementary and EV figures"
  source('~/GitHub/Projects/SEO/Figures.Scripts/supplementary/Fig.EV1.panels.2022.02.07.R')

  source('~/GitHub/Projects/SEO/Figures.Scripts/supplementary/Fig.S1.heterogeneity.in.clustering.R')

  source('~/GitHub/Projects/SEO/Figures.Scripts/supplementary/Fig.S2.panels.2022.02.07.R')

}
# ------------------------
# ------------------------
