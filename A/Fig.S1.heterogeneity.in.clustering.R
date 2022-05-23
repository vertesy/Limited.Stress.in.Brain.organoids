######################################################################
# Fig.S1.heterogeneity.in.clustering.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/supplementary/Fig.S1.heterogeneity.in.clustering.R')

# Setup ------------------------

# OutDirOrig = OutDir = "~/Dropbox (VBC)/Group Folder Knoblich/papers_in_progress/2021 Stress paper/2022.01.28.EMBO.Abel/Figures.01.28/Fig.S4/"
# setup_MarkdownReports(OutDir = OutDir, scriptname = "Fig.S4.heterogeneity.in.clustering.R")


# Setup ------------------------
Dir.Fig.S1 <- p0(OutDirOrig, "/Fig.S1/")
create_set_OutDir(Dir.Fig.S1)

b.save.wplots = F
# (I) Boundaries ------------------------
{

  "(I) Boundaries"
  Cl.Res.Compared = GetClusteringRuns()[1:4]
  pt.size = 0.025
  panels.Fig.S1C <- list(
    'A' = clUMAP(ident = Cl.Res.Compared[1], pt.size = pt.size, save.plot = F) + NoAxes(),
    'B' = clUMAP(ident = Cl.Res.Compared[2], pt.size = pt.size) + NoAxes(),
    'C' = clUMAP(ident = Cl.Res.Compared[3], pt.size = pt.size) + NoAxes(),
    'D' = clUMAP(ident = Cl.Res.Compared[4], pt.size = pt.size) + NoAxes()
  )

  # qA4_grid_plot(plot_list = panels.Fig.S1C, plotname = 'Fig.S4.cluster.boundaries'
  #               , nrow = 2, ncol = 2
  #               , h = hA4/2, w = wA4/1)

  qA4_grid_plot(plot_list = panels.Fig.S1C, plotname = 'Fig.S1C.cluster.boundaries.linear'
                , nrow = 1, ncol = 4
                , h = hA4/3, w = 2*wA4)

}



# (II) Heterogeneity ------------------------
{
  "(II) Heterogeneity"
  pt.size = 0.025
  panels.Fig.S1D <- list(
    'A' = qUMAP(feature = 'percent.mito', cells = StressedNeuronIDs, pt.size = pt.size) + NoAxes(),
    'B' = qUMAP(feature = 'percent.ribo', cells = StressedNeuronIDs, pt.size = pt.size) + NoAxes(),
    'C' = qUMAP(feature = 'GAPDH', cells = StressedNeuronIDs, pt.size = pt.size) + NoAxes(),
    'D' = qUMAP(feature = 'nFeature_RNA', cells = StressedNeuronIDs, pt.size = pt.size) + NoAxes()
  )

  # qA4_grid_plot(plot_list = panels.Fig.S1D, plotname = 'Fig.S4.heterogeneity.in.clustering.block'
  #               , nrow = 2, ncol = 2
  #               , h = hA4/2, w = wA4/1)

  qA4_grid_plot(plot_list = panels.Fig.S1D, plotname = 'Fig.S1D.heterogeneity.in.clustering.linear'
                , nrow = 1, ncol = 4
                , h = hA4/3, w = 2*wA4)

  oo()
}



# End ------------------------
b.save.wplots = T
create_set_OutDir(OutDirOrig)


