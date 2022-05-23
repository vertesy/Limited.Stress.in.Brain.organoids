######################################################################
# Fig.2.panels.2021.12.01.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.2.panels.2021.12.06.R')

# Functions ------------------------
# Setup ------------------------
Dir.Fig.2 <- p0(OutDirOrig, "Fig.2/")
create_set_OutDir(Dir.Fig.2)


# Fig.2A ------------------------------------------------------------------------------------------------------------------------------------------------
{
  "Fig.2A scores"
  create_set_OutDir(Dir.Fig.2, "/A/")

  # Parameters ------------------------
  ident.Fig.2 <- ident.Fig.1.Labeled
  qTHR = 0.09
  UseZscore = T


  # Plot ------------------------
  Scores.V2 <- NULL
  Scores.V2 <- list(
    "Response to ER stress" <- qUMAP(feature = 'Score.GO.0034976', sub = 'GO:0034976-score', title = "Response to ER stress") + NoAxes(),
    "Glycolytic process" <- qUMAP(feature = 'Score.GO.0006096', sub = 'GO:0006096-score', title = "Glycolytic process") + NoAxes(),
    "Neurogenesis" <- qUMAP(feature = 'Score.GO.0022008', sub = 'GO:0022008-score', title = "Neurogenesis") + NoAxes(),
    "Gliogenesis" <- qUMAP(feature = 'Score.GO.0042063', sub = 'GO:0042063-score', title = "Gliogenesis") + NoAxes()
  )

  plg <- plot_grid(plotlist = Scores.V2, ncol = 4) # , rel_heights = c(1.5,1)
  fname <- ppp('Fig.2.A',qTHR, flag.nameiftrue(UseZscore), ident.Fig.2)
  save_plot(filename = ppp(fname, 'png'), plot = plg, base_width = 20, base_height = 4, limitsize = FALSE)
  save_plot(filename = ppp(fname, 'pdf'), plot = plg, base_width = 20, base_height = 4, limitsize = FALSE)
  oo()
}

# Fig.2B ------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(kpps(Dir.Fig.2,"B.is.the.cartoon"))

# Fig.2C-F ------------------------------------------------------------------------------------------------------------------------------------------------
{
  "Fig.1 C-F - Steps of Gruffi"
  create_set_OutDir(Dir.Fig.2, "/C_F/")
  l(unique(combined.obj$integrated_snn_res.120.reassigned))

  combined.obj$integrated_snn_res.120.reassigned <- as.numeric(as.character(combined.obj$integrated_snn_res.120.reassigned))

  combined.obj$'Stressed.Lab' <- translate(combined.obj$'is.Stressed', oldvalues = 0:1, newvalues = c("Normal", "Stressed"))


  ScoreExample <- 'integrated_snn_res.120.reassigned_cl.av_GO:0006096'
  Granule.Av.Score.0006096 <- as.numeric(as.character(combined.obj@meta.data[,ScoreExample]))
  combined.obj$Granule.Av.Score.0006096 <- Granule.Av.Score.0006096
  clUMAP('integrated_snn_res.120.reassigned', axes = T, MaxCategThrHP = Inf, label = F, legend = F)
  gr.distr <- table(combined.obj$integrated_snn_res.120.reassigned)


  Fig.1E <- list(
    'C' = clUMAP('integrated_snn_res.120.reassigned', axes = F, MaxCategThrHP = Inf, label = F, legend = F
                 , title = "Granules", save.plot = F) + NoAxes() + labs(caption = paste(l(gr.distr), "granules, median size:", median(gr.distr))),
    'D' = qUMAP(feature = 'Granule.Av.Score.0006096', axes = F, cols = c('lightgrey','red')
                , title = 'Granule Average Score',  save.plot = F) + NoLegend() + labs(caption = "Glycolysis (GO:0006096)"),
    # 'E' = qhistogram(Granule.Av.Score.0006096, vline = combined.obj@misc$gruffi$thresh.stress.ident1, filtercol = T , add = NULL
    #                  , plotname = "Threshold Estimation", xlab = 'Granule Av. Score', ylab = 'Nr. Granules'
    #                  , palette_use = rev(gg_color_hue(2)), save = F) + labs(caption = "Glycolysis (GO:0006096)"),
    'E' = ggplot() + theme_void(),
    'F' = clUMAP('Stressed.Lab', axes = F, label = F, legend = F, cols = rev(gg_color_hue(2))
                 , title =  "Stress Idenitfication", save.plot = F) + NoAxes() + labs(caption = "Combining Glycolysis, ER-stress and Gliogenesis"),

    # For legend
    'D2' = qUMAP(feature = 'Granule.Av.Score.0006096', axes = F, cols = c('lightgrey','red'), title = 'Granule Average Score',  save.plot = F) + labs(caption = "Glycolysis (GO:0006096)"),
    'F2' = clUMAP('Stressed.Lab', axes = F, label = F, legend = T, cols = rev(gg_color_hue(2)), title =  "Stress Idenitfication", save.plot = F)  + NoAxes()

  )
  nr = 6

  q32vA4_grid_plot(Fig.1E, nrow = nr, ncol = 4, scale = 2)
  q32vA4_grid_plot(Fig.1E, nrow = nr, ncol = 4, scale = 3, suffix = '3', extension = 'pdf')
  q32vA4_grid_plot(Fig.1E, nrow = nr, ncol = 4, scale = 4, suffix = '4', extension = 'png')

}
oo()



# Fig.2E 3D ------------------------------------------------------------------------------------------------------------------------------------------------

{
  "Fig.2E 3D"
  Granule.Scorenames.3D <- c(
    'integrated_snn_res.120.reassigned_cl.av_GO:0034976',
    'integrated_snn_res.120.reassigned_cl.av_GO:0042063',
    'integrated_snn_res.120.reassigned_cl.av_GO:0006096',
    'is.Stressed'
  )
  Granule.Scores.Meta <- as.tibble(combined.obj@meta.data[, Granule.Scorenames.3D])
  colnames(Granule.Scores.Meta) <- c(
    'GO.0034976',
    'GO.0042063',
    'GO.0006096',
    'Stress'
  )

  library(plotly)
  {
    # 'GO:0071456' 	cellular response to hypoxia
    # 'GO:0006096' 	glycolytic process
    # 'GO:0034976' 	response to endoplasmic reticulum stress

    Title.Fig.2E <- kollapse("Gruffi Scores 3D")
    Fig.2E <- plotly::plot_ly(Granule.Scores.Meta, color = ~Stress, colors = rev(gg_color_hue(2))
                              , size = 1
                              , z = ~GO.0006096
                              # , x = ~Score.GO.0022008
                              , y = ~GO.0034976
                              , x = ~GO.0042063
    ) # , colors = c('#BF382A', '#0C4B8E')
    Fig.2E <- Fig.2E %>% add_markers()
    Fig.2E <- Fig.2E %>% plotly::layout(title = Title.Fig.2E, scene = list(
      zaxis = list(title = 'Glycolytic process'), #  - GO.0006096
      yaxis = list(title = 'Response to ER stress'), #  - GO.0034976
      xaxis = list(title = 'Gliogenesis') #  - GO.0042063
      , aspectmode='cube'
    ))
    Fig.2E
    SavePlotlyAsHtml(Fig.2E, category. = Title.Fig.2E)
  }
}




# Fig. Cluster annotation  ------------------------------------------------------------------------------------------------------------------------------------------------
