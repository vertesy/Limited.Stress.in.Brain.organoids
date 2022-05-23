######################################################################
# Fig.1.panels.2022.02.01.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.1.panels.2022.02.01.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
require(pheatmap)

# Setup ------------------------
Dir.Fig.1 <- p0(OutDirOrig, "/Fig.1/")
create_set_OutDir(Dir.Fig.1)


# Metadata ------------------------
# Parameters ------------------------

clUMAP(ident.Fig.1)
clUMAP(ident.Fig.1, aspect.ratio = .6, suffix = '.6')

# QC ------------------------
{
  "Kanton s3 iN's"
  ls.df.CellFractions <- scBarplot.CellFractions(group.by = 'integrated_snn_res.0.3', fill.by = 'ShortNames'
                                              , downsample = F, return_table = T)

  df.CellFractions <- ls.df.CellFractions$values

  (iN.Kan.s3 <- percentage_formatter(df.CellFractions["14","Kan.s3"] / sum(df.CellFractions["14",])))
  (Fraction.in.small <- percentage_formatter( df.CellFractions["14","Kan.s3"] / sum(df.CellFractions[c("14", "5"), ])))

}


# Fig.1A ------------------------------------------------------------------------------------------------------------------------------------------------
"Fig.1A is a table"
dir.create("A.is.a.table")

ident.Fig.1 <- 'integrated_snn_res.0.3.Manual.short'
# Fig.1B ------------------------------------------------------------------------------------------------------------------------------------------------
{
  "Fig.1B"
  create_set_OutDir(Dir.Fig.1, "/B/")
    ident.Fig.1 <- c('integrated_snn_res.0.3')
    ident.Fig.1.Labeled <- ppp(ident.Fig.1, 'Manual.short')

    clUMAP(ident.Fig.1)
    clUMAP(ident.Fig.1.Labeled)
    clUMAP(ident.Fig.1.Labeled, label = T, axes = F)
    clUMAP(ident.Fig.1.Labeled, label = T, axes = F, PNG = F)
    clUMAP(ident.Fig.1.Labeled, label = F, legend = F, axes = F, suffix = "nolabel")
    clUMAP('integrated_snn_res.0.3.Manual.short', highlight.clusters = 'Unclear Scattered')

}; oo(); create_set_OutDir(Dir.Fig.1)

# Fig.1C-D ------------------------------------------------------------------------------------------------------------------------------------------------
{
  "Fig.1C-D"
  "Current version done by O.Eichmueller"

  create_set_OutDir(Dir.Fig.1, "/C-D/")

  Genes.panel.C <- c('PDK1', 'ENO1', 'LDHA', 'VEGFA', 'DDIT3', 'XBP1', 'DDIT4', 'P4HB')

  multiFeaturePlot.A4(list.of.genes = Genes.panel.C, raster = T, suffix = 'raster'
                      , nr.Col = 4, nr.Row = 2, layout = F, w = 16, h = 6
                      , format = 'pdf' )

  multiFeaturePlot.A4(list.of.genes = Genes.panel.C, raster = T, suffix = 'raster.small'
                      , nr.Col = 4, nr.Row = 2, layout = F, w = 10, h = 4
                      , format = 'pdf')

  # noraster
  multiFeaturePlot.A4(list.of.genes = Genes.panel.C, raster = F, cex = 0.025, suffix = "no.raster.0.05"
                      , nr.Col = 4, nr.Row = 2, layout = F, w = 16, h = 8
                      , format = 'png')

  multiFeaturePlot.A4(list.of.genes = Genes.panel.C, raster = F, cex = 0.005, suffix = "no.raster.0.01"
                      , nr.Col = 4, nr.Row = 2, layout = F, w = 16, h = 8
                      , format = 'jpg')


}; create_set_OutDir(Dir.Fig.1)

# Fig.1E,F,G ------------------------------------------------------------------------------------------------------------------------------------------------
dir.create("E.is.STRING.network")
dir.create("F.is.WGNCA.score")
dir.create("G.is.WGNCA.GO")


# Fig.1H ------------------------------------------------------------------------------------------------------------------------------------------------
{
  create_set_OutDir(Dir.Fig.1, "/H/")

  # Which cells are in the stress clusters
  clUMAP('integrated_snn_res.0.3')
  clUMAP('integrated_snn_res.0.3.ordered')


  # combined.obj$'is.in.Stress.Cls' <- combined.obj[['integrated_snn_res.0.3']][,1] %in% c(6,7)

  fill.by1 = 'project'
  fill.by2 = 'ShortNames'

  # combined.obj$'is.HQ' <- (combined.obj$'nFeature_RNA' > 1000)

  # Stress.per.dataset (small groups) ------------------------------------------------------------------------
  {

    ct.cells.sn <-
      combined.obj@meta.data %>%
      group_by( !!as.name(fill.by2) ) %>%
      summarize(n = n())


    ct.LQ.sn <-
      combined.obj@meta.data %>%
      filter( is.HQ == FALSE ) %>%
      group_by( !!as.name(fill.by2) ) %>%
      summarize(n = n())

    (ct.stressed.PROG <-
        combined.obj@meta.data %>%
        filter( integrated_snn_res.0.3 == 6 ) %>%
        group_by( !!as.name(fill.by2) ) %>%
        summarize(n = n()))


    (ct.stressed.NEU <-
        combined.obj@meta.data %>%
        filter( integrated_snn_res.0.3 == 8 ) %>%
        group_by( !!as.name(fill.by2) ) %>%
        summarize(n = n()))

    (ct.stressed.total <-
        combined.obj@meta.data %>%
        filter( is.in.Stress.Cls == TRUE ) %>%
        group_by( !!as.name(fill.by2) ) %>%
        summarize(n = n()))

    ct.discard.sn <-
      combined.obj@meta.data %>%
      filter( is.HQ == FALSE | is.in.Stress.Cls == TRUE ) %>%
      group_by( !!as.name(fill.by2) ) %>%
      summarize(n = n())


    sum(ct.LQ.sn[,2])
    sum(ct.stressed.NEU[,2])
    sum(ct.stressed.PROG[,2])
    sum(ct.discard.sn[,2])

    # clUMAP('integrated_snn_res.0.3')



    {
      left_join(x = ct.LQ.sn, y = ct.stressed.NEU, ct.stressed.PROG, by='ShortNames' )

      Fig.1H.heatmap.data.stress.per.cl <- plyr::join_all(list(ct.LQ.sn, ct.stressed.NEU, ct.stressed.PROG, ct.stressed.total, ct.discard.sn)
                                , by='ShortNames', type='left')  %>% FirstCol2RowNames() %>% t() %>% replace_na(replace = 0)

      write.simple.tsv(Fig.1H.heatmap.data.stress.per.cl)

      rownames(Fig.1H.heatmap.data.stress.per.cl) <- c("Low Qual.", "Str. Neur.", "Str. Prog.", "Str. Total", "LQ & Stress")


      per.sample.ph.rawC.sn <- pheatmap(Fig.1H.heatmap.data.stress.per.cl, angle_col = 45, cluster_cols = F, cluster_rows = F)
      wplot_save_pheatmap(per.sample.ph.rawC.sn, suffix = 'CLUSTER', width = 15, height = 2)

      Fig.1H.Fraction.per.sample.htmp.sn <- 100 * colDivide(mat = Fig.1H.heatmap.data.stress.per.cl, vec = ct.cells.sn$n)
      write.simple.tsv(Fig.1H.Fraction.per.sample.htmp.sn)

      per.sample.ph.Frac.sn <- pheatmap(Fig.1H.Fraction.per.sample.htmp.sn
                             , angle_col = 45, cluster_cols = F, cluster_rows = F
                             , display_numbers = T, number_format = "%.0f", fontsize_number = 15)
      wplot_save_pheatmap(per.sample.ph.Frac.sn, width = 15, height = 3, suffix = 'CLUSTER')



      Fig.1.H.per.sample.ph.Frac.sn.gap <- pheatmap(Fig.1H.Fraction.per.sample.htmp.sn
                                                    , gaps_col = c(3,6,9,11,32)
                                                    , gaps_row = c(3,4)
                                                    , angle_col = 45, cluster_cols = F, cluster_rows = F
                                                    , display_numbers = T, number_format = "%.0f", fontsize_number = 15
                                                    , cellwidth = 25)
      wplot_save_pheatmap(Fig.1.H.per.sample.ph.Frac.sn.gap, width = 15, height = 2.5, suffix = 'CLUSTER')

      # ------------------------------------------------------

      per.sample.ph.Frac.sn.only <- pheatmap(Fig.1H.Fraction.per.sample.htmp.sn[1:2, ]
                                             , angle_col = 45, cluster_cols = F, cluster_rows = F
                                             , display_numbers = T, number_format = "%.0f", fontsize_number = 15)
      wplot_save_pheatmap(per.sample.ph.Frac.sn.only, width = 15, height = 1.5, suffix = 'CLUSTER')


      per.sample.ph.Frac.sn.gap.only <- pheatmap(Fig.1H.Fraction.per.sample.htmp.sn[1:2, ]
                                                 , gaps_col = c(3,6,9,11,32)
                                                 , angle_col = 45, cluster_cols = F, cluster_rows = F
                                                 , display_numbers = T, number_format = "%.0f", fontsize_number = 15)
      wplot_save_pheatmap(per.sample.ph.Frac.sn.gap.only, width = 15, height = 1.5, suffix = 'CLUSTER')


      per.sample.ph.Frac.3.sn.gap.only <- pheatmap(Fig.1H.Fraction.per.sample.htmp.sn[1:3, ]
                                                   , gaps_col = c(3,6,9,11,32)
                                                   , angle_col = 45, cluster_cols = F, cluster_rows = F
                                                   , display_numbers = T, number_format = "%.0f", fontsize_number = 15
                                                   , cellwidth = 25)
      wplot_save_pheatmap(per.sample.ph.Frac.3.sn.gap.only, width = 15, height = 1.75, suffix = 'CLUSTER.25')
      # isave.RDS(obj = Fig.1.H.ph.Frac.3.sn.gap.only)

    } #Plots
  } # Nyuszika groups version


}; create_set_OutDir(Dir.Fig.1)
