######################################################################
# Fig.3.panels.2021.12.01.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.3.panels.2021.12.07.R')

# Functions ------------------------
# Setup ------------------------
Dir.Fig.3 <- p0(OutDirOrig, "Fig.3/")
create_set_OutDir(Dir.Fig.3)


# Fig.3A ------------------------------------------------------------------------------------------------------------------------------------------------
# Mitochondrial.Reads.in.Stressed.Cells ---------------------------------------------------------------------------
{

  create_set_OutDir(Dir.Fig.3, "/A/")

  {
    clUMAP('integrated_snn_res.0.1')
    clUMAP('integrated_snn_res.0.1', highlight.clusters = 1:2)

    combined.obj$is.Neural <- combined.obj$integrated_snn_res.0.1 %!in% 1:2
    clUMAP('is.Neural')

    combined.obj$'Stress.and.Fate' <- translate(vec = combined.obj$is.Neural
                                                , oldvalues = c(T,F)
                                                , newvalues = c('Neurons', 'Prog.'))


    idx.Str.Neurons = which(combined.obj$is.Neural & combined.obj$is.Stressed)
    idx.Str.Progenitors = which(!combined.obj$is.Neural & combined.obj$is.Stressed)

    combined.obj$'Stress.and.Fate'[idx.Str.Neurons] <- 'Str. Neurons'
    combined.obj$'Stress.and.Fate'[idx.Str.Progenitors] <- 'Str. Prog.'

    clUMAP('Stress.and.Fate', palette = 'alphabet')

  }


  # GO:0030433 Score Autophagy of the mitochondrion
  # GO:0000422 Score Ubiquitin dep. ERAD pathway
  # combined.obj$

  ReadFractions.in.Stressed.Cells <-
    combined.obj@meta.data[,c("Stress.and.Fate", "percent.ribo", "percent.mito", 'Score.GO.0030433', 'Score.GO.0000422', 'Score.GO.0006914')] %>%
    tibble() %>%
    mutate('percent.ribo' = 100 * percent.ribo) %>%
    mutate('percent.mito' = 100 * percent.mito) %>%
    dplyr::rename(Stress = Stress.and.Fate)

  Mitochondrial.Reads.in.Stressed.Cells <- ReadFractions.in.Stressed.Cells[, c(1,3)]
  Ribosomal.mRNA.Reads.in.Stressed.Cells <- ReadFractions.in.Stressed.Cells[, c(1,2)]
  Mitophagy.in.Stressed.Cells <- ReadFractions.in.Stressed.Cells[, c(1,4)]
  ERAD.in.Stressed.Cells <- ReadFractions.in.Stressed.Cells[, c(1,5)]
  Autophagy.in.Stressed.Cells <- ReadFractions.in.Stressed.Cells[, c(1,6)]

  colz <- gg_color_hue(2)[c(2,1,1,2)]
  (pn.A1 <- qviolin(Mitochondrial.Reads.in.Stressed.Cells
                   , plotname = "MT. mRNA"
                   , xlab.angle = 45
                   , ylab = "% transcriptome", ylim = c(0,8)
                   , draw_quantiles = 0.5
                   , hide.legend = T
                   , stat.test = F #, stat.label.y = 'top'
                   , palette = colz))

  (pn.A2 <- qviolin(Ribosomal.mRNA.Reads.in.Stressed.Cells
                    , plotname = "Ribosomal mRNA"
                    , xlab = "", xlab.angle = 45
                    , ylab = "% transcriptome"
                    , draw_quantiles = 0.5
                    , stat.test = F
                    , hide.legend = T
                    , palette = colz))


  (pn.A3 <- qviolin(Mitophagy.in.Stressed.Cells
                    , plotname = "Autophagy of the mitochondrion"
                    , xlab = "", xlab.angle = 45
                    , ylab = "GO:0030433 Score", ylim = c(-0.2,0.4)
                    , draw_quantiles = 0.5
                    , stat.test = F
                    , hide.legend = T
                    , palette = colz))

  (pn.A3.sup <- qviolin(Autophagy.in.Stressed.Cells, suffix = "Fig.EV3"
                    , plotname = "Autophagy"
                    , xlab = "", xlab.angle = 45
                    , ylab = "GO:0006914 Score"
                    , draw_quantiles = 0.5
                    , stat.test = F
                    , hide.legend = T
                    , palette = colz
                    , w = 4, h = 4))


  (pn.A4 <- qviolin(ERAD.in.Stressed.Cells
                    , plotname = "Score Ubiquitin dep. ERAD pathway"
                    , xlab = "", xlab.angle = 45
                    , ylab = "GO:0000422 Score", ylim = c(-0.2,0.4)
                    , draw_quantiles = 0.5
                    , stat.test = F
                    , hide.legend = T
                    , palette = colz))

  Fig.3A.panels <- list(
    'A' = pn.A1,
    '1' = ggplot() + theme_void(),
    '2' = ggplot() + theme_void(),
    'B' = pn.A2,
    '3' = ggplot() + theme_void(),
    '4' = ggplot() + theme_void(),
    'C' = pn.A3,
    '5' = ggplot() + theme_void(),
    '6' = ggplot() + theme_void(),
    'D' = pn.A4,
    '7' = ggplot() + theme_void()
  )

  qA4_grid_plot(Fig.3A.panels, nrow = 4, ncol = 3, plotname = sppp('Fig.3A.panels', 4, 'by', 3), scale = 1, labels = NULL)
  qA4_grid_plot(Fig.3A.panels, nrow = 4, ncol = 3, plotname = sppp('Fig.3A.panels', 4, 'by', 3), scale = 1, labels = NULL, extension = 'pdf')
  qA4_grid_plot(Fig.3A.panels, nrow = 4, ncol = 3, plotname = sppp('Fig.3A.panels', 4, 'by', 3), scale = 1.1, labels = NULL, suffix = '110pc')
  qA4_grid_plot(Fig.3A.panels, nrow = 4, ncol = 3, plotname = sppp('Fig.3A.panels', 4, 'by', 3), scale = 1.1, labels = NULL, suffix = '110pc', extension = 'pdf')

}

#  -------------------------------------------------------------


#  -------------------------------------------------------------
{
  (Table.S5X.Stress.Comparison <-
     ReadFractions.in.Stressed.Cells %>%
     group_by(Stress) %>%
     summarise_at(vars( "percent.ribo", "percent.mito", "Score.GO.0030433", "Score.GO.0000422", 'Score.GO.0006914'), mean) %>%
     FirstCol2RowNames.as.df())

  write.simple.tsv(Table.S5X.Stress.Comparison)


  {

    irequire(rstatix)

    # Pairwise comparisons ------------------------------------
    Mitochondrial.Reads.in.Stressed.Cells %>%
      kruskal_test(percent.mito ~ Stress)

    Ribosomal.mRNA.Reads.in.Stressed.Cells %>%
      kruskal_test(percent.ribo ~ Stress)

    Mitophagy.in.Stressed.Cells %>%
      kruskal_test(Score.GO.0030433 ~ Stress)

    ERAD.in.Stressed.Cells %>%
      kruskal_test(Score.GO.0000422 ~ Stress)


    # Effect sizes ------------------------------------
    Mitochondrial.Reads.in.Stressed.Cells %>%
      kruskal_effsize(percent.mito ~ Stress)

    Ribosomal.mRNA.Reads.in.Stressed.Cells %>%
      kruskal_effsize(percent.ribo ~ Stress)

    Mitophagy.in.Stressed.Cells %>%
      kruskal_effsize(Score.GO.0030433 ~ Stress)

    ERAD.in.Stressed.Cells %>%
      kruskal_effsize(Score.GO.0000422 ~ Stress)


    # Pairwise comparisons ------------------------------------
    pwc <- Mitochondrial.Reads.in.Stressed.Cells %>%
      dunn_test(percent.mito ~ Stress, p.adjust.method = "bonferroni")
    view(pwc)

    pwc2 <- Ribosomal.mRNA.Reads.in.Stressed.Cells %>%
      dunn_test(percent.ribo ~ Stress, p.adjust.method = "bonferroni")
    view(pwc2)

    pwc3 <- Mitophagy.in.Stressed.Cells %>%
      dunn_test(Score.GO.0030433 ~ Stress, p.adjust.method = "bonferroni")
    view(pwc3)


    pwc4 <- ERAD.in.Stressed.Cells %>%
      dunn_test(Score.GO.0000422 ~ Stress, p.adjust.method = "bonferroni")
    view(pwc4)

  }

}


# Fig.3B ------------------------------------------------------------------------------------------------------------------------------------------------
{
  create_set_OutDir(Dir.Fig.3, "/3/")
  # BiocManager::install("progeny")
  require(progeny)
  require(pheatmap)

  {
    LabeledIdent <- NewIdent <- combined.obj$'integrated_snn_res.0.3.Manual.short'
    stressed.cellIDs <- which_names(combined.obj$is.Stressed)
    stressed.NeuronIDs <- which_names(LabeledIdent == "Stressed Neurons")
    stressed.ProgenitorsIDs <- which_names(LabeledIdent == "Stressed Prog.")

    True.Stressed.Progenitors <- intersect(stressed.cellIDs, which_names(NewIdent == "Stressed Prog."))
    True.Stressed.Neurons <- intersect(stressed.cellIDs, which_names(NewIdent == "Stressed Neurons"))

    NewIdent[NewIdent == "Stressed Prog."] <- "Stressed Prog. (Clustering)"
    NewIdent[True.Stressed.Progenitors] <- "Stressed Prog. (Gruffi)"

    NewIdent[NewIdent == "Stressed Neurons"] <- "Stressed Neurons (Clustering)"
    NewIdent[True.Stressed.Neurons] <- "Stressed Neurons (Gruffi)"

    Clusters.and.Stress.Split.by.Gruffi <- table(NewIdent)
    qbarplot(Clusters.and.Stress.Split.by.Gruffi, hline = 3103, xlab.angle = 45, subtitle = "Cut at 2%", xlab = '', ylab='Cells')
    combined.obj$'integrated_snn_res.0.3.Manual.Stress.Cl.vs.Gruffi' <- NewIdent
    # combined.obj.sub$'integrated_snn_res.0.3.Manual.Stress.Cl.vs.Gruffi' <- NewIdent[colnames(combined.obj.sub)]

  }


  res.X <- 0.3
  create_set_OutDir(Dir.Fig.3, "/B.", res.X, "proper2")
  identity.4.progeny <- 'integrated_snn_res.0.3.Manual.Stress.Cl.vs.Gruffi'

  ngenes.progeny <- 200
  clUMAP(identity.4.progeny)


  combined.obj.sub <- subsetSeuObj(obj = combined.obj, fraction_ = 0.33)
  # isave.RDS(combined.obj.sub, inOutDir = T)

  Idents(combined.obj.sub) <- identity.4.progeny
  CellsClusters <- data.frame(Cell = names(Idents(combined.obj.sub)),
                              CellType = as.character(Idents(combined.obj.sub)),
                              stringsAsFactors = FALSE)

  tic(); combined.obj.sub <- progeny(combined.obj.sub, scale = FALSE, organism = "Human", verbose = T
                              , top = ngenes.progeny, perm = 1, return_assay = TRUE); toc(); say()

  ## We can now directly apply Seurat functions in our Progeny scores.
  ## For instance, we scale the pathway activity scores.
  combined.obj.sub <- Seurat::ScaleData(combined.obj.sub, assay = "progeny" )

  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <-
    as.data.frame(t(GetAssayData(combined.obj.sub, slot = "scale.data",
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
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


  progeny_scores_per_cluster <- t(summarized_progeny_scores_df)
  Table.SX.progeny.full <- progeny_scores_per_cluster
  write.simple.tsv(Table.SX.progeny.full)


  # Supp Figure ----------------------------------------------------------------------------------
  hm.size = 8
  {
    myColor2 = matlabColors.pheatmap(progeny_scores_per_cluster)
    Fig.S3B.progeny.all.cl = pheatmap(progeny_scores_per_cluster, fontsize = 14,
                                  fontsize_row = 10,
                                  clustering_method = "ward.D2",
                                  cutree_cols = 3,
                                  cutree_rows = 3,
                                  color = myColor2,
                                  # main = p0("PROGENy (", ngenes.progeny, " genes)")
                                  main = p0("Pathway activity scores")
                                  , angle_col = 45
                                  ,  border_color = NA
    ); Fig.S3B.progeny.all.cl

    wplot_save_pheatmap(x = Fig.S3B.progeny.all.cl
                        , suffix = kpp(identity.4.progeny, "genes",ngenes.progeny)
                        , width = hm.size, height = hm.size)
  }

  # Main Figure ----------------------------------------------------------------------------------
  {
    TwoPC <- round(ncol(combined.obj) *0.02)

    Clz <- table(combined.obj@meta.data[,identity.4.progeny])
    LargeClusters <- which_names(namedVec = Clz > TwoPC)
    scBarplot.CellsPerCluster(ident =  identity.4.progeny) + geom_hline(yintercept = TwoPC)

    progeny_scores_per_cluster <- progeny_scores_per_cluster[, LargeClusters]

    Fig.3B.progeny.200 = pheatmap(progeny_scores_per_cluster, fontsize = 14,
                                  fontsize_row = 10,
                                  clustering_method = "ward.D2",
                                  cutree_cols = 3,
                                  cutree_rows = 3,
                                  color = myColor2,
                                  # main = p0("PROGENy (", ngenes.progeny, " genes)")
                                  main = p0("Pathway activity scores")
                                  , angle_col = 45
                                  ,  border_color = NA
    )

    wplot_save_pheatmap(x = Fig.3B.progeny.200
                        , suffix = kpp(identity.4.progeny, "genes",ngenes.progeny)
                        , width = hm.size, height = hm.size)
  }
}



# Fig.3C ------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(kpps(Dir.Fig.3,"C.is.by.Oliver"))


# Fig.3D-H ------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(kpps(Dir.Fig.3,"D-H.is.by.Julia"))

# End ------------------------------------------------------------------------------------------------------------------------------------------------




create_set_Original_OutDir()
isave.RDS(combined.obj, inOutDir = T)
