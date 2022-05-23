######################################################################
# Fig.4.panels.2021.12.01.R
######################################################################
# source('~/GitHub/Projects/SEO/Figures.Scripts/Fig.4.panels.2021.12.07.R')

# Functions ------------------------
# Setup ------------------------
Dir.Fig.4 <- p0(OutDirOrig, "Fig.4/")
create_set_OutDir(Dir.Fig.4)


# Paramters ------------------------------------------------------------------------------------------------------------------------------------------------

DataDir = '~/Dropbox (VBC)/Group Folder Knoblich/Papers_in_progress/2021 Stress paper/Data/RDS/Organoid.integration/'



# Fig.4A ------------------------------------------------------------------------------------------------------------------------------------------------
# Mitochondrial.Reads.in.Stressed.Cells ---------------------------------------------------------------------------
{

  create_set_OutDir(Dir.Fig.4, "/A/")

  {
    # combined.obj.before <- read_rds(p0(DataDir, 'Before.filtering/combined.obj_w.Gruffi_2022.02.11_11.25.Rds.gz'))
    clUMAP(ident = 'integrated_snn_res.0.3.Manual.short', obj = combined.obj.before, label = F, legend = F, axes = F, prefix = 'Fig.4.A1', suffix = "Before.Shortnames"); say(); oo();
    plot3D.umap(category = 'integrated_snn_res.0.3.Manual.short', obj = combined.obj.before, suffix = "before", AutoAnnotBy = 'integrated_snn_res.0.3.Manual.short')

    # combined.obj.after <- read_rds(p0(DataDir, 'After.filtering/combined.obj_After.Gruffi.DGEA_2022.02.12_18.42.Rds.gz'))
    clUMAP(ident = 'integrated_snn_res.0.3.Manual.short.Before.Gruffi', obj = combined.obj.after, label = F, legend = F, axes = F, prefix = 'Fig.4.A2', suffix = "After.ShortNames"); say(); oo();
    clUMAP(ident = 'integrated_snn_res.0.3', obj = combined.obj.after, label = F, legend = F, axes = F, prefix = 'Fig.S4.A2', suffix = "After.NewClustering"); say(); oo();
    plot3D.umap(category = 'integrated_snn_res.0.3.Manual.short.Before.Gruffi', obj = combined.obj.after, suffix = "after", AutoAnnotBy = 'integrated_snn_res.0.3.Manual.short.Before.Gruffi')

    try.dev.off()
  }
}



  # Fig.4B ------------------------------------------------------------------------------------------------------------------------------------------------

# Fig.4C ------------------------------------------------------------------------------------------------------------------------------------------------
dir.create(kpps(Dir.Fig.4,"C-D by next script"))



# End ------------------------------------------------------------------------------------------------------------------------------------------------




create_set_Original_OutDir()
isave.RDS(combined.obj, inOutDir = T)
