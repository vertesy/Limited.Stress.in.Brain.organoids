rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

library(Seurat)
library(ggplot2)
library(Seurat.utils)

OutDir <- OutDirOrig <- "./Vel.Org7.d90.ImpV_RDS"
setwd(OutDirOrig)

combined.obj <- readRDS("./Vel.Org7.d90.ImpV_RDS/Vel.Org7.d90.ImpV.processed.RDS")
combined.obj.filtered <- readRDS("./Vel.Org7.d90.ImpV_RDS/Vel.Org7.d90.ImpV.stressfiltered.RDS")

#-------- EXTENDED VIEW FIGURE 2 - SUBFIGURE A -------#
FeaturePlot(obj =  combined.obj, features = "Score.GO.0006096", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_Score.GO.0006096.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "Score.GO.0034976", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_Score.GO.0034976.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "Score.GO.0042063", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_Score.GO.0042063.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "Score.GO.0022008", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_Score.GO.0022008.png"), width = 10, height = 7)

#-------- EXTENDED VIEW FIGURE 2 - SUBFIGURE B -------#
DimPlot(obj = combined.obj, group.by = "is.Stressed.res.1.reassigned", label = F, pt.size =2) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank()) + scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + 
  ggsave(paste0("UMAP.Stress_cluster.1.reassigned.png"), width = 10, height = 7)

DimPlot(obj = combined.obj, group.by = "is.Stressed.res.6.reassigned", label = F, pt.size =2) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank())+ scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + 
  ggsave(paste0("UMAP.Stress_cluster.6.reassigned.png"), width = 10, height = 7)

DimPlot(obj = combined.obj, group.by = "is.Stressed.res.10.reassigned", label = F, pt.size =2) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank())+ scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + 
  ggsave(paste0("UMAP.Stress_cluster.10.reassigned.png"), width = 10, height = 7)

DimPlot(obj = combined.obj, group.by = "is.Stressed.res.20.reassigned", label = F, pt.size =2) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank())+ scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + 
  ggsave(paste0("UMAP.Stress_cluster.20.reassigned.png"), width = 10, height = 7)

#-------- EXTENDED VIEW FIGURE 2 - SUBFIGURE C -------#
FeaturePlot(obj =  combined.obj, features = "TOP2A", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_TOP2A.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "HOPX", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_HOPX.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "EOMES", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_EOMES.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "KAZN", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_KAZN.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj, features = "SATB2", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_SATB2.png"), width = 10, height = 7)

#-------- EXTENDED VIEW FIGURE 2 - SUBFIGURE D -------#
FeaturePlot(obj =  combined.obj.filtered, features = "TOP2A", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_TOP2A_after_filtering.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj.filtered, features = "HOPX", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_HOPX_after_filtering.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj.filtered, features = "EOMES", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_EOMES_after_filtering.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj.filtered, features = "KAZN", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_KAZN_after_filtering.png"), width = 10, height = 7)

FeaturePlot(obj =  combined.obj.filtered, features = "SATB2", min.cutoff = "q01", max.cutoff = "q99", pt.size =2)+ NoAxes()+ ggtitle(element_blank()) +
  theme(legend.text=element_text(size=20)) + ggsave(paste0("UMAP_SATB2_after_filtering.png"), width = 10, height = 7)
