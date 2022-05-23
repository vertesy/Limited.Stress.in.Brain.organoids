rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

library(Seurat)
library(ggplot2)
library(Seurat.utils)

OutDir <- OutDirOrig <- "./analysis_fetal_organoids_downsampled_RDS/" 
setwd(OutDir)

combined.obj <- readRDS("combined.obj.processed.invivo.SEO.downsampled.24211.RDS")

#-------- FIGURE 3 - SUBFIGURE D -------#
c <- c(gg_color_hue(length(unique(combined.obj$predicted.cell.type))-1), "darkgrey")
c <- c(c[2],c[1],c[3:length(c)])
clUMAP(obj = combined.obj, ident = "predicted.cell.type", save.plot = F, shuffle = TRUE,
       label = F, legend = T, col = c, pt.size = .5, order = F) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank()) + NoLegend() +
  ggsave(paste0("UMAP.predicted.cell.type.png"), width = 10, height = 7)

#-------- FIGURE 3 - SUBFIGURE E -------#
FeaturePlot(obj = combined.obj, features = "Score.GO.0006096", min.cutoff = "q01", max.cutoff = "q99", pt.size = .5)+ NoAxes()+ ggtitle(element_blank()) +
  ggsave(paste0("UMAP_Score.GO.0006096.png"), width = 10, height = 7)

FeaturePlot(obj = combined.obj, features = "Score.GO.0034976", min.cutoff = "q01", max.cutoff = "q99", pt.size = .5)+ NoAxes()+ ggtitle(element_blank()) +
  ggsave(paste0("UMAP_Score.GO.0034976.png"), width = 10, height = 7)

FeaturePlot(obj = combined.obj, features = "Score.GO.0042063", min.cutoff = "q01", max.cutoff = "q99", pt.size = .5)+ NoAxes()+ ggtitle(element_blank()) +
  ggsave(paste0("UMAP_Score.GO.0042063.png"), width = 10, height = 7)

#-------- FIGURE 3 - SUBFIGURE F -------#
clUMAP(obj = combined.obj, ident = "is.Stressed", save.plot = F, label = F, legend = F, splitby = "orig.ident", pt.size = .5) +
  theme(text = element_text(size=20)) + NoAxes() + ggtitle(element_blank())+ scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) +
  ggsave(paste0("UMAP_Stress_cluster_split.png"), width = 18, height = 7)

#-------- FIGURE 3 - SUBFIGURE G -------#
# stress assignment barplot cell numbers
ggplot(combined.obj@meta.data, aes(x=combined.obj$orig.ident, fill=is.Stressed, stat = "count")) +
  scale_fill_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  geom_bar() + xlab("") +
  ggsave(paste0("Barplot_stressed.cells.png"), width = 3, height = 7)
