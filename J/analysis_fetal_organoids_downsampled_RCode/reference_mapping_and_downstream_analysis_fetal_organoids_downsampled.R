rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

library(Seurat)

OutDir <- OutDirOrig <- "./analysis_fetal_organoids_downsampled_RDS/"
setwd(OutDirOrig)
 
# Loading and Splitting Downsampled Organoid Data
SEO.obj.ls <- readRDS("./downsampled_organoid_integration_RDS/combined.obj.processed.SEO.downsampled.24211.RDS")
Idents(SEO.obj.ls) <- SEO.obj.ls$orig.ident
SEO.obj.ls <- SplitObject(SEO.obj.ls)
 
# Loading and Splitting Geschwind Fetal Data
geschwind.obj.ls <- readRDS("./Geschwind_premRNA_CCA__processed_210408.Rds")
Idents(geschwind.obj.ls) <- geschwind.obj.ls$cca_id
geschwind.obj.ls <- SplitObject(geschwind.obj.ls)

# Integrating Organoid and Geschwind Fetal Data
ls.Seurat <- c(SEO.obj.ls, geschwind.obj.ls)
for (j in 1:length(ls.Seurat)) {
  DefaultAssay(ls.Seurat[[j]]) <- "RNA"
  ls.Seurat[[j]] <- NormalizeData(object = ls.Seurat[[j]], normalization.method = "LogNormalize", scale.factor = 10000)
  ls.Seurat[[j]] <- FindVariableFeatures(object = ls.Seurat[[j]], mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR', nfeatures =10000)
}

integration.features <- SelectIntegrationFeatures(object.list = ls.Seurat, nfeatures = 3000)
reference_dataset <- (length(SEO.obj.ls)+1):length(ls.Seurat)
anchors <- FindIntegrationAnchors(object.list = ls.Seurat,
                                  anchor.features = integration.features,
                                  reference = reference_dataset)
combined.obj <- IntegrateData(anchorset = anchors, k.weight = 50) #Org21_Velasco only has 68 sample
DefaultAssay(combined.obj) <- "integrated"

# Downstream Analysis per Seurat on Integrated Object
combined.obj <- FindVariableFeatures(object = combined.obj, mean.function = 'FastExpMean', dispersion.function = 'FastLogVMR')
combined.obj <- ScaleData(combined.obj)
combined.obj <- RunPCA(combined.obj, npcs = 50, features = VariableFeatures(object = combined.obj))
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:50, n.components = 2L) 

saveRDS(combined.obj, "referencemapped.invivo.SEO.downsampled.24211.RDS")
