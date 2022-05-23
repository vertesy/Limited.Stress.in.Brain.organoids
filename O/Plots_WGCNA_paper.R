## scWGCNA on Stress Project #2
## Fri Nov 26 12:42:19 2021 ------------------------------


## Oliver Eichmueller
library(gridExtra)
library(igraph);
library(RColorBrewer);
library(tidyverse)
library(WGCNA)
library(gplots)
library(scWGCNA)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Seurat.utils)

sc.obj <- readRDS('WGCNA_test_7.5k_variable/scobj_scores_WGCNA.Rds')


FeaturePlot(sc.obj, features = "magenta_res03_7.5k")
qUMAP("magenta_res03_7.5k", obj = sc.obj, title = "Stress module", save.plot = F) + 
  scale_color_gradient(low = alpha("grey", 0.5), high = "dark magenta") + NoAxes()


geneInfo <- read.csv('WGCNA_test_7.5k_variable/geneInfoSigned_ds_7.5kgenes_res03.csv')


genes <- geneInfo %>% filter(Initially.Assigned.Module.Color %in% "magenta") %>% magrittr::use_series(GeneSymbol)



library(clusterProfiler)
library(enrichplot)

magenta_module <- enrichGO(geneInfo %>% filter(Initially.Assigned.Module.Color=='magenta') %>%
                             dplyr::arrange(desc(kMEmagenta)) %>% 
                             magrittr::use_series(GeneSymbol),
                           OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "BP")

dotplot(magenta_module, showCategory = 30)+
  scale_color_viridis_c(direction = -1)
ggplot2::labs(title = "Stress Module") +
  theme(plot.title = element_text(hjust= 0.5))

ggplot(magenta_module@result %>% 
         mutate(Count_term = stringr::str_split(GeneRatio, "/", simplify = T)[,2]) %>% 
         mutate(pct = round(as.numeric(Count)/as.numeric(Count_term)*100, 2)) %>%
         filter(Count>6) %>%
         arrange(p.adjust) %>% dplyr::slice(1:24) %>%
         mutate(Description = factor(Description, levels = rev(.$Description))), 
       aes(x = pct, y = Description, size = Count, color = -log(p.adjust)))+
  geom_point() + 
  scale_size(breaks = c(3,5,7,10,15,20), limits = c(3,18)) +
  xlim(1,30) + scale_color_gradient() + theme_linedraw()

pdf('WGCNA_test_7.5k_variable/cnet_magenta.pdf', width = 20, height = 20)
cnetplot(magenta_module, showCategory = 40)
dev.off()
magenta_module_pw <- pairwise_termsim(magenta_module)


treeplot(magenta_module_pw, showCategory = 40, label_format = 40, split = F) 
library(openxlsx)
write.xlsx(magenta_module_pw@result[1:40,],file = "WGCNA_test_7.5k_variable/magenta_module_result.xlsx")