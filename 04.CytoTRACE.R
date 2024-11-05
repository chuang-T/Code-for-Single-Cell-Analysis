library(CytoTRACE)
library(tidyverse)
library(Seurat)
library(harmony)
#
dat <- read_rds("all.identied.rds")
dat <- JoinLayers(dat)
dat <- GetAssayData(dat, assay = "RNA") %>% as.matrix()
results <- CytoTRACE(dat)
save(results, file = "EP_cytotrace.Rdata")

FeaturePlot(dat, "CytoTRACE", label = F)

plotdat <- data.frame(cell_type = as.character(dat$cell_type),
                      CytoTRACE = dat$CytoTRACE)
plotdat %>% group_by(cell_type) %>% summarise(aa = median(CytoTRACE)) %>% arrange(-aa) %>% .[, 1] %>% unlist() -> sort_cell_type
plotdat$cell_type <- factor(plotdat$cell_type, levels = sort_cell_type)
ggplot(plotdat, aes(x = cell_type,
                    y = CytoTRACE, fill = cell_type)) +
  geom_boxplot(outlier.alpha = 0) +
  scale_fill_manual(values = my_color) +
  xlab("Cell type")+ 
  theme_test()+NoLegend()->p
ggsave('cytotrace.boxplot.pdf', p, width = 5, height = 3)
