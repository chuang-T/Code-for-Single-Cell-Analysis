.libPaths("E:/R_lib/4.2.2/")
library(Seurat)
library(optparse)
library(tidyverse)
library(openxlsx)
library(harmony)
library(RColorBrewer)
#
files <- list.files(my_dir, # The address of matrix  file
                    full.names = T,
                    recursive = T) 

my_read <- function(x){
  aa <- Read10X_h5(x)
  newname <- str_split(x, '/', simplify = T)[5]
  aa <- CreateSeuratObject(counts = aa, project = newname)
  aa <- RenameCells(aa, new.names = paste0(aa$orig.ident, "_",Cells(aa)))
  return(aa)
}

all <- lapply(files, my_read)
all <- merge(all[[1]], c(all[[2]], all[[3]], all[[4]], all[[5]], all[[6]]))

# filter
MT_genes <- c("ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6")
all[["percent.mt"]] <- PercentageFeatureSet(all, assay = "RNA", feature = MT_genes)

cut_counts <- c(1000, 30000)
cut_feature <- 300
cut_mt <- 10
min.cells <- 10
{
  all <- subset(all,subset = c(nCount_RNA >= cut_counts[1] & nCount_RNA <= cut_counts[2]))
  all <- subset(all,subset = nFeature_RNA >= cut_feature)
  all <- subset(all, subset = percent.mt <= cut_mt) 
  num.cells <- rowSums(as.matrix(all@assays$RNA@counts) > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  all <- all[rownames(all) %in% genes.use, ]  
  
}

new_name <- str_replace_all(all$orig.ident, '-', '_') %>% 
  str_split("_", simplify = T)

all$sample <- paste0(new_name[, 2], '_',new_name[, 3])
all$sample <- factor(all$sample, levels = c('SPF_JPP', 'GF_JPP', 'SPF_MLN', 'GF_MLN','SPF_SPLE','GF_SPLE'))

all$tissue <- new_name[, 3]
all$tissue <- factor(all$tissue, levels = c('JPP', 'MLN', 'SPLE'))

all$group <- new_name[, 2]
all$group <- factor(all$group, levels = c('SPF', 'GF'))

all <- all %>% 
  NormalizeData(verbose = T) %>%
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
DimPlot(all)

R <- c(0.2) 
all <- FindClusters(all, resolution = R)

# cell type identy
# Tcell
FeaturePlot(all, features = "CD3E", label = T, cols = c('grey','red'))+
  # B cell
  FeaturePlot(all, features = "CD79B", label = T, cols = c('grey','red'))
# Myeloid 
FeaturePlot(all, features = c('CST3','CD68', 'CD163', 'CD14'), label = T, cols = c('grey','red'))

saveRDS(all, file = "all.rds")
