#*********************************************************************************************************
#* README: This RScript contains code to run the published scRNA-seq/CITE-seq analyses from Parigi et al.,
#* Nature Immunology (2026). To prepare the raw counts matrix used as input, raw fastq files were aligned 
#* to the mouse genome (mm10, v1.20 from 10X Genomics) and quantified using the ‘count’ function in 
#* CellRanger (v6.0.0) with default parameters. Code was originally run using Seurat v.5.0.1. 
#* See 'Methods' for additional detail regarding the uLIPSTIC assay and data processing.
#* 
#* This script calls on "SMP_scRNAseq_object_processing_functions.R", also included in this repository.
#*******************************************************************************************************

library(Seurat)
library(sctransform)
library(data.table)

setwd("~/Documents/Research/FuchsLab/Martina_10xscRNA/scSeq_LipsticTelo_reCal_20220224/")

counts <- as.data.frame(fread("count_matrix.csv")) # Available on GEO (GSE314679).
rownames(counts) <- counts$V1
counts$V1 <- NULL

meta.data <- read.csv("meta_dat.csv") # Available on GEO (GSE314679). nCount_HTO: HTO counts per cell.
rows_meta.data <- rownames(meta.data)
rows_meta.data <- gsub("-",".",rows_meta.data)
rownames(meta.data) <- rows_meta.data

### Create Seurat object from all cells.
setwd("~/Documents/Research/FuchsLab/Martina_10xscRNA/scSeq_LipsticTelo_reCal_20220224/SCT_v1/")

obj <- CreateSeuratObject(counts = counts, meta.data = meta.data)

# Rename clustering metadata from Bioinformatics Resource Center. Will redo initialization of the Seurat
# object from the raw counts, but keep the QC metrics provided in the metadata.
colnames(obj@meta.data)[13] <- "Ji-Dung_RNA_snn_res.0.5"
colnames(obj@meta.data)[14] <- "Ji-Dung_seurat_clusters"

# Prepare Seurat object on all cells.
source("~/Documents/Research/FuchsLab/RScripts/SMP_scRNAseq_object_processing_functions.R")
obj <- seurat_init(obj = obj, cell_type = "all", dims = 30, resolution = 0.1, saveRDS = TRUE)
obj <- readRDS("~/Documents/Research/FuchsLab/Martina_10xscRNA/scSeq_LipsticTelo_reCal_20220224/SCT_v1/all/all_Seurat_object.rds")
obj <- FindClusters(obj, resolution = 0.5)

# Genes used to annotate cell types.
features <- c("Cd3e","Cd4","Cd8a","Foxp3","Klrd1","Itgam","Cd207",
              "Krt10","Krt14","Lrig1","Krt24","Sox9",
              "Areg","Ereg","Egf","Egfr","Lta","Ltb","Ltbr","Lag3","Tgfa",
              "H2ab1","Lgals3")
features <- sort(features)

# Plot features. Used in Fig. 2a,b; Fig. 4a; Fig. 5a; Fig. S3a, S3b; Fig.S6c; Fig. S7C.
print_featureplots(obj = obj, cell_type = "all", features = features,
                   pt.size = 1, alpha = 0.8, order = TRUE, cols = c("lightgrey","darkred"),
                   assay = "SCT")

# Plot UMAP
pdf("./all/all_cells_umap_res0.5.pdf")
p <- DimPlot(obj)
print(p)
dev.off()

### Subset Seurat object.
obj <- readRDS("./all/all_Seurat_object.rds") # Object generated above
corr_hto(obj = obj, cell_type = "all")

# Epithelial clusters
epi_clusters <- c("0","2","7","9")
obj_epi <- subset(x = obj, idents = epi_clusters)
obj <- obj_epi
rm(obj_epi)
obj <- seurat_init(obj, cell_type = "pre_epi", resolution = 0.95, dims = 15)
obj <- readRDS("./pre_epi/pre_epi_Seurat_object.rds")

epi_clusters <- c("0","1","2","3","4","5","7","8","9","10","11","12")
obj_epi <- subset(x = obj, idents = epi_clusters)
obj <- obj_epi
rm(obj_epi)
obj <- seurat_init(obj, cell_type = "epi", resolution = 0.95)

# Immune cells (excluding Langerhans cells).
immune_clusters <- c("1","3","4","6","8")
obj_immune <- subset(x = obj, idents = immune_clusters)
obj <- obj_immune
rm(obj_immune)

obj <- seurat_init(obj, cell_type = "pre_immune", resolution = 0.95)

immune_clusters <- c("0","1","2","3","4","5","6","7","8",
                     "9","10","11","12","13","15","16","17")

obj_immune <- subset(x = obj, idents = immune_clusters)
obj <- obj_immune
rm(obj_immune)

obj <- seurat_init(obj, cell_type = "immune", resolution = 0.95)
print_featureplots(obj = obj, cell_type = "immune", features = features)
print_hto_featureplot(obj = obj, cell_type = "immune") # Used in Fig. 1e.

# Immune cell subsets.
setwd("~/Documents/Research/FuchsLab/Martina_10xscRNA/scSeq_LipsticTelo_reCal_20220224/immune/")

new.cluster.ids <- c("CD4 T", "CD8 T", "CD209a myeloid", "DETC", "Pf4 myeloid", "DETC",
                     "Ear2 myeloid", "Treg","Chil3 myeloid","DETC","Malat1 lymphoid","S100a9 myeloid",
                     "Prolif. T","Rorc T", "ILC", "NK")

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)

pdf("immune_cells_umap_labeled.pdf")
toPlot <- DimPlot(obj, reduction = "umap")
print(toPlot)
dev.off()

saveRDS(obj,file = "immune_Seurat_object.rds")

# Analyses on immune cell subsets.
# CD4 T cells. Used in Fig. S2C.
cd4 <- subset(x = obj, idents = "CD4 T")
cd4 <- seurat_init(cd4, "CD4_T")
corr_hto(obj = cd4, cell_type = "CD4_T")

# CD8 T cells. Used in Fig. S2C.
cd8 <- subset(x = obj, idents = "CD8 T")
cd8<- seurat_init(cd8, "CD8_T")
corr_hto(obj = cd8, cell_type = "CD8_T")

# ILCs. Used in Fig. S2C.
ilc <- subset(x = obj, idents = "ILC")
ilc <- seurat_init(ilc, "ILC")
corr_hto(obj = ilc, cell_type = "ILC")

### DETCs. 
detc <- subset(x = obj, idents = "DETC")
detc <- seurat_init(detc, "DETC")

obj <- readRDS("immune/DETC/DETC_Seurat_object.rds")
obj <- FindClusters(obj, resolution = 0.6, random.seed = 1) # Recluster DETCs.
saveRDS(obj, file = "immune/DETC/DETC_Seurat_object_final.rds")

# Used in Fig. 1f.
p <- DimPlot(obj, pt.size = 2, alpha = 1)
print(p)
dev.off()

# Used in Fig. 1g.
print_hto_featureplot(obj = detc, cell_type = "DETC")

# Used in Fig. 1j; Fig. S3C.
features <- c("Ltb","Xcl1","Areg","Lag3","Lta")
print_featureplots(obj = detc, cell_type = "DETC", features = features)

# Used in Fig. 1i.
corr_hto(obj = detc, cell_type = "DETC")
