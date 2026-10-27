#*********************************************************************************************************
#* README: This RScript defines functions used by "SMP_scRNAseq_object_processing_functions.R", 
#* also included in this repository.
#*******************************************************************************************************

# Initialize Seurat objects. Called as seurat_init(obj = obj, cell_type = "name"). Prepares a Seurat object
# normalized with SCTransform (vst.flavor = v1), and runs PCA, UMAP, and clustering, along with basic
# plot generation, identification of cluster-specific markers, and generation of an rds object. Files are
# output into a directory named in the "cell_type" parameter.
seurat_init <- function(obj, cell_type = deparse(substitute(obj)), dims = 30, resolution = 0.8, saveRDS = TRUE){
  
  # Prepare output directory
  working_dir <- getwd()
  new_dir <- paste(working_dir,"/",cell_type, sep = "")
  if(!dir.exists(new_dir)){
    dir.create(new_dir)
  }
  setwd(new_dir)
  
  # Remove SCT transformation if already present
  DefaultAssay(obj) <- "RNA"
  if("SCT" %in% names(obj@assays)){
    obj[["SCT"]] <- NULL
  }
  
  # Normalize uLIPSTIC signal (HTO counts) by total read depth per cell
  if(!("nCount_HTO_scaled" %in% colnames(obj@meta.data))){
    hto <- obj$nCount_HTO
    scale_factor <- obj$nCount_RNA + obj$nCount_HTO
    scale_factor <- t(scale_factor)
    hto_scaled <- t((hto/scale_factor)*10^3)
    obj$nCount_HTO_scaled <- hto_scaled
  }
  
  # Normalize, run PCA
  obj <- SCTransform(obj, vst.flavor = "v1") # Default seed 1448145, n variable features = 3000
  obj <- RunPCA(obj) # Default seed 42
  
  # Elbow plots for QC of dimensions
  filename <- paste(cell_type,"PCA_ElbowPlot.pdf")
  pdf(filename)
  toPlot <- ElbowPlot(obj)
  print(toPlot)
  dev.off()
  
  # Run UMAP and clustering; print plots
  obj <- RunUMAP(obj, dims = 1:dims)
  set.seed(147)
  obj <- FindNeighbors(obj, dims = 1:dims)
  set.seed(147)
  obj <- FindClusters(obj, resolution = resolution) # Default resolution = 0.8
  
  # Plot UMAP
  filename <- paste(cell_type,"_cells_umap.pdf",sep = "")
  pdf(filename)
  toPlot <- DimPlot(obj, reduction = "umap")
  print(toPlot)
  dev.off()
  
  # Get cluster markers
  markers <- FindAllMarkers(obj)
  filename <- paste(cell_type,"_cells_cluster_markers.txt",sep = "")
  write.table(markers, filename, sep = "\t")
  
  # Output rds object
  Idents(obj) <- "seurat_clusters"
  if(saveRDS == TRUE){
    filename <- paste(cell_type,"_Seurat_object.rds",sep = "")
    saveRDS(obj, file = filename)
  }
  
  setwd(working_dir)
  return(obj)
}

# Function used to correlate the RNA expression of each gene within a Seurat 
# object to uLIPSTIC labeling signal defined by CITE-seq (HTO antibodies). 
# Output directory is defined by "cell_type". Outputs Spearman correlations
# with summary statistics.
corr_hto <- function(obj, cell_type = deparse(substitute(obj))){
  
  # Set output directory
  working_dir <- getwd()
  new_dir <- paste(working_dir,"/",cell_type, sep = "")
  if(!dir.exists(new_dir)){
    dir.create(new_dir)
  }
  setwd(new_dir)
  
  # Normalize uLIPSTIC signal (HTO counts) by total read depth
  hto <- obj$nCount_HTO
  scale_factor <- obj$nCount_RNA + obj$nCount_HTO
  
  library(matrixStats)
  norm_counts <- as.data.frame(obj@assays$SCT@data)
  norm_counts <- norm_counts[(rowSds(as.matrix(norm_counts)) != 0),]
  
  hto <- t(hto)
  
  corr_input <- hto
  rownames <- "hto_raw"
  rownames(corr_input) <- rownames
  
  corr_input <- rbind(corr_input, norm_counts)
  corr_input <- t(corr_input)
  mean_exp <- colMeans(corr_input[,2:ncol(corr_input)])
  
  scale_factor <- t(scale_factor)
  hto_scaled <- (hto/scale_factor)*10^3
  
  corr_input_scaled <- hto_scaled
  rownames <- "hto_scaled"
  rownames(corr_input_scaled) <- rownames
  
  corr_input_scaled <- rbind(corr_input_scaled, norm_counts)
  corr_input_scaled <- t(corr_input_scaled)
  
  # Get Spearman correlations between normalized counts for each gene and 
  # normalized uLIPSTIC signal, along with associated p-values.
  corr_matrix_norm_spearman <- t(cor(corr_input_scaled[,1],corr_input_scaled[,-1],method = "spearman"))
  pvals <- apply(corr_input_scaled[,-1], 2, function(x) {
    cor.test(corr_input_scaled[,1], x, method = "spearman")$p.value
  })
  padj <- p.adjust(pvals, method = "fdr")
  
  # Append summary statistics to results.
  corr_matrix_norm_spearman <- cbind(corr_matrix_norm_spearman,pvals)
  corr_matrix_norm_spearman <- cbind(corr_matrix_norm_spearman,padj)
  corr_matrix_norm_spearman <- cbind(corr_matrix_norm_spearman,mean_exp)
  colnames(corr_matrix_norm_spearman) <- c("r","pvalue","fdr","mean_expression")
  corr_matrix_norm_spearman <- as.data.frame(corr_matrix_norm_spearman)
  corr_matrix_norm_spearman <- corr_matrix_norm_spearman[order(corr_matrix_norm_spearman$fdr),]
  
  # Output tables
  filename <- paste(cell_type, "_hto_correlation_spearman_read-normalized.txt", sep = "")
  write.table(corr_matrix_norm_spearman, filename, sep = "\t")

  setwd(working_dir)
}

# Function used to plot scRNA-seq gene expression over UMAP.
print_featureplots <- function(obj, cell_type = deparse(substitute(obj)), features, 
                               palette = "inferno", pt.size = 2, alpha = 1, 
                               order = TRUE, cols = NULL){
  
  working_dir <- getwd()
  new_dir <- paste(working_dir,"/",cell_type, sep = "")
  if(!dir.exists(new_dir)){
    dir.create(new_dir)
  }
  setwd(new_dir)
  
  filename = paste(cell_type, "_cells_featureplots.pdf", sep = "")
  pdf(filename, onefile = TRUE)
  present_features <- c(rownames(GetAssayData(object = obj)), "nCount_HTO")
  for(f in features){
    if(f %in% present_features){
      if(!is.null(cols)){
        toPlot <- FeaturePlot(obj, reduction = "umap", features = f, pt.size = pt.size, alpha = alpha, cols = cols)
      } else{
        toPlot <- FeaturePlot(obj, reduction = "umap", features = f, pt.size = pt.size, alpha = alpha, order = order) & 
          scale_color_viridis_c(alpha = alpha, option = palette)  
      }
      print(toPlot)
    }
  }
  dev.off()
  
  setwd(working_dir)
}

# Function used to plot uLIPSTIC signal (HTO) over UMAP.
print_hto_featureplot <- function(obj, cell_type = deparse(substitute(obj)), cols = c("lightgrey","darkgreen"),
                                  max.cutoff = 300, order = TRUE, pt.size = 2, alpha = 1){
  
  working_dir <- getwd()
  new_dir <- paste(working_dir,"/",cell_type, sep = "")
  if(!dir.exists(new_dir)){
    dir.create(new_dir)
  }
  setwd(new_dir)
  
  filename = paste(cell_type, "_cells_hto_featureplots.pdf", sep = "")
  pdf(filename)
  f <- "nCount_HTO"
  toPlot <- FeaturePlot(obj, reduction = "umap", features = f, cols = cols, max.cutoff = max.cutoff, 
                        order = order, pt.size = pt.size, alpha = alpha)
  print(toPlot)
  dev.off()
  
  setwd(working_dir)
}
