The prepare_MAP function takes 3 inputs and outputs 2 files (a pseudobulked and log2 transformed gene expression matrix that MAP expects (This file will be called "MAP_matrix_input.txt"), as well as the .R file to run MAP (this file will be called MAP_script.R)). These two outputs files are meant to be run on the supercomputer using the instructions found in file 2 of this folder! To load this function to your environment, copy, paste, and run the box below into your R studio. Then, scroll below the box to see usage directions.

```
prepare_MAP <- function(spatial_dir, MAP_output_dir, MAP_input_path){
  library(Seurat)
  library(tidyverse)
  library(edgeR)
  
  num_samples <- length(spatial_dir)
  MAP_genes <- c("LY6G6D", "CYP2W1", "TNNC2", "CTTNBP2", "NKD1", "CAB39L", "MLH1", "EPM2AIP1",
                 "SHROOM4", "RNF43", "PRR15", "ATP9A", "H2AFJ", "FARP1", "TCF7", "MAPRE3",
                 "ZMYND8", "DDX27", "TGFBR2", "PIWIL4", "FECH", "DOCK5", "TYMS", "HPSE",
                 "ASPHD2", "AGR2", "GFI1", "RPL22L1", "RAB27B", "GNLY", "DUSP4")
  
  #Loading first sample
  data.dir <- spatial_dir[1]
  seurat <- Load10X_Spatial(data.dir, filename = "filtered_feature_bc_matrix.h5")
  
  #Filtering
  seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
  seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")
  seurat <- subset(
    seurat, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 &
      nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)
  
  #Pseudobulk
  seurat_sum <- AggregateExpression(
    seurat,
    group.by = "orig.ident",
    normalization.method = NULL,
    verbose = TRUE
  )
  seurat_sum <- as.matrix(seurat_sum$Spatial)
  colnames(seurat_sum) <- basename(data.dir) #assigns sample name
  
  #log normalization via EdgeR for input to MAP
  seurat_edgeR <- DGEList(counts = seurat_sum)
  keep <- filterByExpr(seurat_edgeR, min.count = 50) #default 10, copy paper that puts it to 50
  seurat_edgeR$counts[!keep, ] <- 0
  seurat_log2_cpm <- cpm(seurat_edgeR, log = TRUE)
  
  #Checking that all genes are there and if not, add them in with a value of 0
  for(j in 1:length(MAP_genes)){
    if(MAP_genes[j] %in% rownames(seurat_log2_cpm)){}
    else{
      #Working but apparently inefficient approach
      #seurat_log2_cpm <- rbind(seurat_log2_cpm, 0)
      #rownames(seurat_log2_cpm)[nrow(seurat_log2_cpm)] <- MAP_genes[i]
      
      # 1. Create the new row as a matrix with the correct name and column count
      new_row <- matrix(0, nrow = 1, ncol = ncol(seurat_log2_cpm), 
                        dimnames = list(MAP_genes[j], colnames(seurat_log2_cpm)))
      
      # 2. Bind them together safely
      seurat_log2_cpm <- rbind(seurat_log2_cpm, new_row)
    }
  }
  final_matrix <- seurat_log2_cpm
  
  if(num_samples>1){
    for(i in 2:num_samples){
      #Loading samples
      data.dir <- spatial_dir[i]
      seurat <- Load10X_Spatial(data.dir, filename = "filtered_feature_bc_matrix.h5")
      
      #Filtering
      seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
      seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")
      seurat <- subset(
        seurat, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 &
          nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)
      
      #Pseudobulk
      seurat_sum <- AggregateExpression(
        seurat,
        group.by = "orig.ident",
        normalization.method = NULL,
        verbose = TRUE
      )
      seurat_sum <- as.matrix(seurat_sum$Spatial)
      colnames(seurat_sum) <- basename(data.dir)
      
      #log normalization via EdgeR for input to MAP
      seurat_edgeR <- DGEList(counts = seurat_sum)
      keep <- filterByExpr(seurat_edgeR, min.count = 50) #default 10, copy paper that puts it to 50
      seurat_edgeR$counts[!keep, ] <- 0
      seurat_log2_cpm <- cpm(seurat_edgeR, log = TRUE)
      
      #Checking that all genes are there and if not, add them in with a value of 0
      for(j in 1:length(MAP_genes)){
        if(MAP_genes[j] %in% rownames(seurat_log2_cpm)){}
        else{
          #Working but apparently inefficient approach
          #seurat_log2_cpm <- rbind(seurat_log2_cpm, 0)
          #rownames(seurat_log2_cpm)[nrow(seurat_log2_cpm)] <- MAP_genes[i]
          
          # 1. Create the new row as a matrix with the correct name and column count
          new_row <- matrix(0, nrow = 1, ncol = ncol(seurat_log2_cpm), 
                            dimnames = list(MAP_genes[j], colnames(seurat_log2_cpm)))
          
          # 2. Bind them together safely
          seurat_log2_cpm <- rbind(seurat_log2_cpm, new_row)
        }
      }
      
      # 1. Put all your single-column matrices into a list
      matrix_list <- list(final_matrix, seurat_log2_cpm)
      
      # 2. Merge them all together by row name
      merged_df <- Reduce(function(x, y) merge(x, y, by = 0, all = TRUE), matrix_list)
      
      # 3. Restore the row names and clean up the artifact column
      rownames(merged_df) <- merged_df$Row.names
      merged_df$Row.names <- NULL
      
      # 4. Fill all NA values with 0
      merged_df[is.na(merged_df)] <- 0
      
      # 5. Convert back into a matrix
      final_matrix <- as.matrix(merged_df)
    }
  }
  write.table(final_matrix, 
              "MAP_matrix_input.txt",
              sep = "\t", 
              quote = FALSE, 
              row.names = TRUE, 
              col.names = NA)
  
  my_code <- paste0("library(MAP)", "\ntestP = '", MAP_output_dir, "'", "\nexpData = '", MAP_input_path, "'", "\nrunMAP(expData, testP)")
  writeLines(my_code, "MAP_script.R")
}
```

The first input is is spatial_dir, which can either be a single pathname to a spatial transcriptomics input files folder, or a LIST of pathnames pointing to multiple ST input files folders. This function is meant to accommodate multiple samples at once so that MAP can be run more efficiently. This input will tell the function what samples to pseudobulk and add to the gene expression matrix. This file will be called "MAP_matrix_input.txt"

The second input is MAP_output_dir, which is the directory where you want the MAP results to be deposited. It doesn't have to be an existing folder, as it will be created by the script when it is run later on, but you want this to be a directory on the supercomputer as this is where you will be running the script.

The third input is MAP_input_path, which is the supercomputer directory where you will store the gene expression matrix that you get from this command. You can rename the pseudobulked gene expression file later but recall that it is named MAP_matrix_input.txt by default!

An example of how to run the function:
```
#use c() to create a list with all the sample directories listed
samples <- c("/Users/Me/Downloads/Sample1", "/Users/Me/Downloads/Sample2", "/Users/Me/Downloads/Sample3", "/Users/Me/Downloads/Sample4")

#run the function in just one line! The outputs will appear in your working directory on R studio
prepare_MAP(samples, "/expanse/lustre/projects/csd670/akao1/MAP_output", "/expanse/lustre/projects/csd670/akao1/MAP_matrix_input.txt")
```
