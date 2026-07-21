# Pseudobulking
Spatial transcriptomics (ST) is a technique that maps gene expression at (near) single cell resolution onto an intact tissue sample. This resolution allows us to infer and quantify cell level pathways and communication. It allows for the analysis of cell type and pathway heterogeneity within a single tissue slice, a significant improvement over the bulk RNA sequencing ([bulk RNA-seq](https://www.cd-genomics.com/resource-bulk-rna.html)) technique which measures the average gene expression level across the entire tissue. 

As a newer technology, one limitation of spatial transcriptomics is that many bioinformatics tools and software have been developed for bulk RNA sequencing data, operating at the averaged tissue sample level instead of (near) single cell level. The solution to this problem? [Pseudobulking!](https://www.elucidata.io/blog/navigating-the-single-cell-ship-through-the-pseudo-bulk-route)

Essentially, gene expression data from multiple spots/cells are aggregated (summed) to create a "representative gene expression profile". This can be summed from cells of the same cell type (e.g. summing gene expression counts from only neurons, or only epithelial cells), or from the entire sample. 

The most common reason pseudobulking is done is is for differential expression (DE) analysis of single cell RNA sequencing data. Many DE tools (such as DESeq2 or EdgeR) are specifically developed for Bulk RNA sequencing. As such, pseudo-bulking is used so that we may apply these bulk RNA-seq tools to single cell or spatial transcriptomics data. In addition, the ability to pseudo-bulk within just a specific annotated cell type means that we may apply bulk RNA-seq techniques while preserving some aspect of our finer resolution.

Other benefits of pseudobulking to note:
* Treating each cell as an individual sample artificially inflates p-values. Cells from the same sample/patient are inherently similar/dependent; pseudobulking allows us to avoid this error (known as "pseudoreplication")
* Pseudobulking increases the effective [sequencing depth](https://3billion.io/blog/sequencing-depth-vs-coverage) of the data (number of times a specific base (nucleotide) in the DNA is read), which means the data is more robust. Especially considering the phenomenon of ["dropouts"](https://www.nature.com/articles/s41467-020-14976-9) (a cell has no gene expression detected for gene X even though that gene is known to be expressed in that cell type) due to technical and biological limitations, pseudobulking helps to get a more accurate representative profile of gene expression.

For this project, we want to apply the Microsatellite instability Absolute single sample Predictor (MAP) software to classify or spatial transcriptomics samples as Microsatellite Instable (MSI) or Microsatellite Stable (MSS). This tool was trained on and only works on Bulk RNA sequencing data, meaning that we need to pseudobulk all the spots across the entire sample. This workflow (pseudobulking spatial transcriptomics samples -> applying MAP) was done by [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10781769/) and serves as our precedent.

## Table of Contents

1: Loading Data/Data Entry
* 1.1: filtered_feature_bc_matrix.h5 file (most common!)
* 1.2: .mtx file
* 1.3: .rds file

2: Data Quality Control and Preprocessing

3: Pseudobulk step: Aggregating data (AggregateExpression function)

4: Log-normalization via EdgeR
* 4.1: Create DGEList object
* 4.2: Filter low counts (filterByExpr function)
* 4.3: Generate log2(CPM) matrix

5: Prepare results as an input .txt file for MAP

Note: Steps 1 and 2 are directly copied from my original [spatial transcritpomics notes](https://github.com/alfalfacow/Lab-Notebook/blob/main/Spatial%20Transcriptomics/1-Introductory-Workflow.md)!!

## Step 1: Loading Data/Data Entry
First, you will have to install and initialize the Seurat package on an R studio window, as well as the other necessary packages for the workflow.
```
install.packages("Seurat") #if Seurat not already installed
library(Seurat)
library(tidyverse)
BiocManager::install("edgeR") #if edgeR not already installed)
library(edgeR)
```

Next, you will need to obtain a spatial transcriptomics dataset, either from a publically accessible dataset or one you generated yourself. An example is the [GSE281978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281978) series from the Gene Expression Omnibus (GEO) for Head and Neck Cancer (HNSC). The dataset can be downloaded as a tar file at the bottom of the page.

NOTE: This page will mostly focus on 10x Genomics Visium spatial transcriptomics datasets.

There are many different data formats for spatial transcriptomics that may be present on GEO or other datasets. We will focus on three different ways that gene expression data is stored

### 1.1: filtered_feature_bc_matrix.h5 file
The most common (in my experience) main output file from a Visium dataset is the filtered_feature_bc_matrix.h5, which contains the gene expression data. However, there are other files that make up the dataset that fully distinguish it as a spatial transcriptomics dataset, including the image file of the tissue sample (a .png file) the coordinates of each Visium "cell spot" (tissue_positions_list.csv), and a "scale factor" file that helps to scale the cell spots and map the coordinates onto the image (scalefactors_json.json).

To read the dataset into Seurat as a Seurat object, the files for EACH specific sample need to be organized in a very specific structure. This specific structure is recognized with Seurat's Load10X_Spatial command.

* (Folder) Sample name
* -> (File) (SampleName)_filtered_feature_bc_matrix.h5
* -> (Folder) (SampleName)
* ->-> (File) tissue_lowres_image.png
* ->-> (File) scalefactors_json.json
* ->-> (File) tissue_positions_list.csv

You may need to manually change the contents and names of the files in order to achieve this structure. Now we are ready to read our data on R studio!

```
#Load in Necessary Packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
set.seed(123) #Some steps (like UMAP) depend on random number generators, so we should set a seed in order to make our work reproducible! 123 is an arbitrary number and can be anything

data.dir <- "Path/To/The/SampleName/Folder/Containing/Your/Dataset" #ex: /Users/alfred/Desktop/GSM8633891
NameOfSeuratObject <- Load10X_Spatial(data.dir, filename="(SampleName)_filtered_feature_bc_matrix.h5") #filename should be whatever name is in the folder for the feature matrix file
glimpse(HNSCC) #Allows you to look at the type of data included in a Seurat Object

```
Sometimes, the files may come in a .gz format (compressed file format known as gzip), similar to a .zip file. To read these files, they must be unzipped first. This can be done through the terminal/command line using the gunzip function:
```
cd /path/to/folder/with/files/you/want/to/unzip #cd = change directory, into the folder with the zipped files for easy access
gunzip file_name.format.gz #this gunzips a single file
gunzip file1.format.gz file2.format.gz file3.format.gz #this gunzips multiple files
gunzip *.gz #this gunzips all files in the current directory (from the cd command) that end in .gz
```
If you are curious, for other must-know linux commands for use on command terminal, see [section 0](https://github.com/alfalfacow/Lab-Notebook/blob/main/Supercomputer/0-Must-Know-Linux-Commands) of the Supercomputer "chapter" of this lab notebook!

### 1.2: .mtx file
Another format for the gene expression matrix is the .mtx format (usually named (sample name)_matrix.mtx; [see this example](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9322957). This requires a slightly different workflow for loading in the data into Seurat: since the Load10X_Spatial command recognizes .h5 files, we have to convert the .mtx file into an .h5 file.

To do this, we will need 3 components from the dataset to make up the .h5 file: "matrix.mtx" (gene expression matrix), "features.tsv" (gene names), and "barcodes.tsv" (name of each spot, labeled with a genetic barcode). We will use the [write10x](https://rdrr.io/bioc/DropletUtils/man/write10xCounts.html) function from the DropletUtils package.

Begin by installing and loading the DropletUtils package from Biocondutor!
```
BiocManager::install("DropletUtils") #install package from Bioconductor
library(DropletUtils) #load package into current R session
```
Then, we need to prepare the necessary files. Move the aforementioned 3 components (matrix.mtx, features.tsv, and barcodes.tsv) into the same folder on your computer. You will need to note the path to this folder (example: this format /home/user/desktop, or something similar. You can find this on a mac by right clicking the folder in your file explorer and holding the option button, then clicking "copy (filename) as pathname".)

Format:
* (Folder) arbitrary_folder_name
* -> (File) matrix.mtx
* -> (File) features.tsv
* -> (File) barcodes.tsv

Next, we can load these files as a "sparse matrix" object without any of the spatial data involved YET, using the Read10X command from Seurat. This command expects the files to be gzipped (.gz), so you do NOT have to run the gunzip command like in section 1.1.
```
filter_matrix <- Read10X("path/to/arbitrary_folder_name")
```

Then, we use the write10xCounts command to convert this matrix object into a .h5 file format:
```
write10xCounts(
  path = "/path/to/desired/location/filtered_feature_bc_matrix.h5", #this line specifies the path to the folder/location where the .h5 file                                                                         will be made, as well as the name of the .h5 file
  x = filter_matrix, #this is the matrix object we created earlier
  genome = "GRCh38", #current human reference genome used for gene alignment
  type = "HDF5", #this is the type of file we want to make (HDF5 is same as .h5 file format)
  version = "3", #current version of CellRanger
  gene.id = rownames(filter_matrix),
  gene.symbol = rownames(filter_matrix) 
)
```
Now you can use this file to load the data into Seurat using the steps from section 1.1!

You might be wondering: why are there different formats if there are technically able to be interconverted? According to the [CellRanger website](https://www.10xgenomics.com/support/jp/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices) (a tool from 10X genomics that converts the raw expression data from single cell and spatial transcriptomics pipelines into gene expression matrices), earlier versions of CellRanger (v3 and before) output .mtx, while later versions output the filtered_feature_bc_matrix.h5 file).

### 1.3 .rds file
For GEO datasets this pretty much never appears, but for other public repositories, researchers may choose to save their own Seurat objects as a .rds file. A .rds file stands for "R Data Serialization", with [serialization](https://www.geeksforgeeks.org/r-language/data-serialization-rds-using-r/) defined as "the process of converting complex data structures into an understandable format, suitable for storage and transmission".

To load a .rds file into the environment, simply use the [readRDS()](https://rdrr.io/r/base/readRDS.html) command!
```
readRDS("file_name.rds") #use this if is already within the R studio working directory
readRDS("/path/to/file_name.rds") #use this if the rds file is not within the R studio working directory
```
There is also a .rda format, which can store multiple objects. Load this using the [load() function](https://www.r-bloggers.com/2017/04/load-save-and-rda-files/) (very simple!)
```
load("file.rda")
load("/path/to/file.rda")
```
Again, this format is pretty rare for spatial transcriptomics/visium datasets, but now you will be prepared for it! 

## Step 2: Data Quality Control
To make sure that our data is high quality, the next steps are to quality control and normalize the data. We want to filter out low quality spots (quality control), including spots with with too little or too many genes detected, or spots with a high proportion of unwanted mitochondrial and ribosomal DNA. We also want to normalize the data to account for "variance in molecular counts" among different spots of the tissue (due to technology imperfections and biological differences). 

The four main metrics commonly controlled are:
1. nFeature_Spatial (number of unique genes detected at each spot, part of Seurat object)
2. nCount_Spatial (number of DNA molecules detected at each spot, part of Seurat object)
3. percent.mt (percent of DNA that is from mitochondria, must be manually calculated)
4. percent.ribo (percent of DNA that is from ribosomes, must be manually calculated

The first two have been stored as part of the meta.data of the Seurat object and do not need to be calculated. However, Seurat does not calculate the last two metrics, so we must do that ourselves:

```
#SeuratObject[[meta.data.column]] accesses or creates a new column in the meta.data section of the Seurat object
#PercentageFeatureSet calculates the percent of all counts that include a certain pattern
#The following code calculates the percent of all counts that belong to mitochondrial DNA (named "MT-XX" or "RPS"/"RPL", apparently Regex notation) and puts it into newly defined meta.data columns

SeuratObject[["percent.mt"]] <- PercentageFeatureSet(object = SeuratObject, pattern = "^MT-")
SeuratObject[["percent.ribo"]] <- PercentageFeatureSet(SeuratObject, pattern = "^RP[SL]")

#Visualize the current distribution of each of these metrics with a violin plot (optional)
VlnPlot(
  SeuratObject, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Visualize each metric on the tissue image using the SpatialFeaturePlot() function (optional)
SpatialFeaturePlot(
  SeuratObject, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent-ribo")) &
  theme(legend.position = "bottom")

```
The violin plots should show that each of these metrics are rather widely distributed, with some points close to 0 and others ranging from 10 thousand for nFeature_Spatial to 40 thousand for nCount_Spatial! Our next step is to filter out low quality points of each of these metrics.

Common thresholds for quality control have been taken from this [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11655296/):
1. nFeature_Spatial less than 200 or greater than 7,500
2. nCount_Spatial less than 250 or greater than 50,000
3. percent.mt > 15%
4. percent.ribo > 40%

```
#Updates Seurat object as a subset of the original using the subset() function, keeping the features that are NOT part of the exclusion criteria defined above

SeuratObject <- subset(
  SeuratObject, subset = nFeature_Spatial < 7500 & nFeature_Spatial > 200 &
  nCount_Spatial < 50000 & nCount_Spatial > 250 & percent.mt < 15 & percent.ribo < 40)

#You may plot new violin plots to verify filtering (optional)
VlnPlot(
  SeuratObject, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt", "percent.ribo"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

```
## Once this has been done, you are ready to proceed to the next step!! Do NOT normalize the data (EdgeR requires raw counts without normalization, and we will apply normalization methods using EdgeR directly).

## Step 3: Aggregating data for Pseudobulk
We will be using Seurat's AggregateExpression() function to sum gene expression counts from our sample!
```
SeuratObject_sum <- AggregateExpression(
  SeuratObject,
  group.by = "orig.ident",
  normalization.method = NULL,
  verbose = TRUE
)
```
We will group by orig.ident, which is the metadata section that should defines the tissue sample (thus it is the same for all spots from the same sample). If we wanted to group by cell type, we would use a different metadata section in group.by, but we do not need to worry about that right now. Also, we are NOT normalizing the data, so we set normalization.method as NULL!

Then, we want to update the output as a matrix so it can be read by edgeR in the next step (it is currently a "list" object:
```
SeuratObject_sum <- as.matrix(SeuratObject_sum$Spatial) #The output matrix data we want is within "Spatial" of the "SeuratObject_sum", so we invoke that using the $ operator.
```
At this point, your SeuratObject_sum should be a matrix with rows as gene names and columns as samples (should only be one sample). You should double check the number of genes in this matrix and make sure it matches the original number in your filtered Seurat object!
```
dim(SeuratObject) #this returns the dimensions of the original expression matrix: number of genes and then number of cells
dim(SeuratObject_sum) #this returns the dimensions of the summed matrix: number of genes (should be same as above) and 1 (one sample)
```

## Step 4: Log-normalization via EdgeR
The MAP software asks for the input to be in "log2 transformed expression". However, our aggregated only has raw summed counts. To transform our raw counts, we will use the EdgeR workflow (which is typically used to transform/normalize raw Bulk RNA sequencing data).

### 4.1: Create DGEList object
We will use the matrix we created in step 3 as an input to create a DGEList object (standard input for EdgeR pipeines):
```
edgeR <- DGEList(counts = SeuratObject_sum)
```

### 4.2: Filter low counts (filterByExpr function)
Next, we will filter out genes with lower than a certain threshold of counts. Typically the default is 10 counts minimum, but the [precedent paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10781769/) sets this to 50. Also, typically these genes are removed directly. However for our purposes, instead of filtering the genes out of the data entirely, we want to just set any gene with counts under 50 to 0, because MAP expects expression values of all genes.
```
keep <- filterByExpr(edgeR, min.count = 50) #this function returns a list of all genes ABOVE the threshold set by "min.count"
edgeR$counts[!keep, ] <- 0 #this will set all rows (gene names) that are NOT part of the "keep" list (thus any under 50) to 0
```

### 4.3: Generate log2(CPM) matrix
The penultimate step of our workflow is to go from raw counts to a log-normalized matrix (As mentioned before, MAP expects an input of log2 normalized expression data in a matrix)! Typically this step would be preceeded by something called "Trimmed Mean of M-values" (TMM) normalization (calcNormFactors). However, this is only done when there are multiple samples of interest being compared for differential expression analysis, which we are not interested in. We are only interested in single samples, which means we can skip this step.

We will use edgeR's cpm() function to normalize by Counts Per Million (CPM) (read more about gener expression units [here](https://www.reneshbedre.com/blog/expression_units.html#google_vignette)). CPM is calculated by the following formula, which edgeR will automatically apply for us: gene expression counts * 1,000,000 / total counts for the entire sample. 

```
edgeR_cpm <- cpm(edgeR, log = TRUE) #setting log = TRUE gives us log2(CPM+2). The +2 is a default offset in edgeR to prevent log(0), which will throw an error
```
This should output a matrix object with the log2(CPM+2) normalized data!!!

## Step 5: Prepare results as an input .txt file for MAP
After doing all the hard work to calculate and output our pseudobulked, log2transformed and CPM normalized expression matrix, we are now ready to export it as a .txt file for input into the MAP tool!

First we will verify that the 31 marker genes that MAP relies on are included in the genes list
```
#reference genes list
MAP_genes <- c("LY6G6D", "CYP2W1", "TNNC2", "CTTNBP2", "NKD1", "CAB39L", "MLH1", "EPM2AIP1",
  "SHROOM4", "RNF43", "PRR15", "ATP9A", "H2AFJ", "FARP1", "TCF7", "MAPRE3",
  "ZMYND8", "DDX27", "TGFBR2", "PIWIL4", "FECH", "DOCK5", "TYMS", "HPSE",
  "ASPHD2", "AGR2", "GFI1", "RPL22L1", "RAB27B", "GNLY", "DUSP4")

#verify 31 genes of interest are in the matrix
MAP_genes %in% rownames(edgeR_cpm)
```
Ideally you should see 31 "TRUE" statements printed in the output section of R studio.

Finally, we are ready to export the .txt file!
```
write.table(edgeR_cpm, 
            "MAP_input.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA)
```

Congratulations! You have successfully generated the input file for the MAP tool. The .txt file will be in the same folder that the R studio project you are currently on is located :)

To see how to run this input into MAP, see section 2.
