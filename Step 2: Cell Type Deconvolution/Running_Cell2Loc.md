# Cell2Location
As you may be aware, [Cell2Location](https://cell2location.readthedocs.io/en/latest/index.html) is a tool for "deconvolution" of spatial transcriptomics data. While spatial transcriptomics indeed approaches near single cell resolution (with more and more advanced tools such as 10X Genomics' Atera reaching subcellular resolution), many technologices such as Visium may have more than one cell per spot profiled. Especially considering that Visium spots are 55 micrometers in diameter (while most cells are far below that), there are estimates of up to 10 cells per spot. Thus, labeling each spot as a single cell is not accurate.

Cell2Location bridges this obstacle by estimating how many of each cell type is present within a spot (cell type abundance per spot). This is done by first training Cell2Location using a reference single cell RNA sequencing dataset with cell types already annotated. Cell2Location learns the transcriptomic profile from this reference dataset and applies it to the spatial transcriptomics dataset.

For the purposes of this project, I have already trained Cell2Location on the reference datasets ([GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735) and [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)), both of which have cell types prelabeled. These were chosen due to their use in previous CRC studies using cell2location, and also because it includes cells labeled by consensus molecular subtype (CMS), another CRC-specific feature we may be interested in.

You may have already run Cell2Location for the pancreatic cancer project, so this is mostly instructions for adapting that script for use on the supercomputer. That way, you won't have to leave it running with your computer on all day. Running this on the supercomputer is optional! If you prefer, feel free to continue to run it locally (it may actually be slightly faster since it can utilize GPU, and you may already know how to do it).

These instructions assume familiarity with Expanse/HPC computing, which you have hopefully gained from trying out the [MSI annotation instructions](https://github.com/alfalfacow/Colorectal-Cancer-Spatial-Transcriptomics/tree/main/Step%201%3A%20MSI%20Annotation).

## Table of Contents

1: Preparing the files

2: Running Cell2Location
* 2.1: prepared script on your end
* 2.2: no prepared script

## Step 1: Preparing the files
As input, Cell2Location requires a single cell reference dataset and a spatial transcriptomics dataset, both in the .h5ad file format. I have prepared the single cell reference data already (so don't worry about it), so below is a guide on how to get a .h5ad file from a spatial transcriptomics dataset (Visium)!

We will be using squidpy's [squidpy.read.visium](https://squidpy.readthedocs.io/en/stable/api/squidpy.read.visium.html) function into an AnnData object. Like Seurat's [Load10XSpatial()](https://satijalab.org/seurat/reference/load10x_spatial) function, this function looks for a very specific file structure:

* (Folder) Sample name
* -> (File) filtered_feature_bc_matrix.h5
* -> (Folder) Spatial
* ->-> (File) tissue_lowres_image.png
* ->-> (File) scalefactors_json.json
* ->-> (File) tissue_positions_list.csv

Once you have the command ready, run the following python script:
```
import squidpy as sq

#Reads visium files into AnnData object
adata = sq.read.visium(
    path="/path/to/your/folder/with/Sample_Name",
    counts_file="filtered_feature_bc_matrix.h5",
    library_id="Arbitrary_Sample_ID"
)

#Writes AnnData object into .h5ad format
adata.write_h5ad("processed.h5ad")
```
This will create the .h5ad file in your current working directory! This will be used as input for Cell2Location.

## Step 2: Running Cell2Location
### 2.1: Prepared script on your end
From the pancreatic project, you may already have a working script! The following instructions will teach you how to run this preexisting script as a batch job through the supercomputer. HOWEVER, YOU MAY STILL RUN CELL2LOCATION HOWEVER YOU ARE COMFORTABLE!!!

Make a copy of the following template to create the batch job file: [click this link](https://docs.google.com/document/d/1ZcUjWUp-G6wvQJUmQyNxw2hedvN62fNrn5MarBcxf_8/edit?usp=sharing). The contents are posted below:
```
#!/bin/bash
#SBATCH --job-name=Cell2Location
#SBATCH --output=output.txt
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH --mail-type=all
#SBATCH -A CSD670
#SBATCH --mem=100G

module load singularitypro

echo "Job started"

SIF=/expanse/lustre/projects/csd670/akao1/cell2location/cell2loc.sif
DATA_DIR=/PATH/TO/PROJECT/DIRECTORY/THAT/YOU/MADE
SCRIPT=$DATA_DIR/YOUR_CELL2LOCATION_SCRIPT.py

singularity exec --bind $DATA_DIR:$DATA_DIR $SIF which python
singularity exec --bind $DATA_DIR:$DATA_DIR $SIF python --version
singularity exec --bind $DATA_DIR:$DATA_DIR $SIF python $SCRIPT

echo "Job finished"
```
There are several areas that you will need to edit yourself:
* #SBATCH --mail-user (set this to your own email)
* DATA_DIR= (set this to the project folder/directory that you set up)
* SCRIPT= (set this to $DATA_DIR/the_name_of_your_R_script.R)

Once you have done the edits, save this as a .txt file and transfer it to your supercomputer project directory using the SFTP tab on termius (or whatever preferred method)

Important note: since we unfortunately don't have GPU allocations, we will be resorting to running Cell2Location on CPU. This can be very slow, so we want to allow our batch job to access 8 CPUs at once to speed things up. To your python script, add "torch.set_num_threads(8)" as a line under all the imports.

Then, you will need to prepare all your files in the supercomputer workspace. You will need the following:
* inf_aver.csv file with output from Cell2Location trained on reference scRNA-seq data (I have created this and shared it with you)
* .h5ad file for your spatial transcriptomics sample
* Cell2Location python script (.py file)
* SBATCH input (.txt file)

Place these all in the same project directory. You may follow the same guidelines for creating/choosing a project directory from the MSI Annotations instructions.

Finally, to submit the batch job via command line/terminal:
```
cd /path/to/project/directory
dos2unix whatever_name_sbatch_script.txt
sbatch whatever_name_sbatch_script.txt
```

### 2.2: No prepared script
The following code chunk should be slightly edited (details below) and pasted into a .py file. It contains only the code for cell type abundance inference, as I already did the model training on the reference scRNA-seq dataset.
```
import scanpy as sc
import numpy as np
import pandas as pd
import cell2location
import torch
torch.set_num_threads(8) #pytorch is for deep learning, cell2location calls it during model training. default pytorch may use 1 cpu core so we tell it to use 8

#tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

print("Step1 read data", flush=True)
#1 Read data
inf_aver = pd.read_csv( #trained model posteriors from previous step
    "/expanse/lustre/scratch/akao1/temp_project/STCRC/ref_model/inf_aver.csv",
    index_col=0 #so that gene names read as rownames instead of own col
)
adata_vis = sc.read_h5ad("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563.h5ad")
print(adata_vis.var.columns.tolist(), flush=True) 

#2 filter out mitochondrial genes ("technical artifacts")
adata_vis.var_names_make_unique() #some genes may have non unique names so make them unique
print("Step2 filter MT genes", flush=True)
# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


print("Step3 find shared genes", flush=True)
#3 find shared genes between reference and visium sample
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

print("Step4 prepare visium for model", flush=True)
#4 prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)

print("Step5 create model", flush=True)
#5 create model
mod = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=8,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=200
)

print("Step6 train model", flush=True)
#6 train model
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1
         )
#saving the model data at designated directory
mod.save("/expanse/lustre/scratch/akao1/temp_project/STCRC/sample1_model", overwrite=True)
print("model saved in directory", flush=True)

print("Step export h5ad", flush=True)
#7 export posterior and add onto adata_vis object
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)
#saving anndata file with results (overwrites originally loaded h5ad)
adata_vis.write_h5ad("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_deconvolved.h5ad") #writes the anndata into h5ad
print("h5ad file with model results saved in directory", flush=True)

print("Step8 extract cell type abundances", flush=True)
#8 extract mean cell abundance per spot as dataframe
abundance_df = adata_vis.obsm['means_cell_abundance_w_sf']#takes posterior values from model
abundance_df = pd.DataFrame( #stores values in pandas dataframe
    abundance_df,
    index=adata_vis.obs_names,
    columns=adata_vis.uns['mod']['factor_names']
)

abundance_df_5 = adata_vis.obsm['q05_cell_abundance_w_sf']#takes 5% quantile posterior values from model
abundance_df_5 = pd.DataFrame( #stores values in pandas dataframe
    abundance_df_5,
    index=adata_vis.obs_names,
    columns=adata_vis.uns['mod']['factor_names']
)

# save to csv #saves this dataframe to a csv
abundance_df.to_csv("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_abundance.csv")
abundance_df_5.to_csv("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_conservative_abundance.csv")
print(f"Abundance table shape: {abundance_df.shape}", flush=True)
print("Done", flush=True)
```
For the following lines near the beginning, you need to replace the paths to the inf_aver.csv and .h5ad file with the path to your own project directory where these files should be stored.
```
inf_aver = pd.read_csv( #trained model posteriors from previous step
    "/expanse/lustre/scratch/akao1/temp_project/STCRC/ref_model/inf_aver.csv",
    index_col=0 #so that gene names read as rownames instead of own col
)
adata_vis = sc.read_h5ad("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563.h5ad")
```
For the followng lines near the end, you need to replace the paths with the desired location for the outputs
```
#saving anndata file with results (overwrites originally loaded h5ad)
adata_vis.write_h5ad("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_deconvolved.h5ad") #writes the anndata into h5ad


# save to csv #saves this dataframe to a csv
abundance_df.to_csv("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_abundance.csv")
abundance_df_5.to_csv("/expanse/lustre/scratch/akao1/temp_project/STCRC/GSM8703563_conservative_abundance.csv")
print(f"Abundance table shape: {abundance_df.shape}", flush=True)
print("Done", flush=True)
```
The final output of interest should be two .csv files with mean cell abundance estimates and 5% quantile cell abundance estimates (absolute/minimum cell abundance estimates, according to the cell2location website). You may also see a line that stores the .h5ad file in the output

If you are curious about the function of each line of code, it is adapted from the cell2location vignette [here](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html).

I hope this was helpful! Please let me know at any time if any of the code doesn't work or if there is anything to be troubleshooted :)
