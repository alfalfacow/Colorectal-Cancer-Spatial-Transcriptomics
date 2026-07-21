# Running Cell2Location
As you may be aware, [Cell2Location](https://cell2location.readthedocs.io/en/latest/index.html) is a tool for "deconvolution" of spatial transcriptomics data. While spatial transcriptomics indeed approaches near single cell resolution (with more and more advanced tools such as 10X Genomics' Atera reaching subcellular resolution), many technologices such as Visium may have more than one cell per spot profiled. Especially considering that Visium spots are 55 micrometers in diameter (while most cells are far below that), there are estimates of up to 10 cells per spot. Thus, labeling each spot as a single cell is not accurate.

Cell2Location bridges this obstacle by estimating how many of each cell type is present within a spot (cell type abundance per spot). This is done by first training Cell2Location using a reference single cell RNA sequencing dataset with cell types already annotated. Cell2Location learns the transcriptomic profile from this reference dataset and applies it to the spatial transcriptomics dataset.

For the purposes of this project, I have already trained Cell2Location on the reference datasets ([GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735) and [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)), both of which have cell types prelabeled. These were chosen due to their use in previous CRC studies using cell2location, and also because it includes cells labeled by consensus molecular subtype (CMS), another CRC-specific feature we may be interested in.

You may have already run Cell2Location for the pancreatic cancer project, so this is mostly instructions for adapting that script for use on the supercomputer. That way, you won't have to leave it running with your computer on all day. Running this on the supercomputer is optional! If you prefer, feel free to continue to run it locally (it may actually be slightly faster since it can utilize GPU, and you may already know how to do it).

## Table of Contents

1: Preparing the scripts
* 1.1: Python script
* 1.2: Batch job script

2: Running the scripts

## Step 1: Preparing the scripts
First, you w
