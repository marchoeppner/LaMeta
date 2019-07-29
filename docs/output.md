![](../images/ikmb_bfx_logo.png)

# Interpreting the output

## Folder structure

All outputs are written to the folder specified by `--outdir`. By default, this is called "outputs". 
Inside, you will find the following structure:

`CoAssembly  MultiQC  pipeline_info  refs  Samples`

The two folders "MultiQC" and "pipeline_info" contain statistics and general information about the pipeline run, such as software version numbers
and execution run time.

The folder "Samples" contains sample-specific subfolders with within-sample assemblies and qc data. 

The folder CoAssembly containts the join assembly of reads and relative abundances within a specified group. This is tapically the most important result and can
form the basis for downstream tasks such as functional annotation etc. 
