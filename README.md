![](images/ikmb_bfx_logo.png)

# LaMeta

## UNDER CONSTRUCTION

THE README IS OUTDATED!

## Overview

This pipeline takes metagenomic (paired-end) short-read  as input, as generated
by Illumina sequencing. From this data, the pipeline aims to assemble high-quality
single genomes.

1. All samples are separately quality controlled to remove Illumina library adaptors, low quality sequences and sequence ends, and possible host-genetic (human) and spike in contamination (PhiX). The cleaned reads are used for the subsequently following steps.
2. For all samples a separate metagenomic assembly is performed using Spades in metagenomic mode. The sequences are then mapped back to the resulting scaffolds and binned using the MaxBin2 software.
3. Additionally, a co-assembly with Megahit is performed. Using a groupfile it is possible to split samples into separate groups for this co-assembly. Again the resulting contigs are used as reference for backmapping, followed by two separate binning approaches using MaxBin2 and Metabat2, which for this approach now can also incorporate across-sample abundance differences for the binning procedure.
4. The resulting bins from the single-sample and subgroup co-assembly approaches are finally dereplicated using the dRep package, to achieve the highest-possible quality of single-genome bins combined with low redundancy.
5. All samples are again mapped to the final resulting bins to estimate bin abundance.


## Executing the pipeline

To execute the pipeline with all default settings, do:

`nextflow -c nextflow.config run main.nf --folder /path/to/folder`

Please note: the output will be written where the pipeline is executed, NOT where the input files are located.

### Optional parameters

The pipeline uses several parameters to fine-tune the various pipeline stages. Some of these can be modified during pipeline execution:

## Reading the output

## Dependencies and versions
 * [BBMap (v.37.88)](https://sourceforge.net/projects/bbmap/): QC, Mapping to contigs.
 * [Megahit (v1.1.2)](https://github.com/voutcn/megahit): Groupwise Co-Assemblies.
 * [Spades (v.3.9.0)](http://bioinf.spbau.ru/en/spades): Single-Sample Assemblies.
 * [Samtools (v.1.5)](http://www.htslib.org): Conversion of SAM to BAM.
 * [MaxBin2 (v.2.2.4)](https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html): Binning of Contigs.
 * [Metabat2 (v.2.12.1)](https://bitbucket.org/berkeleylab/metabat/overview): Binning of Contigs.
 * [dRep (v.2.0.5)](https://github.com/MrOlm/drep): Evaluation of Binned Contigs.
 * [CheckM (v.1.0.11)](https://ecogenomics.github.io/CheckM/): Used by dRep.

Many of these tools have additional dependencies that are not listed here. If the tool works properly on its own, these are likely satisfied.

Toe remove human host/lab contamination the database has to be prepared as described [here](http://seqanswers.com/forums/showthread.php?t=42552).

`dRep` and `CheckM` used different versions of Python2 and Python3. Please follow the [instructions](http://drep.readthedocs.io) provided on the dRep website to solve this issue using `pyenv`.

 The version numbers are the software versions used in development/testing of the pipeline.

## Code structure

This pipeline is comprised of the following components:

* main.nf (the actual workflow definition)
* nextflow.config (the top-level configuration file with generic options)
* config/rzcluster.config (the RZcluster specific configuration options)
* README.md (this file)

## Credit

This pipeline was developed by M. Rühlemann.
