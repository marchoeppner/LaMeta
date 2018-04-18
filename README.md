# LaMeta

## UNDER CONSTRUCTION

## Overview

This pipeline takes metagenomic read data and assembles it to recover good quality microbial genomes, followed by annotation.
 
## Executing the pipeline

To execute the pipeline with all default settings, do:

`nextflow -c nextflow.config run main.nf --folder /path/to/folder`

Please note: the output will be written where the pipeline is executed, NOT where the input files are located. 

### Optional parameters

The pipeline uses several parameters to fine-tune the various pipeline stages. Some of these can be modified during pipeline execution:

## Reading the output

## Tools and versions

## Code structure

This pipeline is comprised of the following components:

* main.nf (the actual workflow definition)
* nextflow.config (the top-level configuration file with generic options)
* config/rzcluster.config (the RZcluster specific configuration options)
* README.md (this file)

## Credit

This pipeline was developed by M. Rühlemann.
