![](../images/ikmb_bfx_logo.png)

# LaMeta

## Overview

This pipeline takes metagenomic (paired-end) short-read data as input, as generated by Illumina sequencing. From this data, the pipeline aims to assemble high-quality single genomes.

1. All samples are separately quality controlled to remove Illumina library adaptors, low quality sequences and sequence ends, and possible host-genetic (human) and spike in contamination (PhiX). The cleaned reads are used for the subsequently following steps.
2. For all samples a separate metagenomic assembly is performed using Spades in metagenomic mode. The sequences are then mapped back to the resulting scaffolds and binned using the MaxBin2 software.
3. Additionally, a co-assembly with Megahit is performed. Using a groupfile it is possible to split samples into separate groups for this co-assembly. Again the resulting contigs are used as reference for backmapping, followed by two separate binning approaches using MaxBin2 and Metabat2, which for this approach now can also incorporate across-sample abundance differences for the binning procedure.
4. The resulting bins from the single-sample and subgroup co-assembly approaches are finally dereplicated using the dRep package, to achieve the highest-possible quality of single-genome bins combined with low redundancy.
5. All samples are again mapped to the final resulting bins to estimate bin abundance.

As LaMeta is designed to run on common distributed compute systems by using [Netxflow](https://www.nextflow.io), its throughput is limited only by the available hardware and can theoretically scale to hundreds of samples. 

## Documentation

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Credit

This pipeline was developed by M. Rühlemann and M. Höppner. 
