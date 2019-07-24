![](../images/ikmb_bfx_logo.png)

# Installation and configuration

## Summary

in order to run LaMeta on your system, you have to:

- download Nextflow

- make sure either Conda, Singularity or Docker are available on your system

- fill out a config file that tells the pipeline how to talk to your cluster and where to find some of the external references

## Pre-requisties

LaMeta uses the Nextflow pipeline framework for job distribution on a compute cluster. In addition, software is provisioned either through conda 
or a Docker container. 

### Hardware requirements

LaMeta was written to run on compute clusters running Linux. Mac and Windows are not supported, but you may be able to run the pipeline on your 
Linux Desktop machine - assuming it has sufficient memory (>> 100GB Ram for typical samples).

#### Memory

Most of LaMeta runs on typical HPC compute nodes with 128GB Ram (or less). 
However, the assembly stages - especially the joint assembly of all samples in a given group - can consume large amounts of memory. 
A Node with at least 256GB Ram is recommended. 

### Get Nextflow.

Nextflow can be downloaded [here](https://github.com/nextflow-io/nextflow/releases).
LaMeta has been tested with Nextflow version 19.02 - older versions may work too. 

Note that Nextflow requires Java version 1.8+.

### Install Conda (Option 1)

Conda is an easy way to install a wide range of software packages. For more information, see 
[here](https://conda.io/projects/conda/en/latest/user-guide/overview.html)

Any conda version 4.0+ should work, but try to get the latest version if possible from [here](https://www.anaconda.com/distribution/)

If conda is setup and chosen as privisioning option, LaMeta will make sure that all software is made available during start-up. While this is convenient, 
it is also quite slow. We therefore recommend using Singularity instead (see below).

### Install Singularity (Option 2)

[Singularity](https://www.sylabs.io/singularity/) is a container framework, similar to Docker. Unlike Docker, it can be run on shared compute clusters 
where the user may not have root proviliges. 

Should Singularity not yet be available on your compute cluster, please contact your administrators. 

If Singularity is chosen for software provisioning, LaMeta will pull the container from 
[Dockerhub](https://cloud.docker.com/repository/docker/mhoeppner/lameta). This usually takes less than 30 seconds. 

## Install Docker (Option 3)

[Docker](https://www.docker.com/) is a widely used container-framework. However, it is usually not available on distributed compute clusters.
However, if your setup offers docker (such as in virtualized cloud environment like AWS or OpenStack), you may of course also use Docker. 

### Setting up your own configuration

Translating LaMeta to your compute environment requires that you set some options and fill out a specific config file. 

#### Config files

LaMeta uses three (or more) config files.

`base.config` This is always needed and specifies how much resources each stage of the pipeline requires. 

One of:

`conda.config` This is needed if you wish to use conda for provisioning of software. 

`singularity.config` This is needed if you wish to use Singularity for provisioning of software.

`docker.config` This is needed if your wish to use Eocker for provisioning of software. 

Your cluster:

`your_cluster.config` This file should contain information about your cluster queue. Please see [here](https://www.nextflow.io/docs/latest/executor.html) 
for some information on how to do that or have a look at the other included files. You will need to specify how jobs are to be disributed - either using a
queuing system like Slurm, LSF or SGE - or maybe you want to run this on e Amazon cloud (AWS, which Nextflow supports but has not been tested by us for 
this pipeline). 

An important section in your profile is, in addition to specifying your queuing system, to tell Nextflow how many cores and RAM each node has and what 
the maximum allowed walltime for a job should be. 

```
params {
  max_memory = 128.GB
  max_cpus = 16
  max_time = 240.h
}
```
In this example, each compute note has 128GB of RAM, 16 CPU cores and jobs are allowed to run for 240 hours before they are terminated.

When using either Singularity or Docker, you will have to specify any local paths that should be mounted inside the container and specifically you will 
have to bind the GTDB database into the countainer using the below statement. 

```
        singularity {
                runOptions = " -B /path/to/gtdb/:/refdata/"
        }
```

#### Adding your own profile

Once you have your own cluster config file setup and included under `/conf`, you should add your own execution profile to `nextflow.config`.

```
profiles {

	your_profile {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
		includeConfig 'conf/conda.config'
	}
```


This example assumes that you wish to use conda for provisioning. If you want to use Singularity or Docker, replace conda.config with either of the other 
two config files. 

To use this profile, specify it when starting LaMeta with `-profile your_profile`.

### Reference data

LaMeta requires two larger reference data sets to run and which are not provisioned automatically. Note that all the options listed below can also be set 
from the command line when starting the pipeline (params.some_parameter becomes `--some_parameter`). 

#### GTDB-TK Reference DB

Please download the database [here](https://data.ace.uq.edu.au/public/gtdbtk/release_86/gtdbtk.r86_v2_data.tar.gz)

Once the database has been downloaded and de-compressed (tar -xvf gtdbtk.r86_v2_data.tar.gz), add the full path to the directory to your config file 

`params.gtdb = /path/to/gtdb/`

Make sure to include the terminal "/"; GTDK-TK is unfortunately not smart enought to automatically add it and will fail with a "file not found" error 
otherwise. 

#### Host reference genome

This is optional and can also be provided during pipeline execution; depending on whether you have your reference data in FASTA format or already as a 
BBMap index. 

For instructions on how to prepare a reference dataset, please see [here](http://seqanswers.com/forums/showthread.php?t=42552)

If you wish to provide a genome sequence, you can permanently set it in your config file using:

`params.host = /path/to/genome.fa` 

However, since this will force the pipeline to build the BBMap index every time the pipeline is run, the better option would to build the index and use that. 

If you have a folder container a BBMap `ref`subfolder, you can do:

`params.host_index = /path/to/folder`

Note that "path/to/folder" should point to the folder container the BBMap "ref" folder, not the "ref" folder itself. 

