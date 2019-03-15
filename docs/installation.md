![](../images/ikmb_bfx_logo.png)

## Installation and configuration

### Pre-requisties

LaMeta uses the Nextflow pipeline framework for job distribution on a compute cluster. In addition, software is provisioned either through conda or a singularity container. 

### Get Nextflow.

Nextflow can be downloaded [here](https://github.com/nextflow-io/nextflow/releases).
LaMeta has been tested with Nextflow version 19.02 - older versions may work too. 

Note that Nextflow requires Java version 1.8+.

#### Install Conda

Conda is an easy way to install a wide range of software packages. For more information, see [here](https://conda.io/projects/conda/en/latest/user-guide/overview.html)

Any conda version 4.0+ should work, but try to get the latest version if possible from [here](https://www.anaconda.com/distribution/)

If conda is setup and chosen as privisioning option, LaMeta will make sure that all software is made available during start-up. While this is convenient, it is also quite slow. If this is a concern, consider using Singularity instead (see below).

#### Install Singularity

[Singularity](https://www.sylabs.io/singularity/) is a container framework, similar to Docker. Unlike Docker, it can be run on shared compute clusters where the user may not have root proviliges. 

Should Singularity not yet be available on your compute cluster, please contact your administrators. 

If Singularity is chosen for software provisioning, LaMeta will pull the container from [Singularity hub](https://singularity-hub.org/). This usually takes less than 30 seconds. 

### Setting up your own configuration

Translating LaMeta to your compute environment requires that you set some options and fill out a specific config file. 

#### Config files

LaMeta uses three environment-specific config files.

`base.config` This is always needed and specifies how much resources each stage of the pipeline requires. 

One of:

`conda.config` This is needed if you wish to use conda for provisioning of software. 

`singularity.config` This is needed if you wish to use Singularity for provisioning of software.

Your cluster:

`your_cluster.config` This file should contain information about your cluster queue. Please see [here](https://www.nextflow.io/docs/latest/executor.html) for some information on how to do that or have a look at the other included files. 

#### Adding your own profile

Once you have your own cluster config file setup and included under `/conf`, you should add your own execution profile to `nextflow.config`.

```
profiles {

	your_profile {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
		includeConfig 'conf/conda.config'

```


This example assumes that you wish to use conda for provisioning. 

To use this profile, specify it when starting LaMeta with `-profile your_profile`.


