![](../images/ikmb_bfx_logo.png)

# Usage information

ATTENTION: The information below is focused on the medcluster hosted by the developers. Other systems will need to be [configured](installation.md) first (see installation instructions).

## Executing the pipeline

### On the Kiel Medcluster

`nextflow run marchoeppner/LaMeta --reads '/path/to/*_R{1,2}_001.fastq.gz' --host_index /path/to/index --groupfile /path/to/groupfile`

By default, the pipeline will use the human genome as host index, so that option can be skipped if data was generated from human hosts. 

### On other systems

To execute the pipeline with all default settings, do:

`nextflow -c /path/to/nextflow.config run marchoeppner/LaMeta --reads '/path/to/*_R{1,2}_001.fastq.gz' --host_index /path/to/index --groupfile /path/to/groupfile`

### `--host_index`

The host_index is a folder containing a BBMap index directory ("ref"). 

### `--fasta`

Alternatively, the pipeline can also create a host index for you, if you provide a FASTA file.

For details on how to prepare your reference database, please also see [here](http://seqanswers.com/forums/showthread.php?t=42552).

### `--adapters` 

Path to a (gzipped) FASTA file containing sequencing adapters to be removed (default: Nextera, built-in).

### `--run_name` 

Give this pipeline run a useful name

### `--email` 

If you provide an Email address and your compute system has a configured mail demon, the pipeline can send a report upon completion.

