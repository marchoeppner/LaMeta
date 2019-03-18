![](../images/ikmb_bfx_logo.png)

# Usage information

## Executing the pipeline

To execute the pipeline with all default settings, do:

`nextflow -c nextflow.config run main.nf --reads '/path/to/*_R{1,2}_001.fastq.gz' --host_index /path/to/index -profile your_profile`

Where `host_index` is a folder containing a BBMap index directory ("ref") and `your_profile` refers to the cluster profile you have [configured](installation.md). 

Alternatively, the pipeline can also create an index for you, if you provide a FASTA file:

`nextflow -c nextflow.config run main.nf --reads '/path/to/*_R{1,2}_001.fastq.gz' --host /path/to/genome.fa`

For details on how to prepare your reference database, please also see [here](http://seqanswers.com/forums/showthread.php?t=42552).

### Execution profile

Make sure you read the [Installation and Configuration](installation.md) and specify your profile using `-profile`

### Optional parameters

The pipeline uses several parameters to fine-tune the various pipeline stages. Some of these can be modified during pipeline execution.

`--adapters` Path to a (gzipped) FASTA file containing sequencing adapters to be removed (default: Nextera, built-in).

`--run_name` Give this pipeline run a useful name

`--email` If you provide an Email address and your compute system has a configured mail demon, the pipeline can send a report upon completion.

For a list of all parameters, see:

`nextflow -c nextflow.config run main.nf --help`
