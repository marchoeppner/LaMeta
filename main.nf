/*

LaMeta assembly and annotation pipeline by M. Rühlemann

*/

VERSION = "0.1"

logParams(params, "nextflow_parameters.txt")

// Header log info
log.info "========================================="
log.info "LaMeta assembly and annotation pipeline v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "Author:		Malte Rühlemann"
log.info "========================================="
log.info "Starting at:		$workflow.start"

OUTDIR=file(params.outdir)

TRIMMOMATIC = file(params.trimmomatic)
TRIMMOMATIC_adapters = file(params.trimmomatic_adapters)
TRIMMOMATIC_minlen = file(params.trimmomatic_minlen)

FOLDER=file(params.folder)

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .ifEmpty { exit 1, "Could not find a matching input file" }
  .into { inputTrim }

process runTrim {

  time = { 60.m * task.attempt }
  memory = { 40.GB * task.attempt }
  cpus 20

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Trim"

  input:
  set id, file(left),file(right) from inputTrim

  output:
  set id,file(left_trimmed),file(right_trimmed),file(unpaired) into outputTrim

  script:

  left_trimmed = id + "_R1.trimmed.fastq.gz"
  right_trimmed = id + "_R2.trimmed.fastq.gz"
  left_unpaired = id + "_R1.unpaired.fastq.gz"
  right_unpaired = id + "_R2.unpaired.fastq.gz"
  unpaired = id + "_RU.trimmed.fastq.gz"

  """
  java -jar $TRIMMOMATIC PE -threads ${cpus} -trimlog $s.trimlog.txt -phred33 ${left} ${right} ${left_trimmed} ${left_unpaired} ${right_trimmed} ${right_unpaired} ILLUMINACLIP:${TRIMMOMATIC_adapters}:1:50:30:1:true MINLEN:${TRIMMOMATIC_minlen} SLIDINGWINDOW:15:20
  zcat ${left_unpaired} ${right_unpaired} > ${unpaired}
  rm ${left_unpaired} ${right_unpaired}
  """

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}
//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################
// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"
  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}
