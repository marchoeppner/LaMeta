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
GROUP=file("groupfile.txt")

TEMPLATEDIR=file(params.templatedir)

TRIMMOMATIC = file(params.trimmomatic)
TRIMMOMATIC_adapters = file(params.trimmomatic_adapters)
TRIMMOMATIC_minlen = params.trimmomatic_minlen

BBWRAP = file(params.bbmap)
BBWRAP_decon_hsref = file(params.decon_hsref)
BBWRAP_decon_phixref = file(params.decon_phixref)

MEGAHIT=file(params.megahit)

FOLDER=file(params.folder)

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .ifEmpty { exit 1, "Could not find a matching input file" }
  .set { inputTrim }

process runTrim {

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
  left_unpaired = id + "_R1.unpaired.fastq"
  right_unpaired = id + "_R2.unpaired.fastq"
  unpaired = id + "_RU.trimmed.fastq.gz"
  trimlog = id + ".trimlog.txt"

  """
  module load Java/1.8.0
  java -jar $TRIMMOMATIC PE -threads ${task.cpus} -trimlog ${trimlog} -phred33 ${left} ${right} ${left_trimmed} ${left_unpaired} ${right_trimmed} ${right_unpaired} ILLUMINACLIP:${TRIMMOMATIC_adapters}:1:50:30:1:true MINLEN:${TRIMMOMATIC_minlen} SLIDINGWINDOW:15:20
  cat ${left_unpaired} ${right_unpaired} | gzip -c > ${unpaired}
  rm ${left_unpaired} ${right_unpaired}
  """

}

process runDecon {

  cpus 20

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Decon"

  input:
  set id, file(left_trimmed),file(right_trimmed),file(unpaired) from outputTrim

  output:
  set id,file(left_decon),file(right_decon),file(unpaired_decon) into outputDecon
  set stdout,id,file(left_decon),file(right_decon),file(unpaired_decon) into inputCoAssembly
  script:

  tmp_left_phix = "tmp_" + id + "_R1.nophix.fastq.gz"
  tmp_right_phix = "tmp_" + id + "_R2.nophix.fastq.gz"
  tmp_unpaired_phix = "tmp_" + id + "_RU.nophix.fastq.gz"
  left_decon = id + "_R1.decon.fastq.gz"
  right_decon = id + "_R2.decon.fastq.gz"
  unpaired_decon = id + "_RU.decon.fastq.gz"

  """
  module load Java/1.8.0
  module load BBMap/37.88
  ${BBWRAP} minratio=0.9 threads=${task.cpus} maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 in=${left_trimmed},${unpaired} in2=${right_trimmed},NULL path=${BBWRAP_decon_phixref} outu1=${tmp_left_phix} outu2=${tmp_right_phix} outu=${tmp_unpaired_phix}
  ${BBWRAP} minratio=0.9 threads=${task.cpus} maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 in=${tmp_left_phix},${tmp_unpaired_phix} in2=${tmp_right_phix},NULL path=${BBWRAP_decon_hsref} outu1=${left_decon} outu2=${right_decon} outu=${unpaired_decon}
  rm tmp*
  grep $id $GROUP | cut -f 2 | tr -d '\n'
  """

}

inputCoAssembly.groupTuple().set { inputCoAssemblyByGroup }

process runCoAssembly {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}"

  input:
  set group, id, file(left_decon), file(right_decon), file(unpaired_decon) from inputCoAssemblyByGroup

  output:
  set group, file(left_decon), file(right_decon), file(unpaired_decon), file(megahitlog) into outCoAssembly

  outcontigs = group + ".final_contigs.fasta"
  megahitlog = group + ".megahit.log"

  template "$TEMPLATEDIR/megahit_coassembly.sh"



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
