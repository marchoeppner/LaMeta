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

SPADES=file(params.spades)
SPADES_kmers=params.spades_kmers

SAMTOOLS=file(params.samtools)

JGISUM=file(params.jgisum)
METABAT=file(params.metabat)
CHECKM=file(params.checkm)

FOLDER=file(params.folder)

mode = 'skipdecon'

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

  if( mode == 'skipdecon' )
    """
    mv $left_trimmed $left_decon
    mv $right_trimmed $right_decon
    mv $unpaired $unpaired_decon
    grep $id $GROUP | cut -f 2 | tr -d '\n'
    """

  else
    """
    module load Java/1.8.0
    module load BBMap/37.88
    ${BBWRAP} minratio=0.9 threads=${task.cpus} maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 in=${left_trimmed},${unpaired} in2=${right_trimmed},NULL path=${BBWRAP_decon_phixref} outu1=${tmp_left_phix} outu2=${tmp_right_phix} outu=${tmp_unpaired_phix}
    ${BBWRAP} minratio=0.9 threads=${task.cpus} maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14 in=${tmp_left_phix},${tmp_unpaired_phix} in2=${tmp_right_phix},NULL path=${BBWRAP_decon_hsref} outu1=${left_decon} outu2=${right_decon} outu=${unpaired_decon}
    rm tmp*
    grep $id $GROUP | cut -f 2 | tr -d '\n'
    """

}

outputDecon.into{inputSpades; inputSpadesBackmap}

process runSpades {
  cpus 20
  memory 240.GB

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades"

  input:
  set id, file(left_decon), file(right_decon), file(unpaired_decon) from inputSpades

  output:
  set id, file(outcontigs) into outputSpades

  script:
  outcontigs = id + ".spades_contigs.fasta"

  """
  module load Spades/3.9.0
  $SPADES --meta --pe1-1 $left_decon --pe1-2 $right_decon --pe1-s $unpaired_decon -k $SPADES_kmers -o spades_out -t ${task.cpus}
  mv spades_out/scaffolds.fasta $outcontigs
  """
}

outputSpades.into{inputSpadesBackmapContigs; inputSpadesMetabat}
inputSpadesBackmap.join(inputSpadesBackmapContigs).set { inputSpadesBackmapWithContigs}

process runSpadesBackmap {
  cpus 5
  memory 60.GB

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades"

  input:
  set id, file(left_decon), file(right_decon), file(unpaired_decon), file(spadescontigs) from inputSpadesBackmapWithContigs

  output:
  set id, file(outdepth) into outputSpadesBackmap

  script:
  outdepth = id + ".depth.txt"

  """
  module load Java/1.8.0
  module load BBMap/37.88
  module load Samtools
  ${BBWRAP} -Xmx60g in=$left_decon,$unpaired_decon in2=$right_decon,NULL ref=$outcontigs t=20 out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
  $SAMTOOLS view -u tmp_sam.gz | samtools sort -m 54G -@ 3 -o tmp_final.bam
  $JGISUM --outputDepth $outdepth tmp_final.bam
  rm tmp*
  """
}

inputSpadesMetabat.join(outputSpadesBackmap).set { inputMetabat}

process runMetabat {
  cpus 5
  memory 60.GB

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Metabat"

  input:
  set id, file(spadescontigs), file(depthfile) from inputMetabat

  output:
  set id, file(binfolder), file(checkmout) into outputMetabatSamples

  script:
  binfolder = "metabat_bins"
  checkmout = "checkm_out"

  """
  $METABAT -i $spadescontigs -a $depthfile -o $binfolder/bin
  module load Python/2.7.10
  module load Prodigal/2.6.2
  module load Pplacer/1.1
  mkdir $checkmout
  $CHECKM lineage_wf -f $checkmout/CheckM.txt -t 20 -x fa $binfolder/ $checkmout/SCG
  """
}


inputCoAssembly.groupTuple().into{ inputCoAssemblyByGroup; inputBackmapCoassembly }

inputBackmapCoassembly.transpose() .set { inputBackmapCoassemblyT }

process runCoAssembly {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}"

  input:
  set group, set(id), file(left_decon), file(right_decon), file(unpaired_decon), file(megahitlog) from inputCoAssemblyByGroup

  output:
  set group, file(outcontigs), file(megahitlog) into outCoAssembly

  script:
  outcontigs = group + ".final_contigs.fasta"
  megahitlog = group + ".megahit.log"

  template "$TEMPLATEDIR/megahit_coassembly.sh"
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
