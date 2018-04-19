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
MAXBIN=file(params.maxbin)

FOLDER=file(params.folder)

mode = 'testmode'

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .ifEmpty { exit 1, "Could not find a matching input file" }
  .set { inputTrim }

process runTrim {

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

  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/Samples/${id}/Trim/$left_trimmed $left_trimmed
  cp ${OUTDIR}/Samples/${id}/Trim/$right_trimmed $right_trimmed
  cp ${OUTDIR}/Samples/${id}/Trim/$unpaired $unpaired
  """

  else
  """
  module load Java/1.8.0
  java -jar $TRIMMOMATIC PE -threads ${task.cpus} -trimlog ${trimlog} -phred33 ${left} ${right} ${left_trimmed} ${left_unpaired} ${right_trimmed} ${right_unpaired} ILLUMINACLIP:${TRIMMOMATIC_adapters}:1:50:30:1:true MINLEN:${TRIMMOMATIC_minlen} SLIDINGWINDOW:15:20
  cat ${left_unpaired} ${right_unpaired} | gzip -c > ${unpaired}
  rm ${left_unpaired} ${right_unpaired}
  """

}

process runDecon {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Decon", mode: 'copy'

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

  if( mode == 'testmode' )
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

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

  input:
  set id, file(left_decon), file(right_decon), file(unpaired_decon) from inputSpades

  output:
  set id, file(outcontigs) into outputSpades

  script:
  outcontigs = id + ".spades_contigs.fasta"

  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/Samples/${id}/Spades/$outcontigs $outcontigs
  """

  else
  """
  module load Spades/3.9.0
  $SPADES --meta --pe1-1 $left_decon --pe1-2 $right_decon --pe1-s $unpaired_decon -k $SPADES_kmers -o spades_out -t ${task.cpus}
  mv spades_out/scaffolds.fasta $outcontigs
  """
}

outputSpades.into{inputSpadesBackmapContigs; inputSpadesMaxbin}
inputSpadesBackmap.join(inputSpadesBackmapContigs).set { inputSpadesBackmapWithContigs}

process runSpadesBackmap {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

  input:
  set id, file(left_decon), file(right_decon), file(unpaired_decon), file(spadescontigs) from inputSpadesBackmapWithContigs

  output:
  set id, file(outdepth) into outputSpadesBackmap

  script:
  outdepth = id + ".depth.txt"
  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/Samples/${id}/Spades/$outdepth $outdepth
  """

  else
  """
  module load Java/1.8.0
  module load BBMap/37.88
  module load Samtools/1.5
  ${BBWRAP} -Xmx60g in=$left_decon,$unpaired_decon in2=$right_decon,NULL ref=$spadescontigs t=${task.cpus} out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
  $SAMTOOLS view -u tmp_sam.gz | $SAMTOOLS sort -m 54G -@ 3 -o tmp_final.bam
  $JGISUM --outputDepth $outdepth tmp_final.bam
  rm tmp*
  """
}

inputSpadesMaxbin.join(outputSpadesBackmap).set { inputMaxbin}

process runMaxbin {
  cpus 5
  memory 60.GB

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Maxbin", mode: 'copy'

  input:
  set id, file(spadescontigs), file(depthfile) from inputMaxbin

  output:
  set id, file(binfolder) into outputMaxbinSamples

  script:
  binfolder = "maxbin_bins"
  if( mode == 'testmode' )
  """
  cp -r ${OUTDIR}/Samples/${id}/Maxbin/$binfolder $binfolder
  """

  else
  """
  tail -n+2 $depthfile | cut -f 1,3 > maxbin.cov
  mkdir $binfolder
  $MAXBIN -contig $spadescontigs -abund maxbin.cov -out $binfolder/${id}.bin -thread ${task.cpus}
  mv $binfolder/bin.noclass $binfolder/bin.noclass.fasta
  """
}


inputCoAssembly.groupTuple().into{ inputCoAssemblyByGroup; inputBackmapMegahit }


process runCoAssembly {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

  input:
  set group, id, file(left_decon), file(right_decon), file(unpaired_decon) from inputCoAssemblyByGroup

  output:
  set group, file(outcontigs), file(megahitlog) into outCoAssembly

  script:
  outcontigs = group + ".final_contigs.fasta"
  megahitlog = group + ".megahit.log"

  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/CoAssembly/${group}/$outcontigs $outcontigs
  cp ${OUTDIR}/CoAssembly/${group}/$megahitlog $megahitlog
  """

  else
  template "$TEMPLATEDIR/megahit_coassembly.sh"
}


outCoAssembly.into{ inputContigsBackmapMegahit; inputContigsMegahitMaxbin; inputContigsMegahitMetabat }

inputBackmapMegahit.transpose().combine(inputContigsBackmapMegahit, by: 0) .set { inputBackmapCoassemblyT }


process runCoassemblyBackmap {

  tag "${group}-${id}"
  publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

  input:
  set group, id, file(left_decon), file(right_decon), file(unpaired_decon), file(megahitcontigs), file(megahitlog) from inputBackmapCoassemblyT

  output:
  set group, file(bamout) into outMegahitBackmap

  script:
  bamout = id + ".megahit.final.bam"

  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/CoAssembly/${group}/$bamout $bamout
  """

  else
  """
  module load Java/1.8.0
  module load BBMap/37.88
  module load Samtools/1.5
  ${BBWRAP} -Xmx60g in=$left_decon,$unpaired_decon in2=$right_decon,NULL ref=$megahitcontigs t=${task.cpus} out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
  $SAMTOOLS view -u tmp_sam.gz | $SAMTOOLS sort -m 54G -@ 3 -o $bamout
  rm tmp*
  """
}

outMegahitBackmap.groupTuple().set { inputCollapseBams }

process runCollapseBams {

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

  input:
  set group, file(bams) from inputCollapseBams

  output:
  set group, file(depthfile) into coassemblyDepth
  set group, file(abufolder) into coassemblyAbufolder

  script:
  depthfile = group + "depth.txt"
  abufolder = group + "_abufiles"

  if( mode == 'testmode' )
  """
  cp ${OUTDIR}/CoAssembly/${group}/$depthfile $depthfile
  cp -r ${OUTDIR}/CoAssembly/${group}/$abufolder $abufolder
  """

  else
  """
  $JGISUM --outputDepth $depthfile $bams
  ncol=\$(head -n 1 $depthfile | awk '{print NF}')
  mkdir $abufolder
  for i in \$(seq 4 2 \$ncol); do
  name=\$(head -n 1 $depthfile | cut -f \$i | cut -d "." -f 1)
  cut -f  1,\$i $depthfile | tail -n+2 > $abufolder/\${name}.out
  done
  """
}

coassemblyAbufolder.join(inputContigsMegahitMaxbin).set{ inputMegahitMaxbin }

process runMegahitMaxbin {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}/Maxbin"

  input:
  set group, file(inputfolder), file(megahitcontigs) from inputMegahitMaxbin

  output:
  set id, file(binfolder) into outputMegahitMaxbin

  script:
  binfolder = "maxbin_bins"


  """
  ls ${inputfolder}/*.out > abufiles.txt
  mkdir $binfolder
  $MAXBIN -contig $megahitcontigs -abund_list abufiles.txt -out $binfolder/${group}.maxbin.bin -thread ${task.cpus}
  mv $binfolder/bin.noclass $binfolder/bin.noclass.fasta
  """
}

coassemblyDepth.join(inputContigsMegahitMetabat).set{ inputMegahitMetabat }

process runMegahitMetabat {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}/Metabat"

  input:
  set group, file(inputdepth), file(megahitcontigs) from inputMegahitMetabat

  output:
  set id, file(binfolder) into outputMegahitMetabat

  script:
  binfolder = "metabat_bins"


  """
  mkdir $binfolder
  $METABAT -i $megahitcontigs -a $inputdepth -o $binfolder/${group}.metabat.bin -t ${task.cpus}
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
