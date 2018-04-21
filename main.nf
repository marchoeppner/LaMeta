/*

LaMeta assembly and annotation pipeline by M. Rühlemann

*/

VERSION = "0.2"

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
GROUP=file(params.groupfile)

READMINLEN = params.readminlen

BBWRAP = file(params.bbwrap)
BBDUK = file(params.bbduk)
BBMERGE = file(params.bbmerge)
ADAPTERS = file(params.adapters)
HSREF = file(params.hsref)

MEGAHIT=file(params.megahit)

SPADES=file(params.spades)
SPADES_kmers=params.spades_kmers

SAMTOOLS=file(params.samtools)

JGISUM=file(params.jgisum)
METABAT=file(params.metabat)
CHECKM=file(params.checkm)
MAXBIN=file(params.maxbin)

PYENV3=file(params.pyenv3)
PYENV2=file(params.pyenv2)
DREP=file(params.drep)

FOLDER=file(params.folder)

startfrom = params.startfrom

/*
Channel is created from the files in the input folder given by --folder.
*/

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .ifEmpty { exit 1, "Could not find a matching input file" }
  .set { inputQC }

/*
Mapping agains PhiX and Host genome (defaul:human). Mapped reads/read-pairs (also discordantly)
are discarded.
*/
process runQC {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Decon", mode: 'copy'

  input:
  set id, file(left),file(right) from inputQC

  output:
  set id,file(left_decon),file(right_decon),file(unpaired_decon) into outputQC
  set stdout,id,file(left_decon),file(right_decon),file(unpaired_decon) into inputCoAssembly

  script:

  left_trimmed = "tmp_" + id + "_R1.trimmed.fastq.gz"
  right_trimmed = "tmp_" + id + "_R2.trimmed.fastq.gz"
  unpaired_trimmed = "tmp_" + id + "_RU.trimmed.fastq.gz"

  left_nophix = "tmp_" + id + "_R1.nophix.fastq.gz"
  right_nophix = "tmp_" + id + "_R2.nophix.fastq.gz"
  unpaired_nophix = "tmp_" + id + "_RU.nophix.fastq.gz"

  left_decon = "tmp_" + id + "_R1.decon.fastq.gz"
  right_decon = "tmp_" + id + "_R2.decon.fastq.gz"
  unpaired_decon = "tmp_" + id + "_RU.decon.fastq.gz"

  merged = "tmp_" + id + "_RU.merged.fastq.gz"
  left_clean = id + "_R1.clean.fastq.gz"
  right_clean = id + "_R2.clean.fastq.gz"
  unpaired_clean = id + "_RU.clean.fastq.gz"

  if( startfrom > 0 )
    """
    mv $left_trimmed $left_decon
    mv $right_trimmed $right_decon
    mv $unpaired $unpaired_decon
    grep $id $GROUP | cut -f 2 | tr -d '\n'
    """
  else
    """
    ${BBDUK} threads=${task.cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=${unpaired_trimmed} ref=${ADAPTERS} ktrim=r k=23 mink=11 hdist=1 minlength=${READMINLEN} tpe tbo
    ${BBDUK} threads=${task.cpus} in=${left_trimmed} in2=${right_trimmed} k=31 ref=artifacts,phix ordered cardinality out1=${left_nophix} out2=${right_nophix} minlength=${READMINLEN}
    ${BBDUK} threads=${task.cpus} in=${unpaired_trimmed}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_nophix} minlength=${READMINLEN}
    ${BBWRAP} -Xmx23g threads=${task.cpus} minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=20 minlength=${READMINLEN} in=${left_nophix},${unpaired_nophix} in2=${right_nophix},NULL path=${HSREF} outu1=${left_decon} outu2=${right_decon} outu=${unpaired_decon}
    ${BBMERGE} threads=${task.cpus} in1=${left_decon} in2=${right_decon} out=${merged} outu1=${left_clean} outu2=${right_clean} mininsert=${READMINLEN}
    zcat ${merged} ${unpaired_nophix} | gzip -c > ${unpaired_clean}
    rm tmp*
    grep $id $GROUP | cut -f 2 | tr -d '\n'
    """
}

/*
Co assembly within the groups given in groupfile.
*/
inputCoAssembly.groupTuple().into{ inputCoAssemblyByGroup; inputBackmapMegahit }
process runCoAssembly {

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

  input:
  set group, id, file(left_clean), file(right_clean), file(unpaired_clean) from inputCoAssemblyByGroup

  output:
  set group, file(outcontigs), file(megahitlog) into outCoAssembly

  script:
  outcontigs = group + ".final_contigs.fasta"
  megahitlog = group + ".megahit.log"

  if( startfrom > 1 )
  """
  cp ${OUTDIR}/CoAssembly/${group}/$outcontigs $outcontigs
  cp ${OUTDIR}/CoAssembly/${group}/$megahitlog $megahitlog
  """

  else
  """
  echo $left_clean > out
  echo $right_clean >> out
  echo $unpaired_clean >> out

  awk '
  {
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' out > tmp1

  awk '{printf " -1 " \$1 " -2 " \$2 " -r " \$3}' tmp1 > tmp

  $MEGAHIT \$(cat tmp | tr -d '\n') --num-cpu-threads ${task.cpus} --presets meta-large -o megahit_out --mem-flag 2 --verbose
  mv megahit_out/final.contigs.fa $outcontigs
  mv megahit_out/log $megahitlog
  """
}

/*
Single-samples metagenome assembly with spades
*/
outputQC.into{inputSpades; inputSpadesBackmap}
process runSpades {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

  input:
  set id, file(left_clean), file(right_clean), file(unpaired_clean) from inputSpades

  output:
  set id, file(outcontigs) into outputSpades

  script:
  outcontigs = id + ".spades_contigs.fasta"

  if( startfrom > 1 )
  """
  cp ${OUTDIR}/Samples/${id}/Spades/$outcontigs $outcontigs
  """

  else
  """
  module load Spades/3.9.0
  $SPADES --meta --pe1-1 $left_clean --pe1-2 $right_clean --pe1-s $unpaired_clean -k $SPADES_kmers -o spades_out -t ${task.cpus}
  mv spades_out/scaffolds.fasta $outcontigs
  """
}


/*
Backmapping to spades assembly and contig abundance estimation
*/
outputSpades.into{inputSpadesBackmapContigs; inputSpadesMaxbin}
inputSpadesBackmap.join(inputSpadesBackmapContigs).set { inputSpadesBackmapWithContigs}

process runSpadesBackmap {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

  input:
  set id, file(left_clean), file(right_clean), file(unpaired_clean), file(spadescontigs) from inputSpadesBackmapWithContigs

  output:
  set id, file(outdepth) into outputSpadesBackmap

  script:
  outdepth = id + ".depth.txt"
  if( startfrom > 2 )
  """
  cp ${OUTDIR}/Samples/${id}/Spades/$outdepth $outdepth
  """

  else
  """
  module load Java/1.8.0
  module load BBMap/37.88
  module load Samtools/1.5
  ${BBWRAP} -Xmx60g in=$left_clean,$unpaired_clean in2=$right_clean,NULL ref=$spadescontigs t=${task.cpus} out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
  $SAMTOOLS view -u tmp_sam.gz | $SAMTOOLS sort -m 54G -@ 3 -o tmp_final.bam
  $JGISUM --outputDepth $outdepth tmp_final.bam
  rm tmp*
  """
}

/*
Single-sample binning with Maxbin2
*/
inputSpadesMaxbin.join(outputSpadesBackmap).set {inputMetabat; inputMaxbin}

process runMaxbin {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Maxbin", mode: 'copy'

  input:
  set id, file(spadescontigs), file(depthfile) from inputMaxbin

  output:
  set id, file(binfolder) into outputMaxbinSamples

  script:
  binfolder = id + "_maxbin_bins"
  if( startfrom > 2 )
  """
  cp -r ${OUTDIR}/Samples/${id}/Maxbin/$binfolder $binfolder
  """

  else
  """
  tail -n+2 $depthfile | cut -f 1,3 > maxbin.cov
  mkdir $binfolder
  $MAXBIN -contig $spadescontigs -abund maxbin.cov -out $binfolder/${id}.bin -thread ${task.cpus}
  mv $binfolder/${id}.bin.noclass $binfolder/${id}.bin.noclass.fasta
  """
}

process runMetabat {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Metabat", mode: 'copy'

  input:
  set id, file(spadescontigs), file(depthfile) from inputMetabat

  output:
  set id, file(binfolder) into outputMetabatSamples

  script:
  binfolder = id + "_metabat_bins"

  if( startfrom > 2 )
  """
  cp -r ${OUTDIR}/Samples/${id}/Metabat/$binfolder $binfolder
  """

  else
  """
  mkdir $binfolder
  $METABAT -i $spadescontigs -a $depthfile -o $binfolder/${id}.metabat.bin -t ${task.cpus}
  """
}

/*
Backmapping to Megahit groupwise co-assembly
*/
outCoAssembly.into{ inputContigsBackmapMegahit; inputContigsMegahitMaxbin; inputContigsMegahitMetabat }
inputBackmapMegahit.transpose().combine(inputContigsBackmapMegahit, by: 0) .set { inputBackmapCoassemblyT }

process runCoassemblyBackmap {

  tag "${group}-${id}"
  publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

  input:
  set group, id, file(left_clean), file(right_clean), file(unpaired_clean), file(megahitcontigs), file(megahitlog) from inputBackmapCoassemblyT

  output:
  set group, file(bamout) into outMegahitBackmap

  script:
  bamout = id + ".megahit.final.bam"

  if( startfrom > 2 )
  """
  cp ${OUTDIR}/CoAssembly/${group}/$bamout $bamout
  """

  else
  """
  module load Java/1.8.0
  module load BBMap/37.88
  module load Samtools/1.5
  ${BBWRAP} -Xmx60g in=$left_clean,$unpaired_clean in2=$right_clean,NULL ref=$megahitcontigs t=${task.cpus} out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
  $SAMTOOLS view -u tmp_sam.gz | $SAMTOOLS sort -m 54G -@ 3 -o $bamout
  rm tmp*
  """
}

/*
Contig abundance estimation for co-assemblies
*/
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

  if( startfrom > 2 )
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

/*
Co-assembly binning with Maxbin2
*/
coassemblyAbufolder.join(inputContigsMegahitMaxbin).set{ inputMegahitMaxbin }
process runMegahitMaxbin {
  cpus 20
  memory 240.GB

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}/Maxbin", mode: 'copy'

  input:
  set group, file(inputfolder), file(megahitcontigs) from inputMegahitMaxbin

  output:
  set group, file(binfolder) into outputMegahitMaxbin

  script:
  binfolder = group + "_maxbin_bins"

  if( startfrom > 2 )
  """
  cp -r ${OUTDIR}/CoAssembly/${group}/Maxbin/$binfolder $binfolder
  """

  else
  """
  ls ${inputfolder}/*.out > abufiles.txt
  mkdir $binfolder
  $MAXBIN -contig $megahitcontigs -abund_list abufiles.txt -out $binfolder/${group}.maxbin.bin -thread ${task.cpus}
  mv $binfolder/${group}.maxbin.bin.noclass $binfolder/${group}.maxbin.bin.noclass.fasta
  """
}

/*
Co-assembly binning with Metabat
*/
coassemblyDepth.join(inputContigsMegahitMetabat).set{ inputMegahitMetabat }
process runMegahitMetabat {

  tag "${group}"
  publishDir "${OUTDIR}/CoAssembly/${group}/Metabat", mode: 'copy'

  input:
  set group, file(inputdepth), file(megahitcontigs) from inputMegahitMetabat

  output:
  set group, file(binfolder) into outputMegahitMetabat

  script:
  binfolder = group + "_metabat_bins"

  if( startfrom > 2 )
  """
  cp -r ${OUTDIR}/CoAssembly/${group}/Metabat/$binfolder $binfolder
  """

  else

  """
  mkdir $binfolder
  $METABAT -i $megahitcontigs -a $inputdepth -o $binfolder/${group}.metabat.bin -t ${task.cpus}
  """
}

/*
Dereplication of all bins from single-sample and groupwise co-assemblies
*/
source = Channel.create()
allbinfolders = Channel.create()
outputMaxbinSamples.mix(outputMetabatSamples,outputMegahitMetabat,outputMegahitMaxbin).separate( source, allbinfolders )
allbinfolders.collect().into { inputDrep; test }
test.println()

process runDrep {

  tag "allbins"
  publishDir "${OUTDIR}/Final/dRep", mode: 'copy'

  input:
  file binfolder from inputDrep.collect()

  output:
  file outfolder into outputDrep

  script:
  outfolder = "dRep_out"


  """
  mkdir allbins
  for binf in ${binfolder}; do
  cp \$binf/*.fa* allbins
  done

  pyenv local $PYENV3 $PYENV2
  $DREP bonus testDir --check_dependencies
  $DREP dereplicate $outfolder -g allbins/*.fa* -p ${task.cpus}
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
