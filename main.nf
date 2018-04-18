/* 

16S pre-processing pipeline by M. Rühlemann

*/

VERSION = "1.1"

logParams(params, "nextflow_parameters.txt")

// Header log info
log.info "========================================="
log.info "16S pre-processing pipeline v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "Author:		Malte Rühlemann, Marc Höppner"
log.info "========================================="
log.info "Starting at:		$workflow.start"

OUTDIR=file(params.outdir)

SICKLE = file(params.sickle)
VSEARCH = file(params.vsearch)
FASTQ_QUALITY_FILTER = file(params.fastq_quality_filter)
FASTQ_TO_FASTA = file(params.fastq_to_fasta)
USEARCH = file(params.usearch)
FASOMERECORDS = file(params.fasomerecords)

minlen = params.minlen
mergeminlen = params.mergeminlen
mergemaxlen = params.mergemaxlen
maxee = params.maxee
qcthresh = params.qcthresh
qcpercent = params.qcpercent

chimeradb = file(params.chimeradb)
sintaxdb = file(params.sintaxdb)

FOLDER=file(params.folder)

Channel
  .fromFilePairs(FOLDER + "/*_R{1,2}_001.fastq.gz", flat: true)
  .ifEmpty { exit 1, "Could not find a matching input file" }
  .into { inputSickle }

process runSickle {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Sickle"

  input:
  set id, file(left),file(right) from inputSickle

  output:
  set id,file(left_trimmed),file(right_trimmed) into inputVsearchMerge
  set id,file(left),file(left_trimmed) into count,out

  script:
 
  left_trimmed = id + "_R1.trimmed.fastq.gz"
  right_trimmed = id + "_R2.trimmed.fastq.gz"

  """
	$SICKLE pe -g -f $left -r $right -t sanger -l $minlen -n -o $left_trimmed -p $right_trimmed -s uq.fq
  """

}

process runVsearchMerge {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/VsearchMerge"

  input:
  set id,file(left_read),file(right_read) from inputVsearchMerge

  output:
  set id, file(merged) into inputVsearchFilter
  set id,file(merged) into countMergedFq
  
  script:
  
  merged = id + ".merged.fastq"

  """
	$VSEARCH --fastq_mergepairs $left_read --reverse $right_read --fastq_minmergelen $mergeminlen --fastq_maxmergelen $mergemaxlen --fastqout $merged --threads 1
  """

}


process runVsearchFilter {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/VsearchEE"

  input:
  set id,file(merged) from inputVsearchFilter

  output: 
  set id,file(filtered) into inputFastqQualityFilter
  set id,file(filtered) into countEeFq

  script:

  filtered = id + ".merged.filtered.fastq"

  """
	$VSEARCH --fastq_filter $merged --fastq_maxee $maxee --fastqout $filtered
  """
}


process runFastqQualityFilter {


  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/FastqQualityFilter"

  input:
  set id,file(filtered) from inputFastqQualityFilter

  output:
  set id,file(quality_filtered) into inputFastqToFasta
  set id,file(quality_filtered) into countQcFq

  script:

  quality_filtered = id + ".merged.filtered.quality.fastq"

  """
	$FASTQ_QUALITY_FILTER -Q33 -q $qcthresh -p $qcpercent -i $filtered -o $quality_filtered

  """

}

process runFastqToFasta {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/FastqToFasta"

  input:
  set id,file(fastq) from inputFastqToFasta

  output:
  set id,file(fasta) into inputDerepSample
  set id,file(fasta) into inputBmSample
  
  script:

  fasta = id + ".merged.filtered.quality.fasta"

  """
	$FASTQ_TO_FASTA -i $fastq -o $fasta -Q33
  """

}


process runVsearchDerepSample {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/VsearchDerepSample"
  
  input:
  set id,file(fasta) from inputDerepSample
  
  output:
  file(derepsample) into inputMergeSamples
  
  script:
  
  derepsample = id + ".merged.filtered.quality.derep.fasta"

  """
  $VSEARCH --derep_fulllength $fasta --sizeout --output $derepsample
  """
  
  }


process runMergeSamples {

  tag "ALL"
  publishDir "${OUTDIR}/All/MergedSamples"

  input: 
  file("samples_*") from inputMergeSamples.collect()

  output: 
  file(merged) into inputVsearchDerep

  script:

  merged = "samples_merged.fasta"

  """
	cat samples_* >> $merged
  """

}


process runVsearchDerep {

  tag "ALL"
  publishDir "${OUTDIR}/All/VsearchDerep"

  input:
  file(merged) from inputVsearchDerep

  output:
  file(derep) into inputVsearchCluster

  script: 
  
  derep = "samples_derep.fasta"

  """
	$VSEARCH --derep_fulllength $merged --output $derep --minuniquesize 2 --sizeout --sizein
  """

}



process runVsearchCluster {

  tag "ALL"
  publishDir "${OUTDIR}/All/VsearchCluster"

  input:
  file(derep) from inputVsearchCluster

  output:
  file(clustered) into inputVsearchRemoveChimeraRef

  script:

  clustered = "samples_OTU.fasta"

  """
	$VSEARCH --cluster_size $derep --centroids $clustered --id 0.97 --sizein --sizeorder
  """

}



process runVsearchRemoveChimeraRef {

  tag "ALL"
  publishDir "${OUTDIR}/All/VsearchRemoveChimeraRef"

  input: 
  file(clustered) from inputVsearchRemoveChimeraRef

  output:
  file(nonchimeraref) into inputVsearchRemoveChimeraDeNovo

  script:

  nonchimeraref = "samples_OTU.nochim_ref.fasta"

  """
	$VSEARCH --uchime_denovo $clustered --strand plus --nonchimeras $nonchimeraref --sizein --sizeout
  """
}



process runVsearchRemoveChimeraDeNovo {

  tag "ALL"
  publishDir "${OUTDIR}/All/VsearchRemoveChimeraDeNovo"

  input: 
  file(nonchimeraref) from inputVsearchRemoveChimeraDeNovo

  output:
  file(filtered) into finalOTU

  script:

  filtered = "samples_OTU.nochim_ref.nochim_denovo.fasta"

  """
	$VSEARCH --uchime_denovo $nonchimeraref --strand plus --nonchimeras $filtered --sizein --sizeout --relabel OTU_
  """
}




// Merge each sample channel with the final OTU output for backmapping
inputBmSample
  .combine(finalOTU)
  .into { inputBackmap }

process runBackmap {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/Backmap"

  input:
  set id,file(fasta),file(otu) from inputBackmap
  
  output:
  set id,file(fasta),file(mapped) into inputDeChime

  script:
  
  mapped = id + ".OTU.map"

  """
	$VSEARCH --usearch_global $fasta --db $otu --strand plus --id 0.97 --uc $mapped
  """

}

process runDeChimeSample {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/DeChime"

  input:
  set id,file(fasta),file(mapped) from inputDeChime

  output:
  set id,file(nochimremapped) into inputTaxonomyAssignment
  set id,file(nochimremapped) into countDechimeFa

  script:
   
  nochimremapped = id + ".nochim.fasta"

  """
  grep "^H" $mapped | cut -f 9 | $FASOMERECORDS $fasta /dev/stdin $nochimremapped  
  """

}


process runTaxonomyAssignment {

  tag "${id}"
  publishDir "${OUTDIR}/Samples/${id}/TaxonomyAssignment"

  input:
  set id,file(nochimremapped) from inputTaxonomyAssignment

  output:
  set id,file(taxonomy) into inputFilterTax
  set id,file(nochimremapped) into inputRemoveChloroFa

  script:

  taxonomy = id + ".taxout"

  """
	$USEARCH -sintax $nochimremapped -db $sintaxdb -strand both -sintax_cutoff 0.8 -tabbedout $taxonomy
  """
}


process runFilterTax {

  tag "${id}"
  publishDir "${OUTDIR}/Final"

  input:
  set id,file(taxonomy) from inputFilterTax

  output:
  set id,file(filteredtax) into inputRemoveChloroTax

  script:
  
  filteredtax = id + ".clean.sintax"


  """
	grep -v 'Chloroplast' $taxonomy > $filteredtax
  """

}

inputRemoveChloroFa .phase(inputRemoveChloroTax) .map { p, q -> [p[0], p[1], q[1]] }.into { inputRemoveChloroOut;inputRemoveChloro }

inputRemoveChloroOut .subscribe { println it }

process runRemoveChloroplasts {

  tag "${id}"
  publishDir "${OUTDIR}/Final"

  input:
  set id,file(nochimremapped),file(filteredtax) from inputRemoveChloro

  output:
  set id,file(filteredfa) into countChloroFa

  script:
  
  filteredfa = id + ".clean.fasta"


  """
	$FASOMERECORDS $nochimremapped $filteredtax $filteredfa
  """

}

out .subscribe { println it }
count .phase(countMergedFq) .map { p, q -> [p[0], p[1], p[2], q[1]] }.into { count2;inputSampleStatsOut }
inputSampleStatsOut .subscribe { println it }
count2 .phase(countEeFq) .map { p, q -> [p[0], p[1], p[2], p[3], q[1]] }.into { count3;inputSampleStatsOut2 }
inputSampleStatsOut2 .subscribe { println it }
count3 .phase(countQcFq) .map { p, q -> [p[0], p[1], p[2], p[3], p[4], q[1]] }.into { count4;inputSampleStatsOut3 }
inputSampleStatsOut3 .subscribe { println it }
count4 .phase(countDechimeFa) .map { p, q -> [p[0], p[1], p[2], p[3], p[4], p[5], q[1]] }.into { count5;inputSampleStatsOut4 }
inputSampleStatsOut4 .subscribe { println it }
count5 .phase(countChloroFa) .map { p, q -> [p[0], p[1], p[2], p[3], p[4], p[5], p[6], q[1]] }.into { count6;inputSampleStatsOut5; inputSampleStats}
inputSampleStatsOut5 .subscribe { println it }


process runSampleStats {

  tag "${id}"
  publishDir "${OUTDIR}/Final"

  input:
  set id,file(raw),file(trimmed),file(merged),file(ee),file(qc),file(nochim),file(chloro) from inputSampleStats
 
  output:
  set id,file(stats) into outputSampleStats

  script:

  stats = id + ".clean.stats"

  shell:
  '''
	echo !{id} $(gunzip -c !{raw} | grep '^@' | wc -l) $(gunzip -c !{trimmed} | grep '^@' | wc -l) $(grep '^@' !{merged} | wc -l) $(grep '^@' !{ee} | wc -l) $(grep '^@' !{qc} | wc -l) $(grep '^>' !{nochim} | wc -l) $(grep '^>' !{chloro} | wc -l) > !{stats}
  '''

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
