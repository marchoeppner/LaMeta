#!/usr/bin/env nextflow

/*

||  ||
LaMETA
||  ||

LaMeta assembly and annotation pipeline by M. Rühlemann
Concept: M. Rühlemann
Implementation: M. Rühlemann & M. Hoeppner

*/

def helpMessage() {
  log.info"""
====================================
>> LaMETA Assembly and Annotation <<
====================================

Usage:

A typical command to run this pipeline would be:

nextflow run main.nf --groupfile my_groups.txt --reads 'data/*_R{1,2}_001.fastq.gz' 

Mandatory arguments:

--groupfile		A two-column tab-delimited file that matches fastq libraries to biological groups
--reads			A pattern to define which reads to use

Options:
--host			The location of a genome sequence in FASTA format to use as "host" reference for backmapping
--host_index		The folder containing the BBMap index matching the host fasta file (optional, will be produced otherwise)
--adapters 		A gzip compressed FASTA file with sequencing adapters (default: built-in Nextera adapter file)
--email			Provide an Email to which reports are send. 
--run_name		A name for this run.
--skip_multiqc		Do not generate a QC report
""".stripIndent()
}

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

log.info "=================================================="
log.info "LaMeta metgenomic assembly pipeline v${workflow.manifest.version}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "Authors:		M. Rühlemann & M. Höppner"
log.info "================================================="
log.info "Starting at:		$workflow.start"

/* +++++++++++++++++++++++++++++++++
Specify all input settings and files
+++++++++++++++++++++++++++++++++ */

OUTDIR=file(params.outdir)

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (params.groupfile != false ) {
	GROUP=file(params.groupfile)
	if (!GROUP.exists()) exit 1; "Specified group file could not be found (--groupfile)"
} else {
	exit 1; "Must provide a group file! See documentation for details"
}

READMINLEN = params.readminlen

// Do sanity checks on input files

// Allow providing a custom adapter file, or fall back to the built-in nextera adapters
if(params.adapters) {
	ADAPTERS = file(params.adapters)
} else {
	ADAPTERS = file(workflow.projectDir + "/assets/adapters/nextera.fa.gz")
}
if (!ADAPTERS.exists()) exit 1; "Unable to find specified adapter sequence file (--adapters)"

// Specify a host genome to map against (assume a mapping index exists)
if (params.host) {
	HOST = file(params.host)
	if ( !HOST.exists() ) exit 1; "Unable to find host reference sequence (--host)"
}	

if (params.host_index) {
	HOST_INDEX = file(params.host_index)
	if (!HOST_INDEX.exists() ) exit 1; "Could not located the BBMap index (--host_index)"
	HOST_INDEX_REF = file(params.host_index + "/ref")
	if (!HOST_INDEX_REF.exists() ) exit 1; "The specified host index directory does not contain the BBMap ref folder (--host_index)"
} 

if (params.host !=false && params.host_index != false) {
	println "Specified both a host fasta file and a pre-compiled index. Will ignore the index!"
} else if (params.host == false && params.host_index == false) {
	exit 1; "Provided neither host genome fasta file (--host) nor a BBMap index location (--host_index); cannot proceed without one of the two."
}

if ( params.gtdb != false ) {
	GTDB_FILE = file(params.gtdb)
	if (!GTDB_FILE.exists()) exit 1; "Could not find the specified GTDB-TK Database (--gtdb)"
} else {
	exit 1, "Must prodvide the path to the GTDB-TK Database (--gtdb)"
}

if (params.checkm_db != false ) {
	if (!file(params.checkm_db).exists() ) exit 1; "Could not find the CheckM database (--checkm_db)"
}
	
SPADES_kmers = params.spades_kmers

/* 
Gather information for summary
*/

def summary = [:]

summary['runName'] = run_name
summary['Reads'] = params.reads
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['HostGenome'] = params.host
summary['Adapters'] = params.adapters
summary['SessionID'] = workflow.sessionId
summary['LaMeta Version'] = workflow.manifest.version
summary['StartedAt'] = workflow.start

/*
Channel is created from the file pattern given by --reads.
*/

// Pass the GTDBTK Database as a channel so we can import it as an environment variable into the process
Channel
   .from(params.gtdb + "/")
   .into { inputGTDB; inputGTDBMarkers }

// CheckM database - use existing or retrieve on the fly
if (singurlaty.enabled == false && docker.enabled == false) }

	if (params.checkm_db != false ) {
		CHECKM_DB = Channel.fromPath(params.checkm_db)
	} else {

		CHECKM_DB_URL = Channel.from(params.checkm_db_url)

		// make this a local process for cases where executing nodes may not have access to the web
		process getCheckMDB {

			executor = "local"

			tag "ALL"
			publishDir "${OUTDIR}/refs", mode: 'copy'

			input:
			val(url) from CHECKM_DB_URL

			output:
			file("${checkm_db_path}") into CHECKM_DB

			script:	
			checkm_db_path = "checkm_data"

			"""
				mkdir -p checkm_data 
				cd checkm_data
				wget $url
				tar xzf checkm_data_2015_01_16.tar.gz
				rm *.tar.gz
			"""
		}
	}
} else  {
	// this is where the data lives inside the containter
	CHECKM_DB = Channel.fromPath("/db/checkm_data")
}

// set checkm db location
process runSetCheckmRoot {

	tag "ALL"

	input:
	file(checkm_db_path) from CHECKM_DB

	output:
	file(checkm_db_path) into (checkmdb_out, checkmdb_out_refine)

	script:
	

	"""
		printf "checkm_data\ncheckm_data\n" | checkm data setRoot
	"""

}

// Read fastq files as paired - complain if nothing is found
Channel
  .fromFilePairs(params.reads, flat: true)
  .ifEmpty { exit 1, "Could not find reads matching the specified input format/location" }
  .into { inputQC; inputParseGroupPre }

inputParseGroupPre.map{id, f1, f2 -> id}.set{ inputParseGroup }

// Parse the group file and get the relevant metadata
process parseGroup {

	tag "${id}"

	input:
	val id from inputParseGroup

	output:
	set id, stdout into outputParseGroup1, outputParseGroup2, outputParseGroup3

	script:
	"""
	grep $id $GROUP | awk '{printf \$2}'
	"""
}

/*
Mapping against PhiX and Host genome (default:human). Mapped reads/read-pairs (also discordantly)
are discarded.
*/

if (params.host != false ) {

	Channel
	  .fromPath(params.host)
	  .set { inputBBMapIndex }

	// Create a BBMAP compatible index structure and pass to QC stage
	process runBuildIndex {

		tag "All"
		 publishDir "${OUTDIR}/Host/", mode: 'copy'

		input:
		file(genome_fa) from inputBBMapIndex

		output:	
		set file(index_dir) into BBMapIndex

		script:
		index_dir = "ref"

		"""
		bbmap.sh -Xmx${task.memory.toGiga()}g t=${task.cpus} ref=$genome_fa
		"""
	}
} else if (params.host_index != false ) {
        BBMapIndex = Channel.from(HOST_INDEX)
}

inputQCIndex = inputQC.combine(BBMapIndex)

process runQC {

	// scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Decon", mode: 'copy'

	input:
	set id, file(left),file(right),file(genome_index) from inputQCIndex

	output:
	set id,file(left_clean),file(right_clean),file(unpaired_clean) into inputSpades, inputSpadesBackmap, inputCoAssemblyPre
	file(bbduk_adapter_stats) into outputQCstats

	script:

	left_trimmed = "tmp_" + id + "_R1.trimmed.fastq"
	right_trimmed = "tmp_" + id + "_R2.trimmed.fastq"
	unpaired_trimmed = "tmp_" + id + "_RU.trimmed.fastq"

	left_nophix = "tmp_" + id + "_R1.nophix.fastq"
	right_nophix = "tmp_" + id + "_R2.nophix.fastq"
	unpaired_nophix = "tmp_" + id + "_RU.nophix.fastq"

	left_decon = "tmp_" + id + "_R1.decon.fastq"
	right_decon = "tmp_" + id + "_R2.decon.fastq"
	unpaired_decon = "tmp_" + id + "_RU.decon.fastq"

	merged = "tmp_" + id + "_RU.merged.fastq"
	left_clean = id + "_R1.clean.fastq.gz"
	right_clean = id + "_R2.clean.fastq.gz"
	unpaired_clean = id + "_RU.clean.fastq.gz"

	bbduk_adapter_stats = id + ".bbduk.adapters.stats.txt"
	bbduk_artifact_stats = id + ".bbduk.artifacts.stats.txt"
	
	finalstats = id +".stats.txt"

	"""
	reformat.sh threads=${task.cpus} in=${left} in2=${right} 2>&1 >/dev/null | awk '{print "RAW "\$0}' | tee stats.txt 
    	bbduk.sh stats=$bbduk_adapter_stats threads=${task.cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=${unpaired_trimmed} ref=${ADAPTERS} ktrim=r k=23 mink=11 hdist=1 minlength=${READMINLEN} tpe tbo
    	bbduk.sh stats=$bbduk_artifact_stats threads=${task.cpus} in=${left_trimmed} in2=${right_trimmed} k=31 ref=artifacts,phix ordered cardinality out1=${left_nophix} out2=${right_nophix} minlength=${READMINLEN}
	bbduk.sh threads=${task.cpus} in=${unpaired_trimmed}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_nophix} minlength=${READMINLEN}
	bbwrap.sh -Xmx23g threads=${task.cpus} minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=20 minlength=${READMINLEN} in=${left_nophix},${unpaired_nophix} in2=${right_nophix},NULL path=${genome_index} outu1=${left_decon} outu2=${right_decon} outu=${unpaired_decon} 2>&1 >/dev/null | awk '{print "HOST "\$0}' | tee -a stats.txt 
	bbmerge.sh threads=${task.cpus} in1=${left_decon} in2=${right_decon} out=${merged} outu1=${left_clean} outu2=${right_clean} mininsert=${READMINLEN} 2>&1 >/dev/null | awk '{print "MERGED "\$0}' | tee -a stats.txt
	cat ${merged} ${unpaired_decon} | gzip -c > ${unpaired_clean}
	reformat.sh threads=${task.cpus} in=${unpaired_clean} 2>&1 >/dev/null | awk '{print "UNPAIRED "\$0}' | tee -a stats.txt
	reformat.sh threads=${task.cpus} in1=${left_clean} in2=${right_clean}  2>&1 >/dev/null | awk '{print "PAIRED "\$0}' | tee -a stats.txt
	rm tmp*

	grep "RAW" stats.txt  | grep 'Input:' | awk '{print "READS RAW "\$3/2}' | tee $finalstats
 	grep "HOST" stats.txt | grep "Reads Used:"  | awk '{printf \$4" "}' | awk '{print "READS BIO "\$1/2 + \$2}' | tee -a $finalstats
	egrep "^UNPAIRED" stats.txt  | grep 'Input:' | awk '{print \$3}' | awk '{print "READS CLEAN_UNPAIRED "\$1}' | tee -a $finalstats
	egrep "^PAIRED" stats.txt  | grep 'Input:' | awk '{print \$3}' | awk '{print "READS CLEAN_PAIRED "\$1}' | tee -a $finalstats
	grep "RAW" stats.txt  | grep 'Input:' | awk '{print "BASES RAW "\$5}' | tee -a $finalstats
	egrep "^UNPAIRED" stats.txt  | grep 'Input:' | awk '{print "BASES CLEAN_UNPAIRED "\$5}' | tee -a $finalstats
	egrep "^PAIRED" stats.txt  | grep 'Input:' | awk '{print "BASES CLEAN_PAIRED "\$5}' | tee -a $finalstats

    	"""
}

/*
Co assembly within the groups given in groupfile.
*/

outputParseGroup1.join(inputCoAssemblyPre).map{id, group, left_clean, right_clean, unpaired_clean -> [group, id, left_clean, right_clean, unpaired_clean]}.set{inputCoAssembly}

inputCoAssembly.groupTuple().into{ inputCoAssemblyByGroup; inputBackmapMegahit }

process runCoAssembly {

	// scratch true

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}", mode: 'copy'

	input:
	set group, id, file(left_clean), file(right_clean), file(unpaired_clean) from inputCoAssemblyByGroup

	output:
	set group, file(outcontigs), file(megahitlog) into inputContigsBackmapMegahit, inputContigsMegahitMaxbin, inputContigsMegahitMetabat, inputMegahitMarkergenes, inputMegahitRefine

	script:
	outcontigs = group + ".final_contigs.fasta"
	megahitlog = group + ".megahit.log"

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

  	megahit \$(cat tmp | tr -d '\n') --num-cpu-threads ${task.cpus} --presets meta-large -o megahit_out --mem-flag 2 --verbose
  	cat megahit_out/final.contigs.fa | cut -d ' ' -f 1 | awk -v group=$group '/^>/{print ">Megahit_"group"_contig_"++i; next}{print}' > $outcontigs
  	mv megahit_out/log $megahitlog
  	rm -r megahit_out
  	"""
}

/*
Single-samples metagenome assembly with spades
*/

process runSpades {

	scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

	input:
	set id, file(left_clean), file(right_clean), file(unpaired_clean) from inputSpades

	output:
	set id, file(outcontigs) into inputSpadesBackmapContigs, inputSpadesMaxbin, inputSpadesMarkergenes, inputSpadesRefine

	script:
	outcontigs = id + ".spades_contigs.fasta"

  	"""
  	spades.py --meta --pe1-1 $left_clean --pe1-2 $right_clean --pe1-s $unpaired_clean -k $SPADES_kmers -o spades_out -t ${task.cpus}
  	awk -v id=$id '/^>/{print ">Spades_"id"_contig_"++i; next}{print}' spades_out/scaffolds.fasta > $outcontigs
	rm -r spades_out
  	"""
}


/*
Backmapping to spades assembly and contig abundance estimation
*/

inputSpadesBackmap.join(inputSpadesBackmapContigs).set { inputSpadesBackmapWithContigs}

process runSpadesBackmap {

	scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Spades", mode: 'copy'

	input:
	set id, file(left_clean), file(right_clean), file(unpaired_clean), file(spadescontigs) from inputSpadesBackmapWithContigs

	output:
	set id, file(outdepth) into outputSpadesBackmap

	script:
	def avail_mem = "${(task.memory.toBytes()-1000000000)/task.cpus}"
	outdepth = id + ".depth.txt"
	"""
  	bbwrap.sh -Xmx60g in=$left_clean,$unpaired_clean in2=$right_clean,NULL ref=$spadescontigs t=${task.cpus} out=tmp.sam kfilter=22 subfilter=15 maxindel=80
  	samtools view -u tmp.sam | samtools sort -m ${avail_mem} -@ ${task.cpus} -o tmp_final.bam
  	jgi_summarize_bam_contig_depths --outputDepth $outdepth tmp_final.bam
  	rm tmp*
  	rm -r ref
	"""
}

/*
Single-sample binning with Maxbin2
*/
inputSpadesMaxbin.join(outputSpadesBackmap).into {inputMetabat; inputMaxbin; inputMaxbin40}

process runMaxbin {

	scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Maxbin", mode: 'copy'

	input:
	set id, file(spadescontigs), file(depthfile) from inputMaxbin

	output:
	set id, file(binfolder) into outputMaxbinSamples

	script:
	binfolder = id + "_maxbin_bins"
 
	"""

	tail -n+2 $depthfile | cut -f 1,3 > maxbin.cov
 	mkdir -p $binfolder
	(
	set -Ee
	function _catch {
		touch summary.txt
		echo "exception caught"
		exit 0
	}
	trap _catch ERR
  	run_MaxBin.pl -contig $spadescontigs -abund maxbin.cov -out $binfolder/${id}.bin -thread ${task.cpus}
	)

 	"""
}


process runMaxbin40 {

	scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Maxbin40", mode: 'copy'

	input:
	set id, file(spadescontigs), file(depthfile) from inputMaxbin40

	output:
	set id, file(binfolder) into outputMaxbin40Samples

	script:
	binfolder = id + "_maxbin40_bins"
  	"""

	tail -n+2 $depthfile | cut -f 1,3 > maxbin.cov
  	mkdir -p $binfolder
  	mkdir workfolder
  	mkdir tmp_workfolder
	(
	set -Ee
	function _catch {
		touch summary.txt
		echo "exception caught"
		exit 0
	}
	trap _catch ERR
	run_MaxBin.pl -contig $spadescontigs -abund maxbin.cov -out $binfolder/${id}.bin40 -thread ${task.cpus} -markerset 40
	)

  	"""
}

process runMetabat {

	scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/Metabat", mode: 'copy'

	input:
	set id, file(spadescontigs), file(depthfile) from inputMetabat

	output:
	set id, file(binfolder) into outputMetabatSamples

	script:
	binfolder = id + "_metabat_bins"

	"""
	mkdir $binfolder
	metabat2 -i $spadescontigs -a $depthfile -o $binfolder/${id}.metabat.bin -t ${task.cpus}
	"""
}

process runSpadesMarkergenes {

	// scratch true

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/SpadesMarkergenes", mode: 'copy'
	
	input:
	set id, file(spadescontigs) from inputSpadesMarkergenes
	env GTDBTK_DB_PATH from inputGTDB

	output:
	set id, file(markergenes) into outputSpadesMarkergenes

	script:
	markergenes = id + "_markergenes.txt"

	"""
	mkdir tmp
	cp $spadescontigs tmp
	gtdbtk identify --genome_dir tmp -x fasta --cpus ${task.cpus} --out_dir markers
	cat markers/marker_genes/*/*tophit.tsv | grep -v hits | tr "," "\t" | cut -d ';' -f 1 > $markergenes
	"""
}

inputSpadesRefine.join(outputSpadesMarkergenes).join(outputMaxbinSamples).join(outputMaxbin40Samples).join(outputMetabatSamples).into{ SamplesAllbins; testtest}

process runSpadesRefine {

	tag "${id}"
	publishDir "${OUTDIR}/Samples/${id}/ContigsRefined", mode: 'copy'

	input:
	set id, file(spadescontigs), file(markergenes), file(binmaxbin), file(binmaxbin40), file(binmetabat) from SamplesAllbins
	file(db_path) from checkmdb_out

	output:
	set id, file(refinedcontigsout) into SampleRefinedContigs

	script:
	refinedcontigsout = id + "_refined"

	"""
	R CMD javareconf > /dev/null 2>&1 || true
	mkdir $refinedcontigsout
	grep '>' ${binmaxbin}/*fasta | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 > btc.txt
	grep '>' ${binmaxbin40}/*fasta | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 >> btc.txt
	grep '>' ${binmetabat}/*fa | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 >> btc.txt
	DAS_Tool_ripoff.Rscript ${id} btc.txt ${markergenes} 
  	sort_into_bins.py ${id}.refined.contig_to_bin.out ${spadescontigs}
  	mkdir -p bins
  	mkdir -p $refinedcontigsout/bins
	mv ${id}_cleanbin_*.fasta bins
	chd=\$(readlink -f $db_path)
	printf "\$chd\\n\$chd\\n" | checkm data setRoot
	checkm lineage_wf -t ${task.cpus} -x fasta --nt --tab_table -f ${id}.checkm.out bins checkm_out
	head -n 1 ${id}.checkm.out > $refinedcontigsout/${id}.checkm.out
	for good in \$(awk -F '\t' '{if(\$12 > 50 && \$1!="Bin Id") print \$1}' ${id}.checkm.out); do mv bins/\$good.fasta $refinedcontigsout/bins; grep -w \$good ${id}.checkm.out >> ${refinedcontigsout}/${id}.checkm.out; done
	 mv ${id}.refined* $refinedcontigsout
  	"""
}

/*
Backmapping to Megahit groupwise co-assembly
*/

inputBackmapMegahit.transpose().combine(inputContigsBackmapMegahit, by: 0) .set { inputBackmapCoassemblyT }

process runCoassemblyBackmap {

	scratch true

	tag "${group}-${id}"
	publishDir "${OUTDIR}/CoAssembly/${group}/Backmap", mode: 'copy'

	input:
	set group, id, file(left_clean), file(right_clean), file(unpaired_clean), file(megahitcontigs), file(megahitlog) from inputBackmapCoassemblyT

	output:
	set group, file(bamout) into outMegahitBackmap

	script:
	bamout = id + ".megahit.final.bam"
	mem_per_thread = (task.memory.toGiga() / task.cpus)

	"""
	bbwrap.sh -Xmx60g in=$left_clean,$unpaired_clean in2=$right_clean,NULL ref=$megahitcontigs t=${task.cpus} out=tmp_sam.gz kfilter=22 subfilter=15 maxindel=80
	samtools view -u tmp_sam.gz | samtools sort -m ${mem_per_thread}G -@ ${task.cpus} -o $bamout
	rm tmp*
	rm -r ref
	"""
}

/*
Contig abundance estimation for co-assemblies
*/

outMegahitBackmap.groupTuple().set { inputCollapseBams }

process runCollapseBams {

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/Backmap", mode: 'copy'

	input:
	set group, file(bams) from inputCollapseBams

	output:
	set group, file(depthfile) into coassemblyDepth
	set group, file(abufolder) into coassemblyAbufolder

	script:
	depthfile = group + "_depth.txt"
	abufolder = group + "_abufiles"

	"""
	jgi_summarize_bam_contig_depths --outputDepth $depthfile $bams
	ncol=\$(head -n 1 $depthfile | awk '{print NF}')
	mkdir -p $abufolder
	for i in \$(seq 4 2 \$ncol); do
		name=\$(head -n 1 $depthfile | cut -f \$i | cut -d "." -f 1)
		cut -f  1,\$i $depthfile | tail -n+2 > $abufolder/\${name}.out
	done
	"""
}

/*
Co-assembly binning with Maxbin2
*/

coassemblyAbufolder.join(inputContigsMegahitMaxbin).into{ inputMegahitMaxbin; inputMegahitMaxbin40 }

process runMegahitMaxbin {

	scratch true

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/Maxbin", mode: 'copy'

	input:
	set group, file(inputfolder), file(megahitcontigs), file(megahitlog)  from inputMegahitMaxbin

	output:
	set group, file(binfolder) into outputMegahitMaxbin

	script:
	binfolder = group + "_maxbin_bins"

	"""
	ls ${inputfolder}/*.out > abufiles.txt
	mkdir $binfolder
	run_MaxBin.pl -contig $megahitcontigs -abund_list abufiles.txt -out $binfolder/${group}.maxbin.bin -thread ${task.cpus}

  	"""
}

process runMegahitMaxbin40 {

	scratch true  

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/Maxbin40", mode: 'copy'

	input:
	set group, file(inputfolder), file(megahitcontigs), file(megahitlog)  from inputMegahitMaxbin40

	output:
	set group, file(binfolder) into outputMegahitMaxbin40

	script:
	binfolder = group + "_maxbin40_bins"

	"""
	ls ${inputfolder}/*.out > abufiles.txt
	mkdir $binfolder
	run_MaxBin.pl -contig $megahitcontigs -abund_list abufiles.txt -out $binfolder/${group}.maxbin.bin40 -thread ${task.cpus} -markerset 40
	"""
}


/*
Co-assembly binning with Metabat
*/

coassemblyDepth.join(inputContigsMegahitMetabat).set{ inputMegahitMetabat }

process runMegahitMetabat {

	scratch true

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/Metabat", mode: 'copy'

	input:
	set group, file(inputdepth), file(megahitcontigs), file(megahitlog) from inputMegahitMetabat

	output:
	set group, file(binfolder) into outputMegahitMetabat

	script:
	binfolder = group + "_metabat_bins"

	"""
	mkdir $binfolder
	metabat2 -i $megahitcontigs -a $inputdepth -o $binfolder/${group}.metabat.bin -t ${task.cpus}
	"""
}

process runMegahitMarkergenes {

	scratch true

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/MegahitMarkergenes", mode: 'copy'

	input:
	set group, file(megahitcontigs), file(megahitlog) from inputMegahitMarkergenes
        env GTDBTK_DB_PATH from inputGTDBMarkers

	output:
	set group, file(markergenes) into outputMegahitMarkergenes

	script:
	markergenes = group + "_markergenes.txt"

	"""
	mkdir tmp
	cp $megahitcontigs tmp
	gtdbtk identify --genome_dir tmp -x fasta --cpus ${task.cpus} --out_dir markers
	cat markers/marker_genes/*/*tophit.tsv | grep -v hits | tr "," "\t" | cut -d ';' -f 1 > $markergenes
	rm -r tmp
	"""
}

inputMegahitRefine.join(outputMegahitMarkergenes).join(outputMegahitMaxbin).join(outputMegahitMaxbin40).join(outputMegahitMetabat).into{ MegahitAllbins; testtest2}

process runMegahitRefine {

	tag "${group}"
	publishDir "${OUTDIR}/CoAssembly/${group}/ContigsRefined", mode: 'copy'

	input:
	set group, file(megahitcontigs), file(megahitlog), file(markergenes), file(binmaxbin), file(binmaxbin40), file(binmetabat) from MegahitAllbins
        file(db_path) from checkmdb_out_refine


	output:
	set group, file(refinedcontigsout) into MegahitRefinedContigs

	script:
	refinedcontigsout = group + "_refined"

	"""
	mkdir $refinedcontigsout
	grep '>' ${binmaxbin}/*fasta | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 > btc.txt
	grep '>' ${binmaxbin40}/*fasta | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 >> btc.txt
	grep '>' ${binmetabat}/*fa | tr ':' ' ' | tr -d '>'  | cut -d '/' -f 2 >> btc.txt
	DAS_Tool_ripoff.Rscript ${group} btc.txt ${markergenes}
	sort_into_bins.py ${group}.refined.contig_to_bin.out ${megahitcontigs} 
	mkdir bins
	mkdir -p $refinedcontigsout/bins
	mv ${group}_cleanbin_*.fasta bins
	chd=\$(readlink -f ${db_path})
        printf "\$chd\\n\$chd\\n" | checkm data setRoot
	checkm lineage_wf -t ${task.cpus} -x fasta --nt --tab_table -f ${group}.checkm.out bins checkm_out
	head -n 1 ${group}.checkm.out > $refinedcontigsout/${group}.checkm.out
	for good in \$(awk -F '\t' '{if(\$12 > 50 && \$1!="Bin Id") print \$1}' ${group}.checkm.out); do mv bins/\$good.fasta $refinedcontigsout/bins; grep -w \$good ${group}.checkm.out >> $refinedcontigsout/${group}.checkm.out; done 
	mv ${group}.refined* $refinedcontigsout
	"""
}

process runMultiQCLibrary {

	tag "ALL"
	publishDir "${OUTDIR}/MultiQC/Library", mode: 'copy'

	when:
	params.skip_multiqc != true
	
	input:
	file(reports) from outputQCstats.collect()

	output:
	file(multiqc_library) into multiqc_report

	script:
	multiqc_library = "multiqc_library.html"

	"""
		multiqc -n multiqc_library.html *.txt
	"""

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  summary["finishedAt"] = workflow.complete
  summary["duration"] = workflow.duration

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "LaMeta analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[LaMeta] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[LaMeta] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}
