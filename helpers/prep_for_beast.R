library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

#candidate="g__Collinsella"
candidate=rev(strsplit(getwd(),split="/")[[1]])[1]


checkm<-read.table(paste0(candidate,".checkm.out"),sep="\t",row.names=1)
mummer_files<-list.files("MUMMER")
genomes= unique(sapply(mummer_files, function(x) strsplit(x,"[.]")[[1]][1]))
checkm=checkm[genomes,]
checkm=checkm[grepl("k__Bacteria (UID203)",checkm$V2)==F,]
checkm.good=checkm[checkm$V12>=80 & checkm$V13<10 & checkm$V14<10,]
genomes<-rownames(checkm)

allranges<-vector("list",length(genomes))
names(allranges)<-genomes

for(genome in genomes){
reader<-read.table(paste0("MUMMER/",genome,".filtered.coords"),head=F,stringsAsFactors=F)
raw<-GRanges(seqnames=Rle(reader$V12), ranges=IRanges(reader$V1,reader$V2), strand="+")
masked<-GRanges(slice(coverage(raw),upper=1))
allranges[[genome]]<-masked
}

if(exists("cranges")){rm("cranges")}

first=0
for(genome in rownames(checkm.good)){
if(first==0){cranges=allranges[[genome]]}else{cranges=c(cranges,allranges[[genome]])}
first=1
}

cranges.cov<-coverage(cranges)

sl<-slice(cranges.cov,lower=nrow(checkm.good))

covered<-GRanges(sl)

allsnps<-vector("list",length(genomes))
names(allsnps)<-genomes

for(genome in genomes){
snpfile=paste0("MUMMER/",genome,".filtered.snps")
nrow=system(paste0("wc -l ",snpfile, " | cut -f 1 -d ' ' "),intern=T)
if(nrow=="0"){ref=genome; next}
reader<-read.table(snpfile,head=F,stringsAsFactors=F)
reader<-reader[reader$V3!="." & reader$V2!=".",]
snps<-GRanges(seqnames=Rle(reader$V11), ranges=IRanges(reader$V1,reader$V1), strand="+", ref=reader$V2, genotype=reader$V3)

snpsinrange=data.frame(genome=genome,subsetByOverlaps(snps,covered))
allsnps[[genome]]<-snpsinrange
}

cc<-do.call("rbind", allsnps)
cc$snpid=paste(cc$seqnames, cc$start,sep="_")
head(cc)
cc2<-cc
reference<-unique(cc[,c("snpid","seqnames","start","end","ref")])
rownames(reference)<-reference$snpid
allsnps.ranges<-GRanges(seqnames=Rle(reference$seqnames), ranges=IRanges(reference$start,reference$end), strand="+", ref=reference$ref,snpid=rownames(reference))

reffasta<-getSeq(FaFile(file=paste0(candidate,".ref.fasta")),covered)
constants<-table(strsplit(as.character(unlist(reffasta)),split="")[[1]])[c("A","C","G","T")] - table(reference$ref)

write.table(constants,paste0(candidate,".constants.txt"),quote=F,row.names=F,col.names=F)

allsnps.final<-vector("list",length(genomes))
names(allsnps.final)<-genomes

for(genome in genomes){
print(genome)
this<-data.frame(subsetByOverlaps(allsnps.ranges,allranges[[genome]]))
thissnps<-cc[cc$genome==genome,]
thisref<-this[!this$snpid %in% thissnps$snpid,]
if(nrow(thisref)>0){
thisref$genotype=thisref$ref
thisref$genome=genome
thisref<-thisref[,colnames(thissnps)]}else{thisref=thissnps[0,]}
thisna<-reference[!rownames(reference) %in% this$snpid,]
if(nrow(thisna)>0){thisna$genotype="N"
thisna$genome=genome
thisna$strand="+"
thisna$width=1
thisna<-thisna[,colnames(thissnps)]}else{thisna<-thisref[0,]}
allsnps.final[[genome]]<-rbind(thissnps,thisref,thisna)
}

cc.final<-do.call("rbind", allsnps.final)

aa<-data.frame(reshape2::dcast(genome ~ snpid, data=cc.final, value.var="genotype"),row.names=1)
aa<-aa[,sort(colnames(aa))]

library(seqRFLP)
groups<-read.table("../../../../groupfile.txt",row.names=1,stringsAsFactors=F)

df<-data.frame(names=paste0(groups[gsub("_cleanbin_[0-9]+","",rownames(aa)),],"_",rownames(aa)),seqs=apply(aa,1,paste, collapse=""),stringsAsFactors=F)


df.fasta = dataframe2fas(df, file=paste0(candidate,"_SNPs_for_BEAST.fasta"))

