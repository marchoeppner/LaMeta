import sys
import csv
import pickle
from numpy import median, mean

allkeys=[]
print("Reading sequence keys...")
with open(sys.argv[1],"r") as file:
        for line in file:
                a,b=line.strip().split(" ")
                allkeys.append([a,[float(0)]*int(b)])

mapdict=dict(allkeys)

print("Calculating coverage...")
with open(sys.argv[2],"r") as depthfile:
#with open("H07602-L1_S6_L001.dmnd.out","r") as depthfile:
    reader = csv.reader(depthfile, delimiter='\t')
    for s,row in enumerate(reader):
        if s % 100000 == 0:
            print(s)
        if s == 0:
            thisseq = row[0]
            thisseqcount=0
            thisseqlist=[]
            totseqcount=0
        if thisseq != row[0]:
            for entry in thisseqlist:
                    gene=entry[0]
                    for pos in range(entry[1],entry[2]):
                        mapdict[gene][pos] += float(1/thisseqcount)
            thisseqcount=0
            del thisseqlist[:]
            thisseq = row[0]
            totseqcount=totseqcount+1
        else:
            thisseqlist.append([row[1],int(row[6])-1,int(row[7])-1])
            thisseqcount=thisseqcount+1

pmsf=float(float(totseqcount)/1000000)

print("Writing output...")

outfile= open(sys.argv[3]+".feature_abundance.tsv", 'w')

outfile.write("".join(["\t".join(["Feature","Length","Covered","CoveredGreaterOne","MedianCoverage","MedianCoveragePerMillion","MeanCoverage","MeanCoveragePerMillion"]),"\n"]))
for feature in mapdict.keys():
    cov=mapdict[feature]
    totresid=len(cov)
    residcov0=sum(x>0 for x in cov)
    residcov1=sum(x>1 for x in cov)
    medcov=median(cov)
    meancov=mean(cov)
    cpmmr=float(medcov/pmsf)
    cpmmrmean=float(meancov/pmsf)
    newline= "\t".join([feature,str(totresid), str(residcov0), str(residcov1),str(medcov),str(cpmmr),str(meancov),str(cpmmrmean)])
    outfile.write("".join([newline,"\n"]))

outfile.close()

