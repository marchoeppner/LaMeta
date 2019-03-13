from Bio import SeqIO
import sys

args = sys.argv

id_file = args[1]
input_file = args[2]

sortdict={}
print("Reading contig to bin mapping...")
for i,line in enumerate(open(id_file)):
    if(i>0):
        thisline=line.rstrip("\n").split("\t")
        sortdict[thisline[1]]=thisline[0]

total=len(sortdict)
percent=total/100
counter=0
print(total)
write=False

print("Sorting contigs into bins...")
for r in SeqIO.parse(input_file, "fasta"):
    if r.id in sortdict.keys():
        with open(sortdict[r.id], "a") as outfile:
            towrite=">"+str(r.id)+"\n"+str(r.seq)+"\n"
            outfile.write(towrite)
        counter += 1
	write=True
    if (counter % percent == 0 and counter != 0 and write==True):
	print("Wrote %i of %i records into bins" % (counter, total))
	write=False
    if counter == total:
	break
