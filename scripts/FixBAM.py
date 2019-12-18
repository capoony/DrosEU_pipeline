import pysam 
import sys
import os 
from collections import defaultdict
from optparse import OptionParser, OptionGroup

#########################################################   HELP   #########################################################################
usage="python %prog --input contaminated.bam --contaminated contaminated.bam --prefix sim_ --detect assortative_mapping.sam --output clean"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this. This script reads the Read IDs from the SAM file (--detect) and sort them according to the presence or absence of the prefix for the 'contaminant' chromosomes. According to these two ID lists, reads in the contaminated BAM file (--contaminated) will be split in multiple different BAM files (--output).

""") 

parser.add_option("--contaminated", dest="con", help="A BAM file")
parser.add_option("--detect", dest="det", help="A sam file ")
parser.add_option("--prefix", dest="prefix", help="The prefix added to the chromosome names in the 'contaminant' reference genome")
parser.add_option("--output", dest="out", help="the melanogaster specific fwd FASTQ")

parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the SAM file                                                                            
mel,sim=defaultdict(str),defaultdict(str)
samfile=pysam.Samfile(options.det,"r")
prefix=options.prefix
c=1
for l in samfile.fetch(until_eof=True):
	if c%1000000==0:
		print c,"positions processed"
	if l.tid<0:
		continue
	elif samfile.getrname(l.tid).startswith(prefix):
		sim[l.qname]
	else:
		mel[l.qname]
	c+=1

print "reading SAM done"
	
## index BAM file if necessary
if not os.path.exists(options.con+".bai"):
	print "indexing "+options.con
	os.system("samtools index "+options.con)

samfile=pysam.Samfile(options.con,"rb")
melout=pysam.Samfile(options.out+"-mel.bam","wb",template=samfile)
simout=pysam.Samfile(options.out+"-sim.bam","wb",template=samfile)
missedout=pysam.Samfile(options.out+"-missed.bam","wb",template=samfile)
## split BAM file

c=1
for l in samfile.fetch(until_eof=True):
	if c%1000000==0:
		print c,"positions processed"
	if l.qname.split(" ")[0] in mel:
		melout.write(l)
	elif l.qname.split(" ")[0]  in sim:
		simout.write(l)
	else:
		missedout.write(l)
	c+=1
	
melout.close()
simout.close()
missedout.close()
samfile.close()
