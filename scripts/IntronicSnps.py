import sys
import collections
import gzip
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --sync input.sync --gff transcriptome.gff --target-length 60 > intron60.sync  "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script identifies SNPs (--sync) located in introns (-gff) of a given length (--target-length)
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--gff", dest="gff", help="A GFF file for the transcriptome")
parser.add_option("--sync", dest="input", help="A sync file")
parser.add_option("--target-length", dest="intron", help="The target intron length")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
	''' import data either from a gzipped or or uncrompessed file or from STDIN'''
	import gzip         
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y


gff=load_data(options.gff)
sync=load_data(options.sync)
intron_length_threshold=int(options.intron)

intronhash=collections.defaultdict(lambda: collections.defaultdict(list))

for l in gff:
    if l.startswith("#"):
        continue
    if len(l.split())<8:
        continue
    chrom,source,typ,start,end=l.split()[:5]
    if source!="FlyBase":
        continue
    if typ!="intron":
        continue
    if abs(int(end)-int(start))>intron_length_threshold:
        continue
    if int(end)>int(start):
        for i in range(int(start),int(end)+1):
            intronhash[chrom][i]
    else:
        for i in range(int(end),int(start)+1):
            intronhash[chrom][i]
  
for l in sync:
    if len(l.split("\t"))<2 or l.startswith("#"):
        continue
    chrom,pos=l.split()[:2]
    if int(pos) in intronhash[chrom]:
        print l.rstrip()


