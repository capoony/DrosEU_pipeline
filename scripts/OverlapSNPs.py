import sys
from optparse import OptionParser, OptionGroup
import gzip

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --source input.vcf --target genes.sync > input.sync"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Print rows from target file (--target) that are also present in the source file (--source) based on the info in the first two columns (Chrom and Position )""") 
#########################################################   CODE   #########################################################################

parser.add_option("--source", dest="source", help="The source file with Chromosome and position in the first and second column")
parser.add_option("--target", dest="target", help="The target file with Chromosome and position in the first and second column")

parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
	import gzip			
	if x=="-":
		y=sys.stdin
	elif x.endswith(".gz"):
		y=gzip.open(x,"r")
	else: 
		y=open(x,"r")
	return y

genehash={}
	
pop1=load_data(options.source)
pop2=load_data(options.target)

for l in pop1:
	if len(l.split())>1:
		genehash[l.split()[0]+"_"+l.split()[1]]=l.rstrip()

for l in pop2:
	if len(l.split())>1:
		if l.split()[0]+"_"+l.split()[1] in genehash:
			print l.rstrip()
