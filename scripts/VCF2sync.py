import sys
import collections
import random
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --vcf input.vcf > input.sync  "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script converts a file in VCF format into the SYNC file format.
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--vcf", dest="vcf", help="A vcf file")

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


def sync_format(w):
	al=["A","T","C","G","N","D"]
	sy=[]
	for nuc in al:
		if nuc in w:
			sy.append(str(w[nuc]))
		else:
			sy.append("0")
	return ":".join(sy)

for l in load_data(options.vcf):
	if len(l.split())<9:
		continue
	if l.startswith("#"):
		continue
	chr,pos,ID,ref,alt,qual,Filter,Info,Format=l.split()[:9]
	formID=Format.split(":")
	if "DELETION" in Info or "INSERTION" in Info:
		continue
	ALT=alt.split(",")
	pops=l.split()[9:]
	totalsync=[]
	for n in pops:
		pophash=dict(zip(formID,n.split(":")))
		if pophash["GT"]=="./.":
			totalsync.append("0:0:0:0:0:0")
			continue
		altsum=0
		alleles=collections.defaultdict(int)
		if "ADr" in pophash and "ADf" in pophash:
			popALT1=pophash["ADr"].split(",")
			popALT2=pophash["ADf"].split(",")
		else:        
			popALT=pophash["AD"].split(",")
		for i in range(len(ALT)):
			if "ADr" in pophash and "ADf" in pophash:
				if len(popALT1) <= i:
					continue 
				altsum+=int(popALT1[i])
				altsum+=int(popALT2[i])
				alleles[ALT[i]]+=int(popALT1[i])
				alleles[ALT[i]]+=int(popALT2[i])
			else:
				if len(popALT) <= i:
					continue    
				altsum+=int(popALT[i])
				alleles[ALT[i]]=int(popALT[i])
		alleles[ref]=int(pophash["DP"])-altsum
		totalsync.append(sync_format(alleles))
	print chr+"\t"+pos+"\t"+ref+"\t"+"\t".join(totalsync)
