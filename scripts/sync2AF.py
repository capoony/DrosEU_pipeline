import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input input.sync --output out"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Calculates allele frequencies of the major allele for all libraries and stores it in a matrix, where rows are libraries and columns are SNPs.

""")
#########################################################   CODE   #########################################################################

parser.add_option("--output", dest="out", help="The path and name used for the output-file(s)")
parser.add_option("--input", dest="inp", help="A sync file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def getAlleles(v):
	''' get major and minor allele from all populations'''
	from collections import Counter
	allnuc=""
	for po in v:
		allnuc+=sync2string(po)
	if allnuc=="":
		return "na"
	toptwo=Counter(allnuc).most_common()[:2]
	if len(toptwo)<2:
		return "na"
	return zip(*toptwo)

def sync2string(x):
	''' convert sync format to string of nucleotides  where x is a string in sync format '''
	string=""
	alleles=["A","T","C","G"]
	ah=dict(zip(alleles,map(int,x.split(":")[:4])))
	for k,v in ah.items():
		string+=v*k
	return string

def sync2freqh(x):
	''' convert string in SYNC format to dictionary of freqencies where x is a string in sync format'''
	from collections import defaultdict as d

	nuc=["A","T","C","G"]
	counts=map(int,x.split(":")[:4])
	CO=dict(zip(*[nuc,counts]))
	h=d(float)
	if sum(CO.values())==0:
		return "NA"
	for k,v in CO.items():
		h[k]=v/float(sum(CO.values()))

	return h,sum(CO.values())

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

freqh=d(list)
o1=open(options.out+"_pos.txt","w")
o2=open(options.out+"_freq.txt","w")

for l in load_data(options.inp):
	C,P,R=l.split()[:3]
	Pops=l.split()[3:]
	AA=getAlleles(Pops)[0]
	if AA=="na" or len(AA)==1:
		continue
	A,a=AA
	o1.write(C+"\t"+P+"\t"+A+"/"+a+"\n")
	for i in range(len(Pops)):
		FR=sync2freqh(Pops[i])
		if FR=="NA":
			freqh[i].append("NA")
		else:
			freqh[i].append(FR[0][a])

for k,v in sorted(freqh.items()):
	o2.write("\t".join(map(str,v))+"\n")
