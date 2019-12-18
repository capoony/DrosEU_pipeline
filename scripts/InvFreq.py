import sys
from collections import defaultdict
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="\npython %prog --sync input.sync --inv inversion-markers.txt --meta MetaData.txt > inversion.freq"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script calculates population-specific inversions frequencies based on inversion-specific marker SNPs in Kapun et al. 2014""") 
#########################################################   CODE   #########################################################################

parser.add_option("--sync", dest="sync", help="sync file with all SNPs")
parser.add_option("--inv", dest="inv", help="Inversion specific SNPs with Inversion in first and Chromosome and position in second and third column")
parser.add_option("--meta", dest="pop", help="dataset containing the populations to be analyzed with full name in the first, and position in the sync file in the second column")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def synccounts(x):
	string=""
	alleles=["A","T","C","G"]
	ah=dict(zip(alleles,map(int,x.split(":")[:4])))
	for k,v in ah.items():
		string+=v*k
	return string
def av(x):
	if len(x)==0:
		return 0.0
	return sum(x)/len(x)


sync=open(options.sync,"r")
populations=open(options.pop,"r")
inv=open(options.inv,"r")


pops=[]
poppos=defaultdict(int)
invpos=defaultdict(str)
geodist=defaultdict(float)

## get position of target populations in full dataset

populations.readline()
for l in populations:
	fullid,pos=l.split()[:2]
	pops.append(fullid)
	poppos[fullid]=int(pos)

for l in inv:
	i,c,p,a=l.split()
	invpos[c+"_"+p]=[i,a]

invhash=defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

for l in sync:
	c,p=l.split()[:2]
	allpops=l.split()
	ID=c+"_"+p
	if ID not in invpos.keys():
		continue
	inv,allele=invpos[ID]
	for pop in pops:
		invhash[inv][pop]["inv"]+=synccounts(allpops[poppos[pop]]).count(allele)
		invhash[inv][pop]["tot"]+=len(synccounts(allpops[poppos[pop]]))
	invhash[inv][pop]["inv"]
	invhash[inv][pop]["tot"]

finalhash=defaultdict(list)
latitude=[]
test=0

print "Inv\t"+"\t".join(pops)

for inv,v in invhash.items():
	invfr=[]
	for pop in pops:
		invfr.append(v[pop]["inv"]/float(v[pop]["tot"]))
	print inv+"\t"+"\t".join(map(str,invfr))
