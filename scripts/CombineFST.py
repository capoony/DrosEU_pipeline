import sys
from collections import defaultdict as d
import random
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --diff intron.fst > output.fst  "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Calculate average FST across all loci.
""") 
#########################################################   CODE   #########################################################################

parser.add_option("--diff", dest="diff", help="The output of PopGen_diff.py")
parser.add_option("--stat", dest="stat", help="The statistic calculated by PopGen_diff.py",default=0)

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
stat=int(options.stat)

FST=d(list)
for l in load_data(options.diff):
    a=l.rstrip().split()
    for comp in a[2:]:
        ID,F=comp.split("=")
        f=F.split(",")[stat]
        if f=="NA":
            continue

        FST[ID].append(float(f))
        
for k,v in sorted(FST.items()):
    print "\t".join(k.split(":"))+"\t"+str(sum(v)/len(v))
    
