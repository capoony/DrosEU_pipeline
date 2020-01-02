import sys
from collections import defaultdict as d
from optparse import OptionParser,OptionGroup
import math

#Author: Martin Kapun
#version: 1.0

#########################################################   HELP   #########################################################################
#print
parser = OptionParser()
group=OptionGroup(parser,
"""
H E L P:
____________

Subset a set of SNPs in --input, which are defined by Chromosome and Position in the first two columns. Use information for Recombination Rates of Comeron et al. (2012) and Inversion breakpoints from Corbett-Detig et al. (2014) as criteria and define minimum recombination rates (--r) and minimum Distance (--D) as Criteria to retain SNPs

""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="inp", help="Amplicon AF")
parser.add_option("--inv", dest="inv", help="Inversion breakpoint coordinates")
parser.add_option("--RecRates", dest="RR", help="Recombination Rates from Comeron et al. 2012")
parser.add_option("--r", dest="r", help="minimum recombination Rate/Mb")
parser.add_option("--D", dest="D", help="minimum Distance to inversion breakpoints")

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

dist=int(options.D)
r=int(options.r)

inh=d(lambda:d(str))
for l in open(options.inv,"r"):
    C,P=l.split()[2].split(":")
    S,E=map(int,P.split(".."))
    for i in range(S-dist,E+1+dist):
        inh[C][str(i)]


for l in load_data(options.RR):
    if l.startswith("Chrom"):
        continue
    a=l.rstrip().split()
    if a[2]=="NA":
        inh[a[0]][str(a[1])]
        continue
    if float(a[2])<r:
        inh[a[0]][str(a[1])]

for l in load_data(options.inp):
    C,P=l.split()[:2]
    if P in inh[C]:
        continue
    print l.rstrip()
