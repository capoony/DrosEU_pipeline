import sys
from collections import defaultdict as d
import gzip
from optparse import OptionParser, OptionGroup

#########################################################   HELP   #########################################################################
usage="""python %prog \
      --vcf input.vcf \
      --te TE.gff \
      --indel indel.txt \
      > output.vcf """
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Removes sites that are located in InDels or TE's from VCF input file
""") 

parser.add_option("--vcf", dest="vcf", help="The original VCF file")
parser.add_option("--te", dest="te", help="A GFF file containing the coordinates of TE's (optional) ")
parser.add_option("--indel", dest="indel", help="Text file containing the coordinates of sites to be removed (optional)")

parser.add_option_group(group)
(options, args) = parser.parse_args()


def load_data(x):
    if x=="-":
        y=sys.stdin
    elif x.endswith(".gz"):
        y=gzip.open(x,"r")
    else:
        y=open(x,"r")
    return y


exclude=d(lambda:d(str))

if options.indel is not None:
    for l in load_data(options.indel):
        if len(l.split())<2:
               continue
        C,P=l.split()
        exclude[C][int(P)]

if options.te is not None:
    for l in load_data(options.te):
        if l.startswith("#"):
            continue
        C=l.split()[0]
        S,E=map(int,l.split()[3:5])
        for i in range(S,E+1):
            exclude[C][i]

for l in load_data(options.vcf):
    if l.startswith("#"):
        print l.rstrip()
        continue
    C,P=l.split()[:2]
    if int(P) in exclude[C]:
        continue
    print l.rstrip()


