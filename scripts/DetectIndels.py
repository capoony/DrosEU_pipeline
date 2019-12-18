import sys
import collections
import re
from optparse import OptionParser, OptionGroup


#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --mpileup input.mpileup  --minimum-count 10 --mask 5 > output_indels.txt"
parser = OptionParser(usage=usage)
helptext="""    

H E L P :
_________

Identify sites that are located close to InDels with a minimum count (--minimum-count) across all libraries and print them if they are within a predefined distance (--mask) to an InDel.
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("-m","--mpileup", dest="m", help="An mpileup file")
parser.add_option("--minimum-count", dest="i", help="The minimum count of an indel across all samples pooled")
parser.add_option("--mask", dest="k", help="The number of basepairs masked at an InDel in either direction (up- and downstreams)")
parser.add_option_group(group)
(options, args) = parser.parse_args()

######################################### functions #############################################
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


def keywithmaxvalue(d):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=collections.defaultdict(list)
    for k,v in d.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]


def splitter(l, n):
    ''' This generator function returns equally sized cunks of an list'''
    #credit: Meric Lieberman, 2012
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

def extract_indel(l,sign):
    ''' This function returns an Indel from a sequence string in a pileup'''
    position = l.index(sign)
    numb =""
    i = 0
    while True:
        if l[position+1+i].isdigit():
            numb+= l[position+1+i]
            i+=1
        else:
            break
        
    seqlength = int(numb)
    sequence = l[position:position+i+1+seqlength]
    indel=sequence.replace(numb,"")

    return sequence,indel

##################################### data #################################################

data=load_data(options.m)

##################################### main code #################################################

# parse mpileup and store alternative alleles:
allelehash=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(str)))

for line in data:

    k = line[:-1].split('\t')
    chrom,position,refbase = k[:3]
    div = list(splitter(k,3))
    libraries=div[1:]
    alleles=collections.defaultdict(int)
    # loop through libraries    
    for j in range(len(libraries)):
        
        nuc = libraries[j][1]
        qualities = libraries[j][2]
        
        
        # test if seq-string is empty
        if nuc=="*":
            continue
            
        # find and remove read indices and mapping quality string
        nuc = re.sub(r'\^.',r'',nuc)
        nuc = nuc.replace('$','')
        
        # find and remove InDels
        while "+" in nuc or "-" in nuc:
            if "+" in nuc:
                insertion,ins=extract_indel(nuc,"+")
                alleles[ins.upper()]+=nuc.count(insertion)
                nuc=nuc.replace(insertion,"")
            else:
                deletion,dele=extract_indel(nuc,"-")
                #alleles[dele.upper()]+=nuc.count(deletion)
                nuc=nuc.replace(deletion,"")
    test=0
    if len(alleles)==0:
        continue
    #print alleles
    for v in alleles.values():
        #print v
        if v>5:
            test=1
    if test==1:
        for i in range(int(position)-5,int(position)+6):
            print chrom+"\t"+str(i)
    