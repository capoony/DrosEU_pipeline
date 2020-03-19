import sys
import math
from optparse import OptionParser, OptionGroup
from collections import defaultdict as d
import gzip

#Author: Martin Kapun

#########################################################   HELP   #########################################################################

usage="python %prog --input input.sync --minimum-count 2 --minimum-cov 10 --pool-size 20,30,20 > output.fst"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
'''
H E L P:
____________

Calculate pairwise FST among all populations in the input sync file. Only consider site with an coverage > minimum-cov (--minimum-cov) and an allele-count > minimum-count (--minium-count)
''')
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="sync", help="An input file in the sync file format")
parser.add_option("--minimum-count", dest="b", help="minimum allele count",default=2)
parser.add_option("--minimum-cov", dest="m", help="minimum coverage",default=10)
parser.add_option("--pool-size", dest="s", help="comma separated list of chromosomes in the pools")

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

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na"
    alleles=["A","T","C","G"]
    ah=dict(zip(alleles,map(int,x.split(":")[:4])))
    for k,v in ah.items():
        string+=v*k
    return string

def string2freqh(x):
    ''' convert string of nucleotides to dictionary of freqencies where x is a string of nucleotides'''
    from collections import defaultdict as d
    if x=="" or x=="na":
        return "na","na"
    counts=d(int)
    nuc=["A","T","C","G"]
    for u in nuc:
        if x.count(u)==0:
            continue
        counts[u]=x.count(u)

	h=d(float)
	for k,v in counts.items():
		h[k]=v/float(sum(counts.values()))
    return h,sum(counts.values())

def sync2counth(x):
    ''' convert string in sync format to dictionary of counts where x is a string in sync format'''
    nuc=["A","T","C","G"]
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na"
    counts=map(int,x.split(":")[:4])
    return dict(zip(*[nuc,counts]))

def counth2freqh(x):
	''' calculate freq-list from counts in dictionary'''
	from collections import defaultdict as d

	h=d(float)
	for k,v in x.items():
		h[k]=v/float(sum(x.values()))
	return h,sum(x.values())

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
    return zip(*toptwo)[0]

def FST_wc(p1,p2,n1,n2):
    '''see Bahtia et al. (2013) Genome Research - eq. 6'''
    if n1+n2<=2:
        return "NA"
    nom1=(n1*n2)/(n1+n2)
    nom2=1/(n1+n2-2)
    nom3=n1*p1*(1-p1)+n2*p2*(1-p2)
    nom=2*nom1*nom2*nom3

    denom1=n1*n2/(n1+n2)
    denom2=(p1-p2)**2
    denom3=(2*(n1*n2/(n1+n2)))-1
    denom4=1/(n1+n2-2)
    denom5=n1*p1*(1-p1)+n2*p2*(1-p2)
    denom=denom1*denom2+denom3*denom4*denom5
    
    if denom==0:
        return "NA"

    FST=1-(nom/denom)
    if FST<0:
        return 0.0
    return FST

def mergeDict(x,y):
    return { k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y) }

#outdiff=gzip.open(options.out+".diff.gz","w")
#outpi=gzip.open(options.out+".pi.gz","w")
Sample=map(int,options.s.split(","))
test=0
corrhash={}
Co=1
for l in load_data(options.sync):
    if len(l.split())<2:
        continue
    if Co%100000==0:
        print Co,"Snps processed"
    Co+=1
    C,P,R=l.split()[:3]
    pops=l.split()[3:]
    alleles=getAlleles(pops)
    if len(alleles)==1:
        continue
    A,a=alleles
    poplen=range(len(pops))
    pil=[]
    stat=[]
    if test==1:
        test=0
        continue
    for i in poplen:
        xs=sync2string(pops[i])
        px,nx=string2freqh(xs)
        for j in poplen:
            if j<=i:
                continue
            ys=sync2string(pops[j])
            py,ny=string2freqh(ys)

            ## calculate FST Weir Cockerham
            if px=="na" or py=="na":
                stat.append(str(i+1)+":"+str(j+1)+"=NA")
                continue
            FSTw=FST_wc(px[A],py[A],float(Sample[i]),float(Sample[j]))
            if nx<int(options.m) or ny <int(options.m):
                stat.append(str(i+1)+":"+str(j+1)+"=NA")
            else:
                stat.append(str(i+1)+":"+str(j+1)+"="+str(FSTw))

    print C+"\t"+P+"\t"+"\t".join(map(str,stat))
