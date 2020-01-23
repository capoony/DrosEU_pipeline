import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import math

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="\npython %prog --input input.sync --step 10000 --window 50000 --min-count 2 --min-cov 10 --pool-size 80,80,80 --min-sites-frac 0.75 --site-counts counts_10000-50000.txt --output out"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,"""
H E L P:
____________

This script calculates Pool-Seq corrected Population genetic parameters pi, Theta and Tajima's D using the corrections in Kofler et al. 2011 (PLoS One). This script requires a sync file as an input (--input) where every column represents an individual library. The script produces output files for each statistic jointly. The path and the prefix of these files need to be provided (--output). Average statistics are calculated in windows, were step-sizes are defined by --step and window-sizes by --window. If you want to calculate non-overlapping windows you need to make sure that window-sizes are of the same size as step-sizes.

The scripts needs a mandatory sitecounts file (--sitecounts), which provides the number of "proper sites" (sites, which are not located within TE's, in the proximity to InDels and which did not fail the qualiyt criteria during SNP calling) for every window. The window- and step-sizes in this file must be correspond to the parameters in this script.

Additional optional parameters are "minimum allele count" (--min-count, default=2) where only alleles above this thresholds are considered, "minimum coverage" (--min-cov, default=10) where only sites above this thresholds are considered and "Minimum fraction of covered sites" (--min-sites-frac, default=0.75) where only windows are conssidered whose porportion of "proper sites" (see above) is equal or larger than the threshold.

!!!!!!!

Note: this script is quite memory hungry. Thus, I would run the script for the different chromsomal arms separately unless you have >64 GB of RAM memory

!!!!!!!

""") 
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="sync", help="sync file with all SNPs")
parser.add_option("--step", dest="stp", help="the step-size, e.g. 100000",default="NA")
parser.add_option("--window", dest="win", help="the window-size, e.g. 100000",default="NA")
parser.add_option("--min-count", dest="min", help="minimum allele count, default=2",default=2)
parser.add_option("--min-cov", dest="mincov", help="minimum coverage,default=10",default=10)
parser.add_option("--min-sites-frac", dest="mins", help="minimum fraction of proper sites in a given window, default=75%",default=0.75)
parser.add_option("--pool-size", dest="pools", help="comma-separated list of Pool-size in the same order as in the sync file")
parser.add_option("--sitecount", dest="sites", help="list of sites per window")
parser.add_option("--output", dest="out", help="The path and name used for the output-file(s)")


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

def average(x):
    return sum(x)/float(len(x))

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    alleles=["A","T","C","G"]
    ah=dict(zip(alleles,map(int,x.split(":")[:4])))
    for k,v in ah.items():
        string+=v*k
    
    return string

def string2freqh(x):
    ''' convert string of nucleotides to dictionary of freqencies'''
    from collections import defaultdict as d
    
    counts=d(int)
    alleles=["A","T","C","G"]
    for a in alleles:
        if x.count(a)==0:
            continue
        counts[a]=x.count(a)
    
    h=d(float)
    for k,v in counts.items():
        h[k]=v/float(len(x))
        
    return h,len(x)

def sync2counth(x):
    ''' convert string in sync format to dictionary of counts where x is a string in sync format'''
    nuc=["A","T","C","G"]
    counts=map(int,x.split(":")[:4])
    return dict(zip(*[nuc,counts]))

def binom(x,y):
    ''' calculate the binomial coefficient for x over y'''
    import math
    from decimal import Decimal
    if y == x:
        return 1
    elif y == 1:
        return x
    elif y > x:
        return 0
    else:               
        a = math.factorial(x)
        b = math.factorial(y)
        c = math.factorial(x-y)
    try:
        div = a / (float(b) * float(c))
    except:
        div = Decimal(a) / (Decimal(b) * Decimal(c))
    return float(div)

################ Theta a la Kofler 2011 #######################

def ThetaCorr(M,n,b):
    ''' as in Kofler et al. 2011; Page 7'''
    T=0.0
    for m in range(int(b),int(M-b+1)):
        K=0.0
        for k in range(1,int(n)):
            BE=binom(M,m)
            T1=((k/float(n))**m)
            T2=((((n)-k)/float(n))**(M-m))
            K+=(1/float(k))*BE*T1*T2
        T+=K
    if T==0:
        return "NA"
    return T

################# Pi a la Kofler et al. 2011 #####################

def pi(x,M):
    ''' calculate pi on a SNP-wise basis. where x is a vector of all allelefreqs and n is the samplesize)
    x is a list of allele frequencies
    n is the total coverage'''
    if M==0:
        return "NA"
    else:
        corr=(M)/float(M-1)
        freqsum=sum([y**2 for y in x])
        return (1-freqsum)*corr 

def PiCorr(M,n,b):
    ''' see Kofler et al. 2011; Page 7
    x is a list of allele counts
    n is the poolsize
    b is the minor allele threshold'''
    T=0.0
    for m in range(int(b),int(M-b+1)):
        UT=(2*m*(M-m)/float(M*(M-1)))
        K=0.0
        for k in range(1,int(n)):
            BE=binom(M,m)
            T1=((k/float(n))**m)
            T2=((((n)-k)/float(n))**(M-m))
            K+=(1/float(k))*BE*T1*T2
        T+=UT*K
    if T==0:
        return "NA"
    return T

################# Tajima's D a la Kofler et al. 2011 #############
# N= poolsize
# C= coverage list from SNPs

def an(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 4'''
    toret=0.0
    for i in range(1,int(n)):
        toret+=1/float(i)
    return toret
    
def bn(n):
    '''Kofler et al. 2011 Appendix: P.7 Formula 1'''
    toret=0.0
    for i in range(1,int(n)):
        toret+=1/float((i**2))
    return toret

def fstar(n,AN):
    '''Kofler et al. 2011 Appendix: P.6 Formula 1'''
    return (n - 3)/(AN*(n - 1) - n)
   
def alphas(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 2'''
    AN=an(n)
    FS=fstar(n,AN)
    t1=(FS**2)*(AN-(n/(n-1)))
    st1=AN * ( (4*(n+1)) / ((n-1)**2) )
    st2=2 * ((n+3)/(n-1))
    t2=FS * (st1-st2)
    t3=AN * ( (8*(n+1))/(n*((n-1)**2)) )
    t4= ((n**2)+n+60)/(3*n*(n-1))
    #print "AS:",t1 + t2 - t3 + t4
    return t1 + t2 - t3 + t4

def betas(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 3'''
    AN=an(n)
    BN=bn(n)
    FS=fstar(n,AN)
    
    t1 = (FS**2) * (BN - ((2*(n-1)) /((n-1)**2)))
    st1= BN * (8/(n-1))
    st2= AN * (4/(n*(n-1)))
    st3= ((n**3)+12*(n**2)-35*n+18)/(n*((n-1)**2))
    t2 = FS*(st1-st2-st3)
    t3 = BN * (16/(n*(n-1)))
    t4 = AN * (8/((n**2)*(n-1)))
    st4= 2*(n**4+ 110*(n**2)-255*n+126)
    st5= 9*(n**2)*((n-1)**2)
    t5 = st4/st5
    #print "BS:",t1 + t2 - t3 + t4 + t5
    return t1 + t2 - t3 + t4 + t5

def nbase (np,M):
    pij=pijmatrix(3*np,np)
    nb=0.0
    minj=max([M,np])
    for k in range(1,minj+1):
        nb+=k*pij[M][k]
    #print nb
    return nb

def pijmatrix(MC,np):
    from collections import defaultdict as d
    jb=min([MC,np])
    Matrix=d(lambda:d(float))
    #print np
    Matrix[0][0]=1.0
    for i in range(1,MC+1):
        for j in range(1,jb+1):
            t1= ((1+np-j)/float(np))*Matrix[i-1][j-1]
            t2= (j/float(np))*Matrix[i-1][j]
            pij=t1+t2
            Matrix[i][j]=pij
    
    return Matrix
   
def div(C,theta,avn):
    import math
    AS=alphas(avn)
    BS=betas(avn)
    #print avn,AS,BS,len(C),theta
    div1=(AS/C)*theta + BS*(theta**2)
    return math.sqrt(abs(div1))

def D(pi,theta,C,avn):
    '''Kofler et al. 2011 Appendix: P.7 Formula 3'''
    Div=div(C,theta,avn)
    if pi-theta==0:
        return 0
    else:
        return (pi-theta)/Div
        
################# sliding windows ##############

def slider(x,window,step):
    '''return list of items in sliding windows, x must be a list of numbers'''
    steps={}
    stop=max(x)
    start=0
    end=window
    get=[y for y in set(x) if y>start and y<=end]
    while(end<=stop):
        if get!=[]:
            steps[(start+end)/2.0]=get
        start=start+step
        end=start+window
        get=[y for y in set(x) if y>start and y<=end]
    if get!=[]:
        steps[(start+end)/2.0]=get  
    return steps

################ read parameters ##########


window=map(int,options.win.split(","))
step=map(int,options.stp.split(","))
mincount=int(options.min)
mincov=int(options.mincov)
mins=float(options.mins)
pools=map(int,options.pools.split(","))

OPl=[]
OTl=[]
ODl=[]

for i in range(len(window)):
    OPl.append(open(options.out+"_"+str(window[i])+"_"+str(step[i])+".pi","w"))
    OTl.append(open(options.out+"_"+str(window[i])+"_"+str(step[i])+".th","w"))
    ODl.append(open(options.out+"_"+str(window[i])+"_"+str(step[i])+".D","w"))

PopGendata=d(lambda:d(lambda:d(lambda:d(float))))
avndict={}
Pcorrhash={}
Tcorrhash={}
count=1

for l in load_data(options.sync):
    a=l.split()
    
    ## counter
    if count%10000==0:
        print count,"SNPs processed"
    count+=1

    ## find two major alleles
    pops=a[3:]
    alleles=getAlleles(pops)
    
    for i in range(len(pops)):
        
        ## set up dictionary
        PopGendata[a[0]][int(a[1])][i]["P"]
        PopGendata[a[0]][int(a[1])][i]["T"]
        
        ## calculate pi and theta
        
        ## get frequencies and counts and remove alleles that are not among the two major alleles
        xs=sync2string(pops[i])
        Px=dict([(k,v) for k,v in sync2counth(pops[i]).items() if k in alleles and v>=mincount])
        px=dict([(k,v) for k,v in string2freqh(xs)[0].items() if k in Px.keys()])

        if len(px)>1:
            
            ## get coverage and poolsize
            M=sum(Px.values())
            n=pools[i]
            
            ## calculate pi
            ID=str(M)+"_"+str(n)  
            if ID not in Pcorrhash:
                Pcorrhash[ID]=PiCorr(M,n,mincount)
            pix=pi(px.values(),M)/Pcorrhash[ID]
            
            PopGendata[a[0]][int(a[1])][i]["P"]=pix
            PopGendata[a[0]][int(a[1])][i]["M"]=M

            ## calculate Theta
            if ID not in Tcorrhash:
                Tcorrhash[ID]=ThetaCorr(M,n,mincount)
                
            PopGendata[a[0]][int(a[1])][i]["T"]=1/Tcorrhash[ID]

for Chrom,Values in sorted(PopGendata.items()):
    
    for W in range(len(window)):
        
        ## obtain site counts for windows:
        SITES=options.sites.split(",")
        covh=d(lambda:d(lambda:d(int)))
        for l in open(SITES[W],"r"):
            a=l.split()
            C,P=a[:2]
            popcov=a[2:]
            for I in range(len(popcov)):
                covh[C][int(float(P))][I]=int(popcov[I])
        if Chrom not in covh:
            continue
        
        print Chrom,"chromosome: calculation started for window-size",window[W],"and step-size",step[W]
        ## generate a dictionary with bins according to windows and stepsizes containing the corresponding positions
        bins=slider(Values.keys(),window[W],step[W])
        for Bin,Items in sorted(bins.items()):
            if Items==[]:
                continue
            print Bin,"window: calculation started"
            ph=d(float)
            th=d(float)
            mh=d(list)
            ## loop through positions in bin
            for i in Items:
                ## loop through samples:
                for j in range(len(pools)):
                    ph[j]+=Values[int(i)][j]["P"]
                    th[j]+=Values[int(i)][j]["T"]
                    if "M" not in Values[int(i)][j]:
                        continue
                    mh[j].append(Values[int(i)][j]["M"])
            
            plist=[]
            tlist=[]
            dlist=[]
            for i in range(len(pools)):
                
                ## ignore if no covered sites
                if len(mh[i])==0:
                    dlist.append("NA")
                    plist.append("NA")
                    tlist.append("NA")
                    continue
                
                AvCov=int(average(mh[i]))
                
                ## test if coverage larger than minimum coverage threshold
                if AvCov<mincov:
                    dlist.append("NA")
                    plist.append("NA")
                    tlist.append("NA")
                    continue
                
                ## test if site count at least min-sites-fraction*windowsize:
                if covh[Chrom][Bin][i]<window[W]*mins:
                    dlist.append("NA")
                    plist.append("NA")
                    tlist.append("NA")
                    continue
                
                ID=str(AvCov)+":"+str(pools[i])
                if ID not in avndict:
                    avndict[ID]=nbase(pools[i],AvCov)
                
                ## calculate average pi and Theta   
                Pi=ph[i]/covh[Chrom][Bin][i]
                Theta=th[i]/covh[Chrom][Bin][i]
                
                ## calculate Tajima's D
                dlist.append(D(Pi,Theta,covh[Chrom][Bin][i],avndict[ID]))
                plist.append(Pi)
                tlist.append(Theta)
            
            ## write output
            OPl[W].write(Chrom+"\t"+str(Bin)+"\t"+"\t".join(map(str,plist))+"\n")
            OTl[W].write(Chrom+"\t"+str(Bin)+"\t"+"\t".join(map(str,tlist))+"\n")
            ODl[W].write(Chrom+"\t"+str(Bin)+"\t"+"\t".join(map(str,dlist))+"\n")
                    
