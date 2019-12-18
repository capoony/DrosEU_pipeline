import sys
from collections import defaultdict as d
import gzip
from optparse import OptionParser, OptionGroup
from bisect import bisect_left as BL

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="\npython %prog --badcov SNP_BS.txt --indel indels.txt --te flybase-te.gff --step 10000 --window 50000 --chromosomes 2L:23513712,2R:25286936 --output truewindows"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script generates an output (--out) which contains the number of sites that passed the quality criteria during SNP calling (--badcov), are not located within TE's (--te) and that are not located close to InDels (--indel) in sliding windows with a given window- (--window) and step-size (--step). This file can be used as on input for the calculation of population genetics parameters with PopGen-var.py.
You will need to pass the length of the chromosomes that should be considered in the output (--chromosomes )
Note: this script is quite memory hungry. Thus, I would run the script for the different chromsomal arms separately unless you have >64 GB of RAM memory

""") 
#########################################################   CODE   #########################################################################

parser.add_option("--badcov", dest="bs", help="bad coverage file with the extension BS from SNP calling",default="NA")
parser.add_option("--indel", dest="ind", help="indel file",default="NA")
parser.add_option("--te", dest="te", help="te database in gff fileformat from flybase or from Repeatmasker",default="NA")
parser.add_option("--include", dest="include", help="List of sites to restrict the analysis on",default="NA")
parser.add_option("--step", dest="stp", help="the step-size, e.g. 100000")
parser.add_option("--window", dest="win", help="the window-size, e.g. 100000")
parser.add_option("--output", dest="out", help="output file")
parser.add_option("--chromosomes", dest="chrom", help="the name and length of chromosomes to be considered, Name and length Need to be saparated by a colon and different chromosomes need to be separated by a comma, e.g. 2L:23513712,2R:25286936")

  
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

################# sliding windows ##############

exclude=d(lambda:d(int))
miss=d(lambda:d(list))
#code={"2L":23513712,"2R":25286936,"3L":28110227,"3R":32079331,"X":23542271,"4":32079331,"Y":3667352}
code={}
for Chrom in options.chrom.split(","):
    N,L=Chrom.split(":")
    code[N]=int(L)


WINDOW=map(int,options.win.split(","))
STEP=map(int,options.stp.split(","))

if options.ind!="NA":
    for l in load_data(options.ind):
        if len(l.split())<2:
               continue
        C,P=l.split()
        if C not in code:
            continue
        exclude[C][int(P)]
    
    print "Read InDel file done"

if options.te!="NA":
    for l in load_data(options.te):
        if l.startswith("#"):
            continue
        C=l.split()[0]
        if C not in code:
            continue
        S,E=map(int,l.split()[3:5])
        for i in range(S,E+1):
            exclude[C][i]
    
    print "Read TE file done"

if options.include!="NA":
    INC=d(lambda:d(int))
    for l in load_data(options.include):
        a=l.rstrip().split()
        if a[0] not in code:
            continue
        INC[a[0]][int(a[1])]
    
    for C,v in code.items():
        for i in range(int(v)):
            if i+1 not in INC[C]:
                exclude[C][i]
        
    INC=""
    print "INCLUDE Sites file done"
    
co=1
for l in load_data(options.bs):
    C=l.split()[0]
    if C not in code:
        continue
    if co%1000000==0:
        print co,"Positions done"
    co+=1
    a=l.split()
    for i in range(len(a[2])):
        if int(a[1]) in exclude[a[0]]:
            miss[a[0]][i].append(int(a[1]))
            continue
        if a[2][i]=="1":
            miss[a[0]][i].append(int(a[1]))
    if int(a[1]) in exclude[a[0]]:
        del exclude[a[0]][int(a[1])]

print "Read bad coverage file done"

for k,v in sorted(exclude.items()):
    for p,v1 in sorted(v.items()):
        for i in miss[k]:
            miss[k][i].append(p)

print "appending TE and InDel info done"

for k,v in miss.items():
    for I,v1 in v.items():
        miss[k][I]=sorted(set(v1))

print "sorting dictionary done"

for j in range(len(WINDOW)):
    out=open(options.out+"-"+str(WINDOW[j])+"-"+str(STEP[j])+".txt","w")
    for CHR in sorted(miss.keys()):
        window=WINDOW[j]
        step=STEP[j]
        print window,step,"started"
        steps={}
        stop=code[CHR]
        start=0
        end=window
        get=[]
        for i in range(len(miss[CHR])):
            w=int(WINDOW[j])
            ### count how many sites to be removed are located in the window defined by "start" and "end"
            get.append(w-(BL(miss[CHR][i],end+1)-BL(miss[CHR][i],start)))
        print "window",end,"done"
        while(end<=stop):
            out.write(CHR+"\t"+str((start+end)/2.0)+"\t"+"\t".join(map(str,get))+"\n") 
            start=start+step
            end=start+window
            get=[]
            for i in range(len(miss[CHR])):
                w=int(WINDOW[j])
                get.append(w-(BL(miss[CHR][i],end+1)-BL(miss[CHR][i],start)))
            print "window",end,"done"
        out.write(CHR+"\t"+str((start+end)/2.0)+"\t"+"\t".join(map(str,get))+"\n") 
