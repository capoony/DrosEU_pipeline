import sys
from collections import Counter,defaultdict
from operator import itemgetter
import math
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun 

#########################################################   HELP   #########################################################################
usage="python %prog --fst input.fst --meta coord_EU.txt --output out"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________
Requirements:

R version 2.15
python module rpy2

Calculates geographic euclidean distances taking spherical curvature into account for all samples and test for correlation with genetic distances using Mantel tests in R

""") 
#########################################################   CODE   #########################################################################

parser.add_option("--output", dest="out", help="The path and name used for the output-file(s)")
parser.add_option("--meta", dest="met", help="Text file containing the coordinates and sample ID's")
parser.add_option("--fst", dest="fst", help="Text file containing the aversage pairwise FSTs")
parser.add_option("--color", dest="c", help="optional color",default="blue")


(options, args) = parser.parse_args()
parser.add_option_group(group)

def distance_on_unit_sphere(v,w):
    from math import sin, cos, sqrt, atan2, radians
    
    R = 6373.0
    lat1,lon1,lat2,lon2=map(radians,map(float,[v["lat"],v["long"],w["lat"],w["long"]]))
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c
    return distance


metadata=open(options.met,"r")
data=defaultdict(lambda:defaultdict(float))

r('pdf("'+options.out+'.pdf",width=10,height=5)')

for l in metadata:
    if l.startswith("ID"):
        continue
    name,ID,y,x=l.rstrip().split()[:4]
    data[int(ID)-2]["long"]=x
    data[int(ID)-2]["lat"]=y

#print data  
fstlist=[]
distlist=[]

O=open(options.out+".txt","w")
O.write("A\tB\tFST\tDIST\tC\n")
for l in open(options.fst):
    a,b,fst=l.split()
    fstlist.append(float(fst))
    O.write(a+"\t"+b+"\t"+fst+"\t"+str(distance_on_unit_sphere(data[int(a)],data[int(b)]))+"\t"+options.c+"\n")
    distlist.append(distance_on_unit_sphere(data[int(a)],data[int(b)]))
    
r.assign("fst",robjects.vectors.FloatVector(fstlist))
r.assign("dist",robjects.vectors.FloatVector(distlist))

r('plot(dist,fst,xlab="Kilometers",col=rgb(0,0,0,0.2),ylab="FST",pch=16)')
r('reg=lm(fst~dist)')
r('abline(reg,col="'+options.c+'")')
pval=str(list(r('cbind(coef(summary(reg))[,4][2])[1]'))[0])
rsqu=str(list(r('summary(reg)$r.squared '))[0])
r('legend("topleft",legend=paste("P-val: ",'+pval+',"\nR^2: ",'+rsqu+',sep=""))')
r('dev.off()')