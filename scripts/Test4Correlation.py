import sys
from collections import defaultdict as d
from rpy2.robjects import r
import rpy2.robjects as robjects
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --meta MetaData.txt > Correlation.test  "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Test for correlations of Inversion or TE frequencies with geographic/temporal variables with linear (--stat lm) or generalized linear model with a binomial error structure (--stat glm) and account for potential spatial autocorrelation
""")
#########################################################   CODE   #########################################################################

parser.add_option("--meta", dest="meta", help="A text file containing the dependent and independent variables, see MetaData.txt")
parser.add_option("--stat", dest="stat", help="either linear model (lm) or generalized linear model (glm)",default="glm")

parser.add_option_group(group)
(options, args) = parser.parse_args()


data=open(options.meta,"r")
typ=options.stat

meta=d(list)

## parse input file
for l in data:
    a=l.rstrip().split()
    if l.startswith("ID"):
        header=a
        continue
    for i in range(len(a)):
        meta[header[i]].append(a[i])

r.assign("lat",robjects.vectors.FloatVector([float(x) for x in meta["lat"]]))
r.assign("lon",robjects.vectors.FloatVector([float(x) for x in meta["long"]]))
r.assign("alt",robjects.vectors.FloatVector([float(x) for x in meta["Altitude"]]))
r.assign("sea",robjects.vectors.StrVector(meta["Season"]))
r.assign("ID",robjects.vectors.StrVector(meta["ID"]))
r.assign("cov",robjects.vectors.FloatVector([float(x) for x in meta["coverage"]]))
factors=["lat","lon","alt","sea"]

r('basedat=data.frame("lon"=lon,"lat"=lat,"alt"=alt,"sea"=sea,"cov"=cov)')
plist=[]

## make header
for i in factors:
    plist.append(i+"_DF\t"+i+"_F\t"+i+"_P")
plist.append("MoransI\tM_P")
for i in factors:
    plist.append("EM_"+i+"_F\tEM_"+i+"_P")

print("Factor\t"+"\t".join(plist))

## generate neighbour weights matrix
r('library(spdep)')

## do stats for all data:
var=header[7:]

for I in var:
    plist=[]
    r.assign(I,robjects.vectors.StrVector(meta[I]))
    r('newdat=data.frame(basedat,"'+I+'"=as.numeric('+I+'))')
    r('nd2=na.omit(newdat)')
    if typ=="lm":
        r('Reg=lm('+I+'~lat+lon+alt+sea,data=nd2)')
    else:
        r('Reg=glm('+I+'~lat+lon+alt+sea,weights=cov,family=binomial,data=nd2)')

    ## test for significance of geographic factors
    if typ=="lm":
        r('Reg.a=anova(Reg)')
        Dflist=[str(x) for x in list(r('Reg.a$"Df"'))]
        Pvlist=[str(x) for x in list(r('Reg.a$"Pr(>F)"'))]
        Flist=[str(x) for x in list(r('Reg.a$"F value"'))]
    else:
        r('Reg.a=anova(Reg,test="Chisq")')
        Dflist=[str(x) for x in list(r('Reg.a$"Df"'))[1:]]
        Pvlist=[str(x) for x in list(r('Reg.a$"Pr(>Chi)"')[1:])]
        Flist=[str(x) for x in list(r('Reg.a$"Dev"')[1:])]
    for i in range(len(factors)):
        if typ=="lm":
            plist.append(Dflist[i]+"/"+Dflist[-1]+"\t"+Flist[i]+"\t"+Pvlist[i])
        else:
            plist.append(Dflist[i]+"\t"+Flist[i]+"\t"+Pvlist[i])

    ## test resiudals for autorcorrelation
    r('res=Reg$"residuals"')
    r('neid<-dnearneigh(coordinates(nd2[1:2]), 0, 11)')
    r('neid_w<-tryCatch(nb2listw(neid),error=function(e) e)')
    ## test if Matrix does not work
    if "error" in r('class(neid_w)'):
        Mp="na"
        Mi="na"
    else:
        r('M=moran.test(res,listw=neid_w)')
        Mp=str(list(r('M$"p.value"'))[0])
        Mi=str(list(r('M$"statistic"'))[0])
    plist.append(Mi+"\t"+Mp)

    ## calculate Spatial Error Model when Moran's I is significant
    if Mp=="na":
        plist.append("\t".join(["na"]*len(factors)*2))
    elif float(Mp)>0.05:
        plist.append("\t".join(["na"]*len(factors)*2))
    else:
        r('EM=errorsarlm('+I+'~lat+lon+alt+sea,data=nd2,listw=neid_w)')
        Eflist=[str(x) for x in list(r('summary(EM)$"Coef"[,1]'))[1:]]
        Pvlist=[str(x) for x in list(r('summary(EM)$"Coef"[,4]'))[1:]]
        for i in range(len(factors)):
            plist.append(Eflist[i]+"\t"+Pvlist[i])
    print(I+"\t"+"\t".join(plist))
