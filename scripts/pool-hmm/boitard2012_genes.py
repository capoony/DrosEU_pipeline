#!/usr/bin/python

#### By Maria Bogaerts
#### Overlap between genes presented in Boitard et al. 2012 and the present analysis

import sys
import os

# Genes identified by pool-hmm algorithm and presented in Boitard et al. 2012
genes_boitard = ["CG17636", "RhoGAP1A", "CG17707", "SP71", "CG3038", "CG2995", "cin", "CG13377", "CG13376", "ewg", "CG3777", "CG13375", "CG12470", "Or1a", "CG32816", "y", "ac", "sc", "l(1)sc", "pcl", "ase", "Cyp4g1", "Exp6", "CG13373", "CG18275", "CG32817", "CG18166", "CG3176", "CG18273", "CG3156", "CG17896", "CG17778", "svr", "arg", "elav", "CG4293", "Appl", "su(s)", "CG13367", "Roc1a", "Suv4-20", "skpA", "sdk", "CG13362", "CG13361", "CG5254", "CG5273", "RpL22", "fz3", "eIF4E-7", "CG34320", "CG11378", "CG11384", "CG11379", "CG14627", "CG14626", "CG11380", "CG14625", "CG11381", "CG14624", "CG11382", "CG11398", "CG3638", "CG11403", "A3-3", "CG32812", "DAAM", "CG18091", "fs(1)N", "CG11409", "CG11412", "CG11418", "Tsp2A", "CG12773", "CG11417", "png", "CG14770", "CG3056", "SNF1A", "CG3719", "CG32813", "CG11448", "futsch", "futsch", "Gr2a", "CG14785", "CG14786", "CG14787", "l(1)G0431", "O-fut2", "CG14777", "CG32808", "CG14778", "pck", "CG14780", "Rab27", "CG14782", "sta", "Nmdar2", "CG14795", "CG32810", "Adar", "CG32806", "CG14801", "CG14812", "deltaCOP", "CG14814", "MED18", "CG14815", "CG14803", "CG14816", "CG14804", "CG14817", "CG14805", "CG14818", "CG14806", "trr", "mRpL16", "arm", "CG32803", "CG32801", "Edem1", "mip130", "CG17766", "csw", "ph-d", "ph-p", "CG3835", "Pgd", "bcn92", "wapl", "Cyp4d1", "CG3630", "CG3621", "Cyp4d14", "Mct1", "CG18031", "msta", "Vinc", "CG14052", "Tlk", "CG3033", "mof", "CG3016", "CG16721", "Ca-alpha1T", "CG1958", "CG1677", "CG2059", "unc-119", "CG11368", "CG32719", "CG10777", "CG10778", "RpS14a", "RpS14b", "CG1530", "l(1)G0193", "CG1531", "CG15332", "CG17255", "CG2889", "CG2887", "PPP4R2r", "CG32687", "Cyp4g15", "CG1749", "Spase25", "CG33235", "CG32666", "CG32666", "CG1572", "PGRP-SA", "RpII215", "CG11699", "l(1)G0237", "CG11697", "CG11696", "e(y)2", "CG11695", "nod", "CG1561", "rho-4", "CG2533", "cac", "gd", "tsg", "CG18130", "fw", "sno", "REG", "mew", "hiw", "CG5541", "PGRP-LE", "sd", "CG8509", "Ranbp16", "Stim", "CG8924", "CG8928", "CG15603", "CG15604", "CG15814", "CG6506", "CG32554", "CG32557", "CG6762", "Arp8", "CG6769", "mnb", "Sh", "CG15373", "l(1)G0003", "CG6540", "CG6617", "Ing3", "CG6659", "fu", "CG6696", "Grip84", "car", "Tao-1", "CG14218", "CG14204", "CG11566", "stg1", "unc", "CG15445", "CG34120", "waw", "bbx", "slgA", "Hlc", "mst"]

# File with the FlyBase code (FBgn) for each gene
translation_file = open("boitard_2012_genes_FBgn.txt", 'r').readlines()

# Final file
end_file=open(sys.argv[1],'w')

for t in translation_file:
    fbgn = t.split("\t")[0]
    gene_name = t.split("\t")[1]
    gene = gene_name.split("\r")[0]
    gene_list = []
    for arg in sys.argv[2:]:
        print(arg)
        name = arg.split(".")[0]
        strain = name.split("_")[2]
        with open(arg) as f:
            for line in f:
                gene_strain = line.split("\n")[0]
                if fbgn == gene_strain:
                    if strain not in gene_list:
                        gene_list.append(strain)
    print >> end_file, gene + "\t" + ",".join(gene_list)
