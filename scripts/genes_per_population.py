 #!/usr/bin/python

#### By Maria Bogaerts
#### Genes identified by pool-hmm per population
#### Previously, by using awk and bash commands, we separated one gene by file with all populations identified for that specific gene

import sys
import os

# Final file
summary_file = open("genes_per_pop.txt", 'w')

# Define samples per population
pop1 = ["11"]
pop2 = ["12"]
pop3 = ["13"]
pop4 = ["14","15"]
pop5 = ["16","18"]
pop6 = ["19","20","21","22"]
pop7 = ["23","24"]
pop8 = ["25"]
pop9 = ["26"]
pop10 = ["27"]
pop11 = ["28","29"]
pop12 = ["30"]
pop13 = ["31","32"]
pop14 = ["33"]
pop15 = ["34","35"]
pop16 = ["36","37"]
pop17 = ["38"]
pop18 = ["39","41"]
pop19 = ["42"]


for arg in sys.argv[1:]:
    populations = []
    te_name = arg.split("_")[0]
    with open(arg) as f:
        for line in f:
            s = line.split("-")[0]
            print(s)
            if s in pop1 and "pop1" not in populations:
                populations.append("pop1")
            elif s in pop2 and "pop2" not in populations:
                populations.append("pop2")
            elif s in pop3 and "pop3" not in populations:
                populations.append("pop3")
            elif s in pop4 and "pop4" not in populations:
                populations.append("pop4")
            elif s in pop5 and "pop5" not in populations:
                populations.append("pop5")
            elif s in pop6 and "pop6" not in populations:
                populations.append("pop6")
            elif s in pop7 and "pop7" not in populations:
                populations.append("pop7")
            elif s in pop8 and "pop8" not in populations:
                populations.append("pop8")
            elif s in pop9 and "pop9" not in populations:
                populations.append("pop9")
            elif s in pop10 and "pop10" not in populations:
                populations.append("pop10")
            elif s in pop11 and "pop11" not in populations:
                populations.append("pop11")
            elif s in pop12 and "pop12" not in populations:
                populations.append("pop12")
            elif s in pop13 and "pop13" not in populations:
                populations.append("pop13")
            elif s in pop14 and "pop14" not in populations:
                populations.append("pop14")
            elif s in pop15 and "pop15" not in populations:
                populations.append("pop15")
            elif s in pop16 and "pop16" not in populations:
                populations.append("pop16")
            elif s in pop17 and "pop17" not in populations:
                populations.append("pop17")
            elif s in pop18 and "pop18" not in populations:
                populations.append("pop18")
            elif s in pop19 and "pop19" not in populations:
                populations.append("pop19")
    print >> summary_file, str(te_name) + "\t" + str(len(populations)) + "\t" + ",".join(populations)
