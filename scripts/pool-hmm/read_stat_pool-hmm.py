#!/usr/bin/python

#### By Maria Bogaerts
#### Transform Pool-HMM (Boitard et al. 2013) stat files results into bed files and remove heterochromatic regions (centromeres and telomeres)

# Usage: python read_stat_pool-hmm.py *.stat

import sys
import os

chromosomes = ["2L", "2R", "3L", "3R", "4", "X"]

final_file = open("summary_all.txt", 'w')

for arg in sys.argv[1:]:
    print(arg)
    name = arg.split(".")[0]
    strain = arg.split("-")[0]
    chromosome = name.split("-")[1]
    arg_final = open("summary_" + chromosome + "_" + strain + ".bed", 'w')
    with open(arg) as f:
        next(f)
        for line in f:
            start = int(line.split(" ")[0])
            end = int(line.split(" ")[1])
            score = line.split(" ")[2]
            score_final = score.split("\n")[0]
            if chromosome == "2L":
                if start >= 530000 and end <= 18870000:
                    print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                    print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
            if chromosome == "2R":
                if start >= 5982495 and end <= 24972477:
                    print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                    print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
            if chromosome == "3L":
                if start >= 750000 and end <= 19026900:
                    print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                    print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
            if chromosome == "3R":
                if start >= 6754278 and end <= 31614278:
                    print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                    print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
            if chromosome == "4":
                print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
            if chromosome == "X":
                if start >= 1325967 and end <= 21338973:
                    print >> arg_final, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)
                    print >> final_file, chromosome + "\t" + str(start) + "\t" + str(end) + "\t" + strain + "\t" + str(score_final)