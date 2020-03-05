#!/usr/bin/python

#### By Maria Bogaerts
#### Selecting 30 out of the 48 total samples, and calculating Tajima's D mean among all the samples

import sys
import os

# Tajima's D initial file
tajimasd = open("DrosEU_50000_10000.D", 'r')

# Final file
mean_file = open("window_mean_30pops.txt", 'w')

for t in tajimasd:
    chromosome = t.split("\t")[0]
    position = float(t.split("\t")[1])
    populations = []
    for x in range(11,40): #(2,49)
        if t.split("\t")[x] == "NA":
            print "NA"
        else:
            populations.append(float(t.split("\t")[x]))

    print(populations)
    if len(populations) != 0:
        mean = (sum(populations))/len(populations)
    else:
        mean = "NA"
    print >> mean_file, chromosome + "\t" + str(int(position)-25000) + "\t" + str(int(position)+25000) + "\t" + str(mean)