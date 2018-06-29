# coding: utf-8
# # Script for calculating the occurence of double-crossover events in simulated African American data

import numpy as np
import re
import pprint

#load in an example A and B bed file for one individual, indiv 1
Atractsfile = open('/Users/elizabeth/Dropbox/Projects/MGH-Broad/PGC/PTSD/LAI/FixSwitchErrors/Grady/simualtedTruth/SimGrady1.1.A.bed', 'r')
Btractsfile = open('/Users/elizabeth/Dropbox/Projects/MGH-Broad/PGC/PTSD/LAI/FixSwitchErrors/Grady/simualtedTruth/SimGrady1.1.B.bed', 'r')
out = open('/Users/elizabeth/Dropbox/Projects/MGH-Broad/PGC/PTSD/LAI/FixSwitchErrors/Grady/DoubleRecombs_truthGrady.txt', 'w')
Atracts = np.genfromtxt(Atractsfile, dtype=None)
Btracts = np.genfromtxt(Btractsfile, dtype=None)
#convert the tracts files to np dataframes


Atracts[:2]


# ### want to see if there are double-crossover events within a given window size in the A and B 
#  #just bother with the start and end in cM, not bp

for line in range(0, len(Atracts)):
    chromA = Atracts[line][0]
    startA = Atracts[line][5]
    endA = Atracts[line][7]
    ancA = Atracts[line][3]
    
    for item in range(0, len(Btracts)):
        chromB = Btracts[item][0]
        startB = Btracts[item][5]
        endB = Btracts[item][7]
        ancB = Btracts[item][3]

        count=0

        if chromA == chromB:
            if ancA != ancB:
                #match to ensure that the ancestries are opposing. Right now, only modeling 2 way admixed scenario
                
                #if the start and stop positions are within 10 cM of each other, consider this to be the same location
                if (startA-10) < startB < (startA+10):
                    if (endA-10) < endB < (endA+10):
                        count = count +1
                    if count > 0:
                        output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chromA), str(startA), str(endA), str(ancA), str(startB), str(endB),str(ancB))
                        #output the chromosome number, start and end pos of A and B and their respective ancs for sanity check
                        #these are the double recombination events - the ones where there are matching tracts in opposing ancestries
                        out.write(output)
out.close()
